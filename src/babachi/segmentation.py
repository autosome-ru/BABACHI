from abc import ABC, abstractmethod
import numpy as np
import logging
import multiprocessing as mp
from scipy.special import logsumexp
import time
from typing import List

from babachi.stats import fast_find_optimal_borders
from babachi.logging import root_logger, set_logger_config
from babachi.models import BADSegmentsContainer, filter_segments, BADSegment, GenomeSNPsHandler
from babachi.helpers import pack


class Segmentation(ABC):
    def __init__(self):
        self.dtype = np.float64
        self.segments_container = BADSegmentsContainer()
        self.S = None  # snp-wise log-likelihoods for each BAD
        self.C = None  # snp-wise marginal log-likelihoods
        self.P = None  # segment-wise log-likelihoods for each BAD
        self.L = None  # segment-wise marginal log-likelihoods

        self.inf_score = 100000000  # Log-likelihood penalty for Z and R

        self.boundaries_indexes = None
        self.boundary_numbers = None
        self.snps_positions = []
        self.total_cover = None

        self.total_snps_count = None
        self.last_snp_number = None
        self.candidate_numbers = None
        self.candidates_count = None
        self.sub_chromosome: 'SubChromosomeSegmentation' = None
        self.score = None  # score[i] = best log-likelihood among all segmentations of snps[0,i]
        self.has_boundary_cache = None  # bool boundaries, len=total_snps_count * (total_snps_count + 1).
        self.best_boundaries_count = None  # best_boundaries_count[i] = number of boundaries
        # before ith snp in best segmentation

    @staticmethod
    def get_norm(p, N, trim_cover):
        result = 0
        current_multiplier = 1
        denominator_multiplier = 1
        for k in range(trim_cover):
            result += current_multiplier * np.power(p, N - k) * np.power(1 - p, k) / denominator_multiplier
            current_multiplier *= int(N - k)
            denominator_multiplier *= k + 1

        return -result


    def log_likelihood(self, N, X, BAD):
        """
        allele_reads_tr <= X <= N/2
        """
        p = 1.0 / (1.0 + BAD)
        #FIXME
        log_norm = np.log1p(
            self.get_norm(p, N, self.sub_chromosome.gs.allele_reads_tr) +
            self.get_norm(1 - p, N, self.sub_chromosome.gs.allele_reads_tr)
        )
        if (self.sub_chromosome.gs.individual_likelihood_mode in (
                'corrected',
                'bayesian') and N == 2 * X) or self.sub_chromosome.gs.individual_likelihood_mode == 'binomial':
            return X * np.log(p) + (N - X) * np.log(1 - p) + np.log(self.sub_chromosome.gs.prior[BAD]) - log_norm
        elif self.sub_chromosome.gs.individual_likelihood_mode == 'corrected':
            return X * np.log(p) + (N - X) * np.log(1 - p) + np.log(self.sub_chromosome.gs.prior[BAD]) - log_norm \
                   + np.log1p(float(BAD) ** (2 * X - N))
        elif self.sub_chromosome.gs.individual_likelihood_mode == 'bayesian':
            return X * np.log(p) + (N - X) * np.log(1 - p) + np.log(self.sub_chromosome.gs.prior[BAD]) - log_norm \
                   - np.log1p(float(BAD) ** (2 * X - N)) + np.log1p(float(BAD) ** (4 * X - 2 * N))

    def get_P(self, first, last):
        if last - first == 0:
            return self.sub_chromosome.P_initial[:, last]
        else:
            return np.sum(self.sub_chromosome.P_initial[:, first: last + 1], axis=1)

    def construct_likelihood_matrices(self):
        P = np.zeros((len(self.sub_chromosome.gs.BAD_list), self.candidates_count + 1,
                      self.candidates_count + 1),
                     dtype=self.dtype)
        S = np.zeros((len(self.sub_chromosome.gs.BAD_list), self.candidates_count + 1),
                     dtype=self.dtype)
        self.L = np.zeros((self.candidates_count + 1, self.candidates_count + 1),
                          dtype=self.dtype)  # if in init -> -memory
        for j in range(0, self.candidates_count + 1):
            if j == 0:
                first = 0
            else:
                first = self.candidate_numbers[j - 1] + 1
            if j == self.candidates_count:
                last = self.last_snp_number
            else:
                last = self.candidate_numbers[j]

            S[:, j] = self.get_P(first, last)
        self.S = S
        self.P = P

    def modify_P(self):
        self.C = np.cumsum(self.S, axis=1)
        B = self.C.transpose()
        self.P = (B[:, np.newaxis, :] - np.insert(B, 0, 0, axis=0)).transpose()[:, :-1, :]

    def modify_L(self):
        if self.sub_chromosome.gs.scoring_mode == 'marginal':
            self.L[:, :] = logsumexp(self.P, axis=0)
        elif self.sub_chromosome.gs.scoring_mode == 'maximum':
            self.L[:, :] = self.P.max(axis=0)

    def get_parameter_penalty(self, boundaries):
        if isinstance(self, AtomicRegionSegmentation):
            N = self.total_snps_count
        else:
            N = self.sub_chromosome.unique_snp_positions

        return -1 / 2 * boundaries * (np.log(N) + 1) * self.sub_chromosome.gs.b_penalty

    def initialize_boundaries_arrays(self):
        self.score = [0] * (self.candidates_count + 1)
        self.has_boundary_cache = [[False] * self.candidates_count for _ in range(self.candidates_count + 1)]
        self.best_boundaries_count = [0] * (self.candidates_count + 1)

    @abstractmethod
    def find_optimal_boundaries(self):
        pass

    def estimate(self):
        self.construct_likelihood_matrices()
        self.modify_P()
        self.modify_L()
        if self.sub_chromosome.gs.fast and isinstance(self, AtomicRegionSegmentation):
            self.boundaries_indexes = fast_find_optimal_borders(
                self.candidates_count,
                self.L,
                self.sub_chromosome.gs.min_seg_snps,
                self.sub_chromosome.gs.min_seg_bp,
                self.snps_positions,
                np.array(self.candidate_numbers, dtype=np.int_),
                self.first_snp_number,
                np.float64(self.inf_score),
                self.total_snps_count,
                self.sub_chromosome.gs.b_penalty
            )
        else:
            self.initialize_boundaries_arrays()
            self.find_optimal_boundaries()


class AtomicRegionSegmentation(Segmentation):
    def __init__(self, sub_chromosome: 'SubChromosomeSegmentation', start, end):
        super().__init__()
        self.sub_chromosome = sub_chromosome
        self.start_snp_index = start
        self.end_snp_index = end
        self.total_snps_count = end - start + 1
        if self.end_snp_index == sub_chromosome.candidates_count:
            self.last_snp_number = sub_chromosome.total_snps_count - 1
        else:
            self.last_snp_number = sub_chromosome.candidate_numbers[end]
        if self.start_snp_index == 0:
            self.first_snp_number = 0
        else:
            self.first_snp_number = sub_chromosome.candidate_numbers[start - 1] + 1
        self.snps_positions = sub_chromosome.snps_positions[
            self.first_snp_number: self.last_snp_number + 1
        ]
        self.total_cover = sub_chromosome.allele_read_counts_array[
            self.first_snp_number:
            self.last_snp_number + 1
        ].sum()
        self.candidate_numbers = sub_chromosome.candidate_numbers[start:end]
        self.candidates_count = end - start

    def find_optimal_boundaries(self):
        for i in range(self.candidates_count + 1):
            self.score[i] = self.L[0, i]

            kf = -1
            current_optimal = self.score[i]

            check_optimal = True
            if self.sub_chromosome.gs.min_seg_snps or self.sub_chromosome.gs.min_seg_bp:
                piece_positions = self.snps_positions[0: self.candidate_numbers[i] + 1 - self.first_snp_number] \
                    if i != self.candidates_count else self.snps_positions
                if (self.sub_chromosome.gs.min_seg_snps and len(
                        np.unique(piece_positions)) < self.sub_chromosome.gs.min_seg_snps) or \
                        (self.sub_chromosome.gs.min_seg_bp and piece_positions[-1] - piece_positions[
                            0] < self.sub_chromosome.gs.min_seg_bp):
                    self.score[i] -= self.inf_score
                    check_optimal = False
            if check_optimal:
                for k in range(i):
                    parameter_penalty = self.get_parameter_penalty(self.best_boundaries_count[k] + 1)

                    z_penalty = 0
                    if self.sub_chromosome.gs.min_seg_snps or self.sub_chromosome.gs.min_seg_bp:
                        piece_positions = self.snps_positions[
                                          self.candidate_numbers[k] + 1 - self.first_snp_number: self.candidate_numbers[i] + 1 - self.first_snp_number
                            ] \
                            if i != self.candidates_count else self.snps_positions[
                                                               self.candidate_numbers[k] + 1 - self.first_snp_number:]
                        if len(np.unique(piece_positions)) < self.sub_chromosome.gs.min_seg_snps or \
                                piece_positions[-1] - piece_positions[0] < self.sub_chromosome.gs.min_seg_bp:
                            z_penalty = -self.inf_score

                    likelihood = self.score[k] + self.L[k + 1, i]
                    candidate = likelihood + parameter_penalty + z_penalty
                    if candidate > current_optimal:
                        current_optimal = candidate
                        self.score[i] = likelihood + z_penalty
                        kf = k
            if kf != -1:
                for j in range(0, kf):
                    self.has_boundary_cache[i][j] = self.has_boundary_cache[kf][j]
                self.has_boundary_cache[i][kf] = 1
                for j in range(kf + 1, i):
                    self.has_boundary_cache[i][j] = 0

            self.best_boundaries_count[i] = sum(self.has_boundary_cache[i])

        self.boundaries_indexes = [self.candidate_numbers[j] for j in range(self.candidates_count) if
                                   self.has_boundary_cache[-1][j]]


class SubChromosomeSegmentation(Segmentation):  # sub_chromosome
    def __init__(
            self,
            genome_segmentator: 'GenomeSegmentator',
            chromosome_segmentation: 'ChromosomeSegmentation', 
            allele_read_counts_array,
            snps_positions,
            part
        ):
        super().__init__()

        self.gs = genome_segmentator
        self.chromosome_segmentation = chromosome_segmentation
        self.sub_chromosome = self
        self.index_in_chromosome = part
        self.allele_read_counts_array = allele_read_counts_array

        self.total_snps_count = len(self.allele_read_counts_array)
        self.total_cover = self.allele_read_counts_array.sum()
        self.snps_positions = snps_positions
        assert len(self.snps_positions) == self.total_snps_count
        self.unique_snp_positions = len(np.unique(self.snps_positions))
        self.end_snp_index = (self.total_snps_count - 1) - 1  # index from 0, and #boundaries = #snps - 1
        self.candidate_numbers = [i for i in range(self.total_snps_count - 1) if
                                  self.snps_positions[i] != self.snps_positions[i + 1]]
        self.candidates_count = len(self.candidate_numbers)
        self.last_snp_number = self.total_snps_count - 1

        self.P_initial = None  # snp-wise log-likelihoods for each BAD

    @staticmethod
    def split_into_overlapping_regions(total_length, length_of_region, overlap):
        result = []
        if total_length <= length_of_region:
            result.append((0, total_length - 1))
            return result
        total_length -= overlap
        div, _ = divmod(total_length - 1, length_of_region - overlap)
        new_l, num = divmod(total_length - 1, div)
        # --- ----- ------ ----- never mind, just my colleague code ---- -------
        for i in range(div):
            if i < num:
                result.append(((new_l + 1) * i, (new_l + 1) * (i + 1) + overlap))
            else:
                result.append((new_l * i + num, new_l * (i + 1) + overlap + num))
        return result

    def set_candidates(self, candidate_set):
        self.candidate_numbers = sorted(list(candidate_set))
        self.candidates_count = len(self.candidate_numbers)
        self.end_snp_index = self.candidates_count - 1

    def construct_initial_likelihood_matrices(self):
        vector_likelihood = np.vectorize(self.log_likelihood, excluded=['BAD'])
        S = np.zeros((len(self.gs.BAD_list), self.total_snps_count), dtype=self.dtype)
        if self.sub_chromosome.gs.individual_likelihood_mode == 'binomial':
            X = self.allele_read_counts_array[:, 0]
        else:
            X = self.allele_read_counts_array.min(axis=1)

        N = self.allele_read_counts_array.sum(axis=1)
        for i in range(len(self.gs.BAD_list)):
            S[i, :] = vector_likelihood(N, X, BAD=self.gs.BAD_list[i])
        self.P_initial = S

    def find_optimal_boundaries(self):
        for i in range(self.candidates_count + 1):
            self.score[i] = self.L[0, i]

            last_boundary_index = -1
            current_optimal = self.score[i]

            check_optimal = True
            if self.sub_chromosome.gs.min_seg_snps or self.sub_chromosome.gs.min_seg_bp:
                piece_positions = self.snps_positions[0: self.candidate_numbers[i] + 1] \
                    if i != self.candidates_count else self.snps_positions
                if (self.sub_chromosome.gs.min_seg_snps and len(
                        np.unique(piece_positions)) < self.sub_chromosome.gs.min_seg_snps) or \
                        (self.sub_chromosome.gs.min_seg_bp and piece_positions[-1] - piece_positions[
                            0] < self.sub_chromosome.gs.min_seg_bp):
                    self.score[i] -= self.inf_score
                    check_optimal = False
            if check_optimal:
                for k in range(i):
                    parameter_penalty = self.get_parameter_penalty(self.best_boundaries_count[k] + 1)

                    z_penalty = 0
                    if self.gs.min_seg_snps or self.gs.min_seg_bp:
                        piece_positions = self.snps_positions[
                                          self.candidate_numbers[k] + 1: self.candidate_numbers[i] + 1] \
                            if i != self.candidates_count else self.snps_positions[self.candidate_numbers[k] + 1:]
                        if len(np.unique(piece_positions)) < self.sub_chromosome.gs.min_seg_snps or \
                                piece_positions[-1] - piece_positions[0] < self.sub_chromosome.gs.min_seg_bp:
                            z_penalty = -self.inf_score

                    likelihood = self.score[k] + self.L[k + 1, i]
                    candidate = likelihood + parameter_penalty + z_penalty
                    if candidate > current_optimal:
                        current_optimal = candidate
                        self.score[i] = likelihood + z_penalty
                        last_boundary_index = k
            if last_boundary_index != -1:
                for j in range(0, last_boundary_index):
                    self.has_boundary_cache[i][j] = self.has_boundary_cache[last_boundary_index][j]
                self.has_boundary_cache[i][last_boundary_index] = 1
                for j in range(last_boundary_index + 1, i):
                    self.has_boundary_cache[i][j] = 0

            self.best_boundaries_count[i] = sum(self.has_boundary_cache[i])

        self.boundaries_indexes = [self.candidate_numbers[j] for j in range(self.candidates_count) if
                                   self.has_boundary_cache[-1][j]]
        self.boundary_numbers = [-1] + [j for j in range(self.candidates_count) if self.has_boundary_cache[-1][j]] + [
            self.candidates_count]
        cumulative_counts = [0] + [x + 1 for x in self.boundaries_indexes] + [self.last_snp_number + 1]
        self.segments_container.snps_counts = [cumulative_counts[i + 1] - cumulative_counts[i] for i in
                                               range(len(cumulative_counts) - 1)]

        self.segments_container.snp_id_counts = [
            len(set(self.snps_positions[cumulative_counts[i]: cumulative_counts[i + 1]])) for i in
            range(len(cumulative_counts) - 1)
        ]
        self.segments_container.covers = [
            sum(ref_count + alt_count for ref_count, alt_count in
                self.allele_read_counts_array[cumulative_counts[i]:cumulative_counts[i + 1]]) for i in
            range(len(cumulative_counts) - 1)]

        for i in range(len(self.has_boundary_cache[-1])):
            if self.has_boundary_cache[-1][i]:
                self.segments_container.boundaries_positions.append(
                    (self.snps_positions[self.candidate_numbers[i]] + self.snps_positions[
                        self.candidate_numbers[i] + 1]) / 2)

    def estimate_BAD(self):
        for n in range(len(self.boundary_numbers) - 1):
            first = self.boundary_numbers[n] + 1
            last = self.boundary_numbers[n + 1]
            likelihoods = self.P[:, first, last] - self.L[first, last]
            estimated_BAD_index = np.argmax(likelihoods)
            self.segments_container.BAD_estimations.append(self.gs.BAD_list[estimated_BAD_index])
            self.segments_container.likelihoods.append(list(likelihoods))

    def estimate_sub_chr(self):
        if self.total_snps_count == 0:
            return

        self.gs.logger.debug('Constructing initial SNP-wise likelihoods...')
        self.construct_initial_likelihood_matrices()
        atomic_regions_limits = self.split_into_overlapping_regions(
            self.candidates_count + 1,
            self.gs.atomic_region_length,
            self.gs.overlap
        )
        boundary_set = set()
        counter = 0
        self.gs.logger.debug('Segmenting atomic regions:')
        for first, last in atomic_regions_limits:
            counter += 1
            self.gs.logger.debug(
                'Processing {}/{} atomic regions from SNP {} to {} for {} (subchromosome {} of {}).'.format(
                    counter, len(atomic_regions_limits),
                    first, last, self.chromosome_segmentation.chromosome,
                    self.index_in_chromosome,
                    len(self.chromosome_segmentation.get_sub_chromosomes_slices())))
            atomic_region_segmentation = AtomicRegionSegmentation(self, first, last)
            atomic_region_segmentation.estimate()
            boundary_set |= set(atomic_region_segmentation.boundaries_indexes)
        self.candidate_numbers = sorted(list(boundary_set))
        self.candidates_count = len(self.candidate_numbers)
        self.gs.logger.debug('Unique SNPs positions in subchromosome {}: {}'.format(self.index_in_chromosome,
                                                                                    self.unique_snp_positions))

        self.estimate()
        self.estimate_BAD()


class ChromosomeSegmentation:  # chromosome
    def __init__(self, genome_segmentator: 'GenomeSegmentator', chromosome, length):
        self.gs = genome_segmentator

        self.chromosome = chromosome  # name
        self.length = length  # length, bp

        self.allele_read_counts_array = self.gs.snps_collection[self.chromosome].read_counts
        self.snps_positions = self.gs.snps_collection[self.chromosome].positions
        self.total_snps_count = len(self.allele_read_counts_array)
        self.segments_container = BADSegmentsContainer()
        if self.total_snps_count == 0:
            return
        self.total_read_coverage = self.allele_read_counts_array.sum()
        self.critical_gap_factor = 1 - 10 ** (- 1 / np.sqrt(self.total_snps_count))
        self.effective_length = self.snps_positions[-1] - self.snps_positions[0]


    def adjust_critical_gap(self):
        condition = True
        length_difference = 0
        black_list_i = set()
        while condition:
            self.effective_length -= length_difference
            length_difference = 0

            for i in range(self.total_snps_count - 1):
                if i in black_list_i:
                    continue
                difference = self.snps_positions[i + 1] - self.snps_positions[i]
                if difference > (self.effective_length - difference) * self.critical_gap_factor:
                    length_difference += difference
                    black_list_i.add(i)

            condition = length_difference != 0

    def get_sub_chromosomes_slices(self):
        sub_chromosome_slice_indexes = []
        current_tuple_start = 0
        for i in range(self.total_snps_count - 1):
            if self.snps_positions[i + 1] - self.snps_positions[i] > self.critical_gap_factor * self.effective_length:
                sub_chromosome_slice_indexes.append((current_tuple_start, i + 1))
                current_tuple_start = i + 1
        sub_chromosome_slice_indexes.append((current_tuple_start, self.total_snps_count))
        return sub_chromosome_slice_indexes

    def estimate_chr(self) -> BADSegmentsContainer:
        self.gs.logger.info('Processing SNPs in {}'.format(self.chromosome))
        if not self.total_snps_count or self.total_snps_count < self.gs.snp_per_chr_tr:
            return self.segments_container

        start_t = time.perf_counter()
        self.adjust_critical_gap()

        #  boundary for first snp
        if self.snps_positions[0] <= self.critical_gap_factor * self.effective_length:
            self.segments_container.boundaries_positions.append(self.snps_positions[0])
        else:
            self.segments_container.boundaries_positions.append((1, self.snps_positions[0]))
        self.gs.logger.debug(
            'Stage 1 subchromosomes (start SNP index, end SNP index): {}'.format(
                self.get_sub_chromosomes_slices()))

        for part, (start, end) in enumerate(self.get_sub_chromosomes_slices(), 1):
            # check
            unique_positions = len(np.unique(self.snps_positions[start: end]))
            if unique_positions < self.gs.min_subchr_length:
                self.segments_container += BADSegmentsContainer(
                    boundaries_positions=[],
                    BAD_estimations=[0],
                    likelihoods=[[0] * len(self.gs.BAD_list)],
                    snps_counts=[end - start],
                    covers=[0],
                    snp_id_counts=[unique_positions]
                )
            else:
                sub_chromosome = SubChromosomeSegmentation(
                    self.gs,
                    self, 
                    self.allele_read_counts_array[start: end],
                    self.snps_positions[start: end],
                    part
                )
                start_t = time.perf_counter()
                sub_chromosome.estimate_sub_chr()
                self.gs.logger.debug('Subchromosome time: {}, subchromosome SNPs: {}'.format(
                    time.perf_counter() - start_t, unique_positions
                ))

                self.segments_container += sub_chromosome.segments_container
            if end != self.total_snps_count:
                self.segments_container.boundaries_positions += [(self.snps_positions[end - 1], self.snps_positions[end])]

        #  boundary for last snp
        if self.length - self.snps_positions[-1] <= self.critical_gap_factor * self.effective_length:
            self.segments_container.boundaries_positions.append(self.length)
        else:
            self.segments_container.boundaries_positions.append((self.snps_positions[-1] + 1, self.length))

        boundaries_print = [
            '{:.2f}Mbp'.format(x / 1000000)
            if np.issubdtype(type(x), np.number) else
            '{:.2f}Mbp[{:.2f}Mbp]'.format(x[0] / 1000000, (x[1] - x[0]) / 1000000)
            for x in self.segments_container.boundaries_positions
        ]
        self.gs.logger.debug(
            '\nEstimated BADs: {}\nSNP counts: {}\nSNP IDs counts: {}\nCritical gap: {:.0f}bp'
            '\nBoundaries positions (location[deletion length (if any)]): {}'.format(
                '[' + ', '.join('{:.2f}'.format(BAD) for BAD in self.segments_container.BAD_estimations) + ']',
                self.segments_container.snps_counts, self.segments_container.snp_id_counts,
                self.critical_gap_factor * self.effective_length,
                '[' + ', '.join(boundaries_print) + ']')
        )
        self.gs.logger.debug('{} time: {} s'.format(self.chromosome, time.perf_counter() - start_t))
        return self.segments_container


class GenomeSegmentator:  # gs
    def __init__( # TODO move to config
        self, 
        snps_collection: GenomeSNPsHandler,
        chrom_sizes: dict,
        segmentation_mode='corrected', 
        scoring_mode='marginal',
        states=None, b_penalty=4, prior=None,
        allele_reads_tr=5, min_seg_snps=3, min_seg_bp=0,
        post_seg_filter=0,
        jobs=1,
        atomic_region_size=600, chr_filter=100, subchr_filter=3,
        logger=root_logger,
        logger_level=logging.INFO,
    ):
        self.logger = logger
        self.logger_level = logger_level
        self.individual_likelihood_mode = segmentation_mode  # 'corrected', 'binomial' or 'bayesian'
        self.scoring_mode = scoring_mode  # marginal or maximum
        self.b_penalty = b_penalty  # boundary penalty coefficient k ('CAIC' * k)
        self.allele_reads_tr = allele_reads_tr  # "minimal read count on each allele" snp filter
        self.post_seg_filter = post_seg_filter  # min number of SNPs in the segment to be included in the output
        self.jobs = jobs

        if states is None or len(states) == 0:
            self.BAD_list = [1, 2, 3, 4, 5, 6]
        else:
            self.BAD_list = sorted(states)
        if prior is None:
            prior = dict(zip(self.BAD_list, [1] * len(self.BAD_list)))
        else:
            if len(prior) != len(self.BAD_list):
                raise AssertionError('Length of prior = {} is not equal to number of states = {}'.format(
                    len(prior), len(self.BAD_list)))
        self.prior = prior

        self.snps_collection = snps_collection  # input snp collection

        self.snp_per_chr_tr = chr_filter  # minimal number of snps in chromosome to start segmentation
        self.atomic_region_length = atomic_region_size  # length of an atomic region in snps
        self.overlap = atomic_region_size // 2  # length of regions overlap in snps
        self.min_subchr_length = subchr_filter  # minimal subchromosome length in snps
        self.fast = True  # use numba to optimize execution speed
        self.min_seg_snps = min_seg_snps  # minimal BAD segment length in SNPs
        self.min_seg_bp = min_seg_bp  # minimal BAD segment length in bp
        self.chrom_sizes = chrom_sizes  # ['chr1', 'chr2', ...]

        self.chr_segmentations = []  # list of ChromosomeSegmentation instances

        for chromosome, length in self.chrom_sizes.items():
            chr_segmentation = ChromosomeSegmentation(self, chromosome, length)
            self.logger.debug('{} total SNP count: {}'.format(chromosome, chr_segmentation.total_snps_count))
            self.logger.debug('{} length: {}'.format(chromosome, length))
            self.chr_segmentations.append(chr_segmentation)
        else:
            self.logger.debug('-----------------------------------------')

    def __getstate__(self):
        state = self.__dict__
        if 'logger' in state:
            del state['logger']
        return state

    def __setstate__(self, state):
        for name, value in state.items():
            setattr(self, name, value)
        assert not hasattr(self, 'logger')
        self.logger = root_logger
        set_logger_config(self.logger, self.logger_level)

    def start_chromosome(self, segmentation: ChromosomeSegmentation):
        return segmentation.estimate_chr()
    
    def estimate_BAD(self) -> List[BADSegment]:
        results = []
        jobs = min(self.jobs, len(self.chr_segmentations), max(1, mp.cpu_count()))

        if jobs == 0:
            return results
        
        elif jobs == 1:
            # Run serially
            for i, segmentation in enumerate(self.chr_segmentations):
                segments_container = self.start_chromosome(segmentation)
                # awful coding practice, need to rewrite at some point
                res = segments_container.get_BAD_segments(segmentation)
                results.extend(filter_segments(res, self.post_seg_filter))
                self.chr_segmentations[i] = None
        else:
            # Run in parallel using multiprocessing
            ctx = mp.get_context("forkserver")
            with ctx.Pool(jobs) as pool:
                for i, res in zip(
                    np.arange(len(self.chr_segmentations)),
                    pool.map(self.start_chromosome, self.chr_segmentations)):
                    res = res.get_BAD_segments(self.chr_segmentations[i])
                    results.extend(filter_segments(res, self.post_seg_filter))
                    self.chr_segmentations[i] = None
        return results

    def write_BAD(self, segments, outpath):
        with open(outpath, 'w') as outfile:
            outfile.write(
                pack(['#chr', 'start', 'end', 'BAD', 'SNP_count', 'SNP_ID_count', 'sum_cover'] + ['Q{:.2f}'.format(BAD) for BAD in self.BAD_list]))
            for segment in segments:
                outfile.write(str(segment))


