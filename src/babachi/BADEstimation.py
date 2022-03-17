"""
Usage:
    babachi (<file> | --test) [options]
    babachi visualize <file> (-b <badmap>| --badmap <badmap>)

Arguments:
    <file>            Path to input file in tsv format with columns:
                      chr pos ID ref_base alt_base ref_read_count alt_read_count
                      Expected to be sorted be (chr, pos)
    <badmap>          Path to badmap .bed format file
    <int>             Non negative integer
    <float>           Non negative number
    <string>          String of states separated with "," (to provide fraction use "/", e.g. 4/3). Each state must be >= 1


Required arguments:
    --test                                  Run segmentation on test file
    -b <badmap>, --badmap <badmap>          Input badmap file
    -O <path>, --output <path>              Output directory or file path. [default: ./]

Optional arguments:
    -h, --help                              Show help.
    -q, --quiet                             Suppress log messages.
    --force-sort                            Chromosomes will be sorted in numerical order

    --allele-reads-tr <int>                 Allelic reads threshold. Input SNPs will be filtered by ref_read_count >= x and
                                            alt_read_count >= x. [default: 5]
    -B <float>, --boundary-penalty <float>  Boundary penalty coefficient [default: 4]
    --states <string>                       States string [default: 1,2,3,4,5,6]
    -Z <int>, --min-seg-snps <int>          Only allow segments containing Z or more unique SNPs (IDs/positions) [default: 3]
    -R <int>, --min-seg-bp <int>            Only allow segments containing R or more base pairs [default: 1000]
    -P <int>, --post-segment-filter <int>   Remove segments with less than P unique SNPs (IDs/positions) from output [default: 0]
    -A <int>, --atomic-region-size <int>    Atomic region size in SNPs [default: 600]
    -C <int>, --chr-min-snps <int>          Minimum number of SNPs on a chromosome to start segmentation [default: 100]
    -S <int>, --subchr-filter <int>         Exclude subchromosomes with less than C unique SNPs  [default: 3]

Visualization:
    --visualize                             Perform visualization of SNP-wise AD and BAD for each chromosome.
                                            Will create a directory in output path for the .svg visualizations.
"""

import math
import re

import numpy as np
import vcf
from numba import njit
import os.path

from schema import Schema, And, Use, SchemaError, Const, Or
import time
from .helpers import ChromosomePosition, pack
from abc import ABC, abstractmethod
from docopt import docopt
from .visualize_segmentation import init_from_snps_collection

# TODO check if input sorted


class BADSegmentsContainer:
    allowed_fields = ['boundaries_positions', 'BAD_estimations', 'likelihoods', 'snps_counts', 'covers',
                      'snp_id_counts']

    def __init__(self, **kwargs):
        self.BAD_estimations = []  # estimated BADs for split segments
        self.likelihoods = []  # likelihoods of split segments for each BAD
        self.snps_counts = []  # number of snps in segments
        self.snp_id_counts = []  # numer of unique snp positions in segments
        self.covers = []  # sums of covers for each segment
        self.boundaries_positions = []  # boundary positions between snps at x1 and x2:
        # (x1+x2)/2 for x2-x1<=CRITICAL_GAP, (x1, x2) else
        for arg in self.allowed_fields:
            if arg in kwargs:
                assert isinstance(kwargs[arg], list)
                setattr(self, arg, kwargs[arg])

    def __add__(self, other):
        if not isinstance(other, BADSegmentsContainer):
            raise NotImplementedError
        return BADSegmentsContainer(**{arg: getattr(self, arg) + getattr(other, arg) for arg in self.allowed_fields})

    def __iadd__(self, other):
        if not isinstance(other, BADSegmentsContainer):
            raise NotImplementedError
        for arg in self.allowed_fields:
            setattr(self, arg, getattr(self, arg) + getattr(other, arg))
        return self

    def get_BAD_segments(self, chromosome_segmentation):
        current_position = None
        if chromosome_segmentation.total_snps_count >= chromosome_segmentation.gs.snp_per_chr_tr:
            for counter, boundary in enumerate(self.boundaries_positions, -1):
                if current_position is None:
                    if isinstance(boundary, tuple):
                        current_position = boundary[1]
                    else:
                        current_position = 1
                elif isinstance(boundary, tuple):
                    yield BADSegment(
                        chromosome_segmentation.chromosome,
                        current_position,
                        boundary[0] + 1,
                        self.BAD_estimations[counter],
                        self.likelihoods[counter],
                        self.snps_counts[counter],
                        self.covers[counter],
                        self.snp_id_counts[counter],
                    )
                    current_position = boundary[0] + 1
                    yield BADSegment(
                        chromosome_segmentation.chromosome,
                        current_position,
                        boundary[1],
                        0,
                        [0] * len(chromosome_segmentation.gs.BAD_list),
                        0,
                        0,
                        0,
                    )
                    current_position = boundary[1]
                else:
                    yield BADSegment(
                        chromosome_segmentation.chromosome,
                        current_position,
                        math.floor(boundary) + 1,
                        self.BAD_estimations[counter],
                        self.likelihoods[counter],
                        self.snps_counts[counter],
                        self.covers[counter],
                        self.snp_id_counts[counter],
                    )
                    current_position = math.floor(boundary) + 1


class BADSegment:
    def __init__(self, chr, start, end, BAD, likelihoods, snps_count, total_cover, snp_id_count):
        self.chr = chr
        self.start = start
        self.end = end
        self.BAD = BAD
        self.likelihoods = likelihoods
        self.snps_count = snps_count
        self.total_cover = total_cover
        self.snp_id_count = snp_id_count

    def __repr__(self):
        return pack([self.chr, self.start, self.end, self.BAD, self.snps_count, self.snp_id_count, self.total_cover,
                     *self.likelihoods])


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
        self.sub_chromosome = None
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
        log_norm = np.log1p(self.get_norm(p, N, self.sub_chromosome.gs.allele_reads_tr) +
                            self.get_norm(1 - p, N, self.sub_chromosome.gs.allele_reads_tr))
        if (self.sub_chromosome.gs.individual_likelihood_mode in ('corrected', 'bayesian') and N == 2 * X) or self.sub_chromosome.gs.individual_likelihood_mode == 'binomial':
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
            Q = np.sort(self.P, axis=0)
            self.L[:, :] = Q[-1, :, :] + np.log1p(np.sum(np.exp(Q[:-2, :, :] - Q[-1, :, :]), axis=0))
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
        self.has_boundary_cache = [[False] * self.candidates_count for i in range(self.candidates_count + 1)]
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
    def __init__(self, sub_chromosome, start, end):
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
        self.snps_positions = sub_chromosome.snps_positions[self.first_snp_number: self.last_snp_number + 1]
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
                                          self.candidate_numbers[k] + 1 - self.first_snp_number: self.candidate_numbers[
                                                                                                     i] + 1 - self.first_snp_number] \
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
    def __init__(self, genome_segmentator, chromosome_segmentation, allele_read_counts_array, snps_positions, part):
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

        if self.gs.verbose:
            print('Constructing initial SNP-wise likelihoods...')
        self.construct_initial_likelihood_matrices()
        atomic_regions_limits = self.split_into_overlapping_regions(self.candidates_count + 1,
                                                                    self.gs.atomic_region_length,
                                                                    self.gs.overlap)
        boundary_set = set()
        counter = 0
        if self.gs.verbose:
            print('Segmenting atomic regions:')
        for first, last in atomic_regions_limits:
            counter += 1
            if self.gs.verbose:
                print(
                    'Making {} out of {} atomic region{} from SNP {} to {} for {} (subchromosome {} of {}).'.format(
                        counter, len(atomic_regions_limits), 's' * bool((len(atomic_regions_limits) - 1)),
                        first, last, self.chromosome_segmentation.chromosome,
                        self.index_in_chromosome,
                        len(self.chromosome_segmentation.get_sub_chromosomes_slices())))
            atomic_region_segmentation = AtomicRegionSegmentation(self, first, last)
            atomic_region_segmentation.estimate()
            boundary_set |= set(atomic_region_segmentation.boundaries_indexes)
        self.candidate_numbers = sorted(list(boundary_set))
        self.candidates_count = len(self.candidate_numbers)
        if self.gs.verbose:
            print('Unique SNPs positions in subchromosome {}: {}'.format(self.index_in_chromosome, self.unique_snp_positions))

        self.estimate()
        self.estimate_BAD()


class ChromosomeSegmentation:  # chromosome
    def __init__(self, genome_segmentator, chromosome, length=0):
        self.gs = genome_segmentator

        self.chromosome = chromosome  # name
        self.length = length  # length, bp

        self.allele_read_counts_array, self.snps_positions = self.unpack_snp_collection()  # unpack
        self.total_snps_count = len(self.allele_read_counts_array)
        self.segments_container = BADSegmentsContainer()
        if self.total_snps_count == 0:
            return
        self.total_read_coverage = self.allele_read_counts_array.sum()
        self.critical_gap_factor = 1 - 10 ** (- 1 / np.sqrt(self.total_snps_count))
        self.effective_length = self.snps_positions[-1] - self.snps_positions[0]

    def unpack_snp_collection(self):  # TODO: redundant reshaping
        try:
            positions, snps = zip(
                *((pos, (ref_count, alt_count)) for pos, ref_count, alt_count in
                  self.gs.snps_collection[self.chromosome])
            )
        except ValueError:
            positions, snps = [], []
        return np.array(snps, dtype=np.int_), np.array(positions, dtype=np.int_)

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

    def estimate_chr(self):
        print('Processing SNPs in {}'.format(self.chromosome))
        if not self.total_snps_count or self.total_snps_count < self.gs.snp_per_chr_tr:
            return

        start_t = time.perf_counter()
        self.adjust_critical_gap()

        #  boundary for first snp
        if self.snps_positions[0] <= self.critical_gap_factor * self.effective_length:
            self.segments_container.boundaries_positions.append(self.snps_positions[0])
        else:
            self.segments_container.boundaries_positions.append((1, self.snps_positions[0]))
        if self.gs.verbose:
            print(
                'Stage 1 subchromosomes (start SNP index, end SNP index): {}'.format(
                    self.get_sub_chromosomes_slices()))

        for part, (st, ed) in enumerate(self.get_sub_chromosomes_slices(), 1):
            # check
            unique_positions = len(np.unique(self.snps_positions[st: ed]))
            if unique_positions < self.gs.min_subchr_length:
                self.segments_container += BADSegmentsContainer(
                    boundaries_positions=[],
                    BAD_estimations=[0],
                    likelihoods=[[0] * len(self.gs.BAD_list)],
                    snps_counts=[ed - st],
                    covers=[0],
                    snp_id_counts=[unique_positions]
                )
            else:
                sub_chromosome = SubChromosomeSegmentation(self.gs, self, self.allele_read_counts_array[st: ed],
                                                           self.snps_positions[st: ed], part)
                start_t = time.perf_counter()
                sub_chromosome.estimate_sub_chr()
                if self.gs.verbose:
                    print('Subchromosome time: {}, subchromosome SNPs: {}'.format(
                        time.perf_counter() - start_t, unique_positions
                    ))

                self.segments_container += sub_chromosome.segments_container
            if ed != self.total_snps_count:
                self.segments_container.boundaries_positions += [(self.snps_positions[ed - 1], self.snps_positions[ed])]

        #  boundary for last snp
        if self.length - self.snps_positions[-1] <= self.critical_gap_factor * self.effective_length:
            self.segments_container.boundaries_positions.append(self.length)
        else:
            self.segments_container.boundaries_positions.append((self.snps_positions[-1] + 1, self.length))

        if self.gs.verbose:
            print(
                '\nEstimated BADs: {}\nSNP counts: {}\nSNP IDs counts: {}\nCritical gap: {:.0f}bp'
                '\nBoundaries positions (location[deletion length (if any)]): {}'.format(
                    '[' + ', '.join('{:.2f}'.format(BAD) for BAD in self.segments_container.BAD_estimations) + ']',
                    self.segments_container.snps_counts, self.segments_container.snp_id_counts,
                    self.critical_gap_factor * self.effective_length,
                    '[' + ', '.join(map(
                        lambda x: '{:.2f}Mbp'.format(x/1000000) if isinstance(x, (int, float, np.int_, np.float_)) else (
                            '{:.2f}Mbp[{:.2f}Mbp]'.format(x[0]/1000000, (x[1] - x[0])/1000000)),
                        self.segments_container.boundaries_positions)) + ']'))
            print('{} time: {} s\n\n'.format(self.chromosome, time.perf_counter() - start_t))


class GenomeSegmentator:  # gs
    def __init__(self, snps_collection, out, chromosomes_order, segmentation_mode='corrected', scoring_mode='marginal', states=None,
                 b_penalty=4, prior=None, verbose=False, allele_reads_tr=5, min_seg_snps=3, min_seg_bp=0, post_seg_filter=0,
                 atomic_region_size=600, chr_filter=100, subchr_filter=3):

        self.verbose = verbose
        self.individual_likelihood_mode = segmentation_mode  # 'corrected', 'binomial' or 'bayesian'
        self.scoring_mode = scoring_mode  # marginal or maximum
        self.b_penalty = b_penalty  # boundary penalty coefficient k ('CAIC' * k)
        self.allele_reads_tr = allele_reads_tr  # "minimal read count on each allele" snp filter
        self.post_seg_filter = post_seg_filter  # min number of SNPs in the segment to be included in the output

        if states is None or len(states) == 0:
            self.BAD_list = [1, 2, 3, 4, 5]
        else:
            self.BAD_list = sorted(states)
        if prior is None:
            self.prior = dict(zip(self.BAD_list, [1] * len(self.BAD_list)))
        else:
            if len(prior) != len(self.BAD_list):
                raise AssertionError('Length of prior = {} is not equal to number of states = {}'.format(
                    len(prior), len(self.BAD_list)))
            self.prior = prior

        self.snps_collection = snps_collection  # input snp collection
        self.out = out  # output file in .bed format

        self.snp_per_chr_tr = chr_filter  # minimal number of snps in chromosome to start segmentation
        self.atomic_region_length = atomic_region_size  # length of an atomic region in snps
        self.overlap = atomic_region_size // 2  # length of regions overlap in snps
        self.min_subchr_length = subchr_filter  # minimal subchromosome length in snps
        self.fast = True  # use numba to optimize execution speed
        self.min_seg_snps = min_seg_snps  # minimal BAD segment length in SNPs
        self.min_seg_bp = min_seg_bp  # minimal BAD segment length in bp
        self.chromosomes_order = chromosomes_order  # ['chr1', 'chr2', ...]

        self.chr_segmentations = []  # list of ChromosomeSegmentation instances

        for chromosome in self.chromosomes_order:
            chr_segmentation = ChromosomeSegmentation(self, chromosome, ChromosomePosition.chromosomes[chromosome])
            if self.verbose:
                print('{} total SNP count: {}'.format(chromosome, chr_segmentation.total_snps_count))
            self.chr_segmentations.append(chr_segmentation)
        else:
            if self.verbose:
                print('\n')

    # noinspection PyTypeChecker
    def estimate_BAD(self):
        with open(self.out, 'w') as outfile:
            outfile.write(
                pack(['#chr', 'start', 'end', 'BAD', 'SNP_count', 'SNP_ID_count', 'sum_cover'] + ['Q{:.2f}'.format(BAD)
                                                                                                  for BAD in
                                                                                                  self.BAD_list]))
            for j in range(len(self.chr_segmentations)):
                chromosome = self.chr_segmentations[j]
                chromosome.estimate_chr()
                self.write_BAD_to_file(chromosome, outfile)
                self.chr_segmentations[j] = None

    def write_BAD_to_file(self, chromosome_segmentation, outfile):
        segments_generator = chromosome_segmentation.segments_container.get_BAD_segments(chromosome_segmentation)
        for segment in self.filter_segments(segments_generator):
            outfile.write(str(segment))

    def filter_segments(self, segments):
        for segment in segments:
            if not self.post_seg_filter:
                if segment.BAD != 0:
                    yield segment
            else:
                if segment.BAD != 0 and segment.snp_id_count >= self.post_seg_filter:
                    yield segment


def parse_input_file(opened_file, allele_reads_tr=5, force_sort=False):
    snps_collection = {chromosome: [] for chromosome in ChromosomePosition.chromosomes}
    chromosomes_order = []
    vcfReader = vcf.Reader(opened_file)
    if len(vcfReader.samples) != 1:
        return False
    for line_number, record in enumerate(vcfReader, 1):
        if len(record.samples) != 1:
            return False
        else:
            sample = record.samples[0]
        if record.CHROM not in ChromosomePosition.chromosomes:
            print('Invalid chromosome name: {} in line #{}'.format(record.CHROM, line_number))
            return False
        if len(record.ALT) > 1 or record.INFO['MAF'] < 0.05:
            continue
        if sample.data[0] != '0/1':
            continue
        ref_read_count, alt_read_count = sample.data[3]
        if min(ref_read_count, alt_read_count) < 5:
            continue
        if record.CHROM not in chromosomes_order:
            chromosomes_order.append(record.CHROM)
        if ref_read_count < allele_reads_tr or alt_read_count < allele_reads_tr:
            continue

        snps_collection[record.CHROM].append((record.start, ref_read_count, alt_read_count))
    if force_sort:
        chromosomes_order = ChromosomePosition.sorted_chromosomes
    return snps_collection, chromosomes_order, opened_file.name


def convert_frac_to_float(string):
    if re.match(r"^[1-9]+[0-9]*/[1-9]+[0-9]*$", string):
        num, denom = string.split('/')
        if int(denom) <= 0:
            return False
        else:
            value = int(num) / int(denom)
    elif re.match(r"^[1-9]+[0-9]*\.[1-9]+[0-9]*$", string):
        try:
            value = float(string)
        except ValueError:
            return False
    elif re.match(r"^[1-9]+[0-9]*$", string):
        try:
            value = int(string)
        except ValueError:
            return False
    else:
        return False
    if value >= 1:
        return value
    else:
        return False


def check_states(string):
    if not string:
        return False
    string = string.strip().split(',')
    ret_val = list(map(convert_frac_to_float, string))
    if not all(ret_val):
        return False
    else:
        return ret_val


@njit(cache=True)
def fast_find_optimal_borders(
        candidates_count,
        L,
        min_seg_snps,
        min_seg_bp,
        snps_positions,
        candidate_numbers,
        first_snp_number,
        inf_score,
        total_snps_count,
        b_penalty
):
    score = np.zeros(candidates_count + 1, dtype=np.float64)
    best_boundaries_count = np.zeros(candidates_count + 1, dtype=np.int_)
    has_boundary_cache = np.zeros((candidates_count + 1, candidates_count), dtype=np.bool_)
    if min_seg_snps:
        unique_positions = np.zeros(candidates_count + 1)
        current_index = 0
        for i in range(1, candidates_count + 1):
            if snps_positions[i] != snps_positions[i - 1]:
                current_index += 1
            unique_positions[i] = current_index
    for i in range(candidates_count + 1):
        score[i] = L[0, i]

        kf = -1
        current_optimal = score[i]

        check_optimal = True
        if min_seg_snps or min_seg_bp:
            last_index = candidate_numbers[i] + 1 - first_snp_number if i != candidates_count else -1
            first_index = 0
            if (min_seg_snps and unique_positions[last_index] - unique_positions[first_index] < min_seg_snps) or \
                    (min_seg_bp and snps_positions[last_index] - snps_positions[first_index] < min_seg_bp):
                score[i] -= inf_score
                check_optimal = False
        if check_optimal:
            for k in range(i):
                parameter_penalty = -1 / 2 * (best_boundaries_count[k] + 1) * (np.log(total_snps_count) + 1) * b_penalty

                z_penalty = 0
                if min_seg_snps or min_seg_bp:
                    last_index = candidate_numbers[i] + 1 - first_snp_number if i != candidates_count else -1
                    first_index = candidate_numbers[k] + 1 - first_snp_number
                    if (min_seg_snps and unique_positions[last_index] - unique_positions[first_index] < min_seg_snps) or \
                            (min_seg_bp and snps_positions[last_index] - snps_positions[first_index] < min_seg_bp):
                        z_penalty = -inf_score

                likelihood = score[k] + L[k + 1, i]
                candidate = likelihood + parameter_penalty + z_penalty
                if candidate > current_optimal:
                    current_optimal = candidate
                    score[i] = likelihood + z_penalty
                    kf = k
        if kf != -1:
            for j in range(0, kf):
                has_boundary_cache[i, j] = has_boundary_cache[kf, j]
            has_boundary_cache[i, kf] = 1
            for j in range(kf + 1, i):
                has_boundary_cache[i, j] = 0

        best_boundaries_count[i] = has_boundary_cache[i, :].sum()

    boundaries_indexes = [candidate_numbers[j] for j in range(candidates_count) if
                          has_boundary_cache[-1, j]]
    return boundaries_indexes


def segmentation_start():
    args = docopt(__doc__)
    if args['--test']:
        args['<file>'] = os.path.join(os.path.dirname(__file__), 'tests/test.tsv')

    schema = Schema({
        '<file>': And(
            Const(os.path.exists, error='Input file should exist'),
            Use(open),
            Use(lambda x: parse_input_file(x, int(args['--allele-reads-tr']), args["--force-sort"]))
        ),
        '--boundary-penalty': And(
            Use(float),
            Const(lambda x: x >= 0), error='Boundary penalty coefficient should be non negative number'
        ),
        '--min-seg-snps': And(
            Use(int),
            Const(lambda x: x >= 0), error='Min SNPs in segment parameter coefficient should be non negative integer'
        ),
        '--min-seg-bp': And(
            Use(int),
            Const(lambda x: x >= 0), error='Min segment length coefficient should be non negative integer'
        ),
        '--badmap': Or(
            Const(lambda x: x is None and not args['visualize']),
            And(
                Const(os.path.exists, error='Badmap file should exist'),
                Const(lambda x: os.access(x, os.R_OK), error='No read permission for badmap file')
            )),
        '--output': Or(
            And(
                Const(os.path.exists),
                Const(lambda x: os.access(x, os.W_OK), error='No write permissions')
            ),
            And(
                Const(lambda x: not os.path.exists(x)),
                Const(lambda x: os.access(os.path.dirname(x) if os.path.dirname(x) != '' else '.', os.W_OK),
                      error='No write permissions')
            ),
        ),
        '--states': Use(
            check_states, error='''Incorrect value for --states.
            Must be "," separated list of numbers or fractions in the form "x/y", each >= 1'''
        ),
        '--allele-reads-tr': And(
            Use(int),
            Const(lambda x: x >= 0), error='Allelic reads threshold must be a non negative integer'
        ),
        '--post-segment-filter': And(
            Use(int),
            Const(lambda x: x >= 0), error='Segments length filter (in SNPs) must be a non negative integer'
        ),
        '--subchr-filter': And(
            Use(int),
            Const(lambda x: x >= 0), error='Subchromosome post filter (in SNPs) must be a non negative integer'
        ),
        '--chr-min-snps': And(
            Use(int),
            Const(lambda x: x >= 0), error='Subchromosome pre filter (in SNPs) must be a non negative integer'
        ),
        '--atomic-region-size': And(
            Use(int),
            Const(lambda x: x >= 200), error='Atomic region size (in SNPs) must be a non negative integer'
        ),
        str: bool
    })
    try:
        args = schema.validate(args)
    except SchemaError as e:
        print(__doc__)
        exit('Error: {}'.format(e))

    snps_collection, chromosomes_order, full_name = args['<file>']
    file_name = os.path.splitext(os.path.basename(full_name))[0]
    if not args['visualize']:
        badmap_file_path = args['--output']
        if os.path.isdir(badmap_file_path):
            badmap_file_path += file_name + '.bed'

        verbose = not args['--quiet']
        mode = 'corrected'
        t = time.perf_counter()
        GS = GenomeSegmentator(snps_collection=snps_collection,
                               chromosomes_order=chromosomes_order,
                               out=badmap_file_path,
                               segmentation_mode=mode,
                               states=args['--states'],
                               b_penalty=args['--boundary-penalty'],
                               verbose=verbose,
                               allele_reads_tr=args['--allele-reads-tr'],
                               min_seg_snps=args['--min-seg-snps'],
                               min_seg_bp=args['--min-seg-bp'],
                               post_seg_filter=args['--post-segment-filter'],
                               atomic_region_size=args['--atomic-region-size'],
                               chr_filter=args['--chr-min-snps'],
                               subchr_filter=args['--subchr-filter']
                               )
        try:
            GS.estimate_BAD()
        except Exception as e:
            raise e
        print('Total time: {} s'.format(time.perf_counter() - t))
    else:
        badmap_file_path = args['--badmap']
    if args['--visualize'] or args['visualize']:
        init_from_snps_collection(snps_collection=snps_collection, BAD_file=badmap_file_path)
