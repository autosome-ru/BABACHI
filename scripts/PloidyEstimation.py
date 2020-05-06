"""

Usage:
    segmentation <file> [options]

Arguments:
    <file>     Path to input file in vcf format (chr pos ref_base alt_base).

Options:
    -h, --help                  Show help.
    -V, --version               Show version.
    -v, --verbose               Print additional messages during work time.
    -a, --add                   Create additional file with SNPs intersection.
    --log <path>                Path to a verbose appending log.
    -O <path>, --output <path>  Output directory or file path. [default: ./]
"""

import math
import numpy as np
import os.path
from schema import Schema, And, Use, SchemaError, Const, Optional, Or
import time
from scripts.helpers import ChromPos, pack
from abc import ABC, abstractmethod
from docopt import docopt


class BADSegmentsContainer:
    allowed_fields = ['boundaries_positions', 'BAD_estimations', 'likelihoods', 'snps_counts', 'total_cover']

    def __init__(self, **kwargs):
        self.BAD_estimations = []  # estimated BADs for split segments
        self.likelihoods = []  # likelihoods of split segments for each BAD
        self.snps_counts = []  # number of snps in segments
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


class Segmentation(ABC):
    def __init__(self):
        self.dtype = np.float64
        self.segments_container = BADSegmentsContainer()
        self.S = None  # snp-wise log-likelihoods for each BAD
        self.C = None  # snp-wise marginal log-likelihoods
        self.P = None  # segment-wise log-likelihoods for each BAD
        self.L = None  # segment-wise marginal log-likelihoods

        self.boundaries_indexes = None
        self.boundary_numbers = None
        self.snps_positions = []
        self.total_cover = None

        self.total_snps = None
        self.last_snp_number = None
        self.candidate_numbers = None
        self.candidates_count = None
        self.sub_chromosome = None
        self.score = None  # score[i] = best log-likelihood among all segmentations of snps[0,i]
        self.has_boundary = None  # bool boundaries, len=total_snps.
        self.best_boundaries_count = None  # best_boundaries_count[i] = number of boundaries
        # before ith snp in best segmentation

    @staticmethod
    def get_norm(p, N, trim_cover):
        result = 0
        current_multiplier = 1
        denominator_multiplier = 1
        for k in range(trim_cover):
            result += current_multiplier * np.power(p, N - k) * np.power(1 - p, k) / denominator_multiplier
            current_multiplier *= (N - k)
            denominator_multiplier *= k + 1

        return -result

    def log_likelihood(self, N, X, i):
        """
        allele_reads_tr <= X <= N/2
        """
        p = 1.0 / (1.0 + i)
        log_norm = np.log1p(self.get_norm(p, N, self.sub_chromosome.gs.allele_reads_tr) +
                            self.get_norm(1 - p, N, self.sub_chromosome.gs.allele_reads_tr))
        if (
                self.sub_chromosome.gs.mode == 'corrected' and N == 2 * X) or self.sub_chromosome.gs.mode == 'binomial':
            return X * np.log(p) + (N - X) * np.log(1 - p) + np.log(
                self.sub_chromosome.gs.prior[i]) - log_norm
        elif self.sub_chromosome.gs.mode == 'corrected':
            return X * np.log(p) + (N - X) * np.log(1 - p) + np.log(self.sub_chromosome.gs.prior[i]) \
                   + np.log1p(i ** (2 * X - N)) - log_norm

    def get_P(self, first, last):
        if last - first == 1:
            return self.sub_chromosome.P_initial[:, last]
        else:
            return np.sum(self.sub_chromosome.P_initial[:, first + 1:last + 1], axis=1)

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
                first = -1
            else:
                first = self.candidate_numbers[j - 1]
            if j == self.candidates_count:
                last = self.last_snp_number
            else:
                last = self.candidate_numbers[j]

            S[:, j] = self.get_P(first, last)
        self.S = S
        self.P = P

    def modify_P(self):
        # print('Constructing p-matrix')
        self.C = np.cumsum(self.S, axis=1)
        for j in range(self.candidates_count + 1):
            if j == 0:
                subtract = 0
            else:
                subtract = self.C[:, j - 1]
            for k in range(j, self.candidates_count + 1):
                self.P[:, j, k] = self.C[:, k] - subtract

    def modify_L(self):
        # print('Constructing L')
        Q = np.sort(self.P, axis=0)
        self.L[:, :] = Q[-1, :, :] + np.log1p(np.sum(np.exp(Q[:-2, :, :] - Q[-1, :, :]), axis=0))

    def get_parameter_penalty(self, boundaries, alphabet):
        k = boundaries * alphabet
        if isinstance(self, AtomicRegionSegmentation):
            N = self.total_snps
        else:
            N = self.sub_chromosome.gs.unique_snp_positions

        if self.sub_chromosome.gs.b_penalty == 'CAIC':
            return -1 / 2 * k * (np.log(N) + 1)
        elif self.sub_chromosome.gs.b_penalty == 'SQRT':
            return -1 / 2 * k * (np.sqrt(N) + 1)
        elif self.sub_chromosome.gs.b_penalty == 'CBRT':
            return -1 / 2 * k * (N ** (1 / 3) + 1)
        else:
            raise ValueError(self.sub_chromosome.b_penalty)

    def initialize_boundaries_arrays(self):
        self.score = [0] * (self.candidates_count + 1)
        self.has_boundary = [False] * self.candidates_count
        self.best_boundaries_count = [0] * (self.candidates_count + 1)

    @abstractmethod
    def find_optimal_boundaries(self):
        pass

    def estimate(self):
        self.construct_likelihood_matrices()
        self.modify_P()
        self.modify_L()
        self.find_optimal_boundaries()


class AtomicRegionSegmentation(Segmentation):
    def __init__(self, sub_chromosome, start, end):
        super().__init__()
        self.sub_chromosome = sub_chromosome
        self.start_snp_index = start
        self.end_snp_index = end
        self.total_snps = end - start + 1
        self.positions = sub_chromosome.snps_positions[
                         sub_chromosome.candidate_numbers[start]:sub_chromosome.candidate_numbers[end - 1] + 2]
        self.total_cover = sum(x[1] + x[2] for x in sub_chromosome.snps_array[
                                                    sub_chromosome.candidate_numbers[start]:
                                                    sub_chromosome.candidate_numbers[end - 1] + 2])
        self.candidate_numbers = sub_chromosome.candidate_numbers[start:end + 1]
        self.candidates_count = end - start
        if self.end_snp_index == sub_chromosome.candidates_count:
            self.last_snp_number = sub_chromosome.total_snps_count - 1
        else:
            self.last_snp_number = sub_chromosome.candidate_numbers[end + 1] - 1

        self.initialize_boundaries_arrays()

    # print(self.start, self.end, self.candidates_count, len(self.positions))

    def find_optimal_boundaries(self):
        # print('Constructing boundaries')
        for i in range(self.candidates_count + 1):
            self.score[i] = self.L[0, i]

            kf = -1
            current_optimal = self.score[i]

            for k in range(i):
                parameter_penalty = self.get_parameter_penalty(self.best_boundaries_count[k] + 1, len(
                    self.sub_chromosome.chromosome_segmentation.BAD_list))

                likelihood = self.score[k] + self.L[k + 1, i]
                candidate = likelihood + parameter_penalty
                if candidate > current_optimal:
                    current_optimal = candidate
                    self.score[i] = likelihood
                    kf = k
            if kf != -1:
                self.has_boundary[kf] = True
            for j in range(kf + 1, i):
                self.has_boundary[j] = False

            self.best_boundaries_count[i] = self.has_boundary.count(True)

        self.boundaries_indexes = [self.candidate_numbers[i] for i in range(self.candidates_count) if
                                   self.has_boundary[i]]


class SubChromosomeSegmentation(Segmentation):  # sub_chromosome
    def __init__(self, genome_segmentator, chromosome_segmentation, snps_array, part):
        super().__init__()

        self.gs = genome_segmentator
        self.chromosome_segmentation = chromosome_segmentation
        self.sub_chromosome = self
        self.index_in_chromosome = part
        self.snps_array = snps_array

        self.total_snps = len(self.snps_array)
        self.total_cover = sum(x[1] + x[2] for x in self.snps_array)
        self.start_snp_index = 0
        self.end_snp_index = (self.total_snps - 1) - 1  # index from 0, and #boundaries = #snps - 1
        self.candidate_numbers = [i for i in range(self.total_snps - 1) if
                                  self.snps_array[i][0] != self.snps_array[i + 1][0]]
        self.candidates_count = len(self.candidate_numbers)
        self.last_snp_number = self.total_snps - 1

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
        current_snip = -1

        S = np.zeros((len(self.chromosome_segmentation.BAD_list), self.total_snps), dtype=self.dtype)
        for j in range(0, self.total_snps):

            pos, ref_c, alt_c = self.snps_array[self.start_snp_index + j]
            current_snip += 1
            N = ref_c + alt_c
            X = min(ref_c, alt_c)

            self.snps_positions.append(pos)

            for i in range(len(self.chromosome_segmentation.BAD_list)):
                assert (self.chromosome_segmentation.BAD_list[i] > 0)
                S[i, current_snip] = self.log_likelihood(N, X, self.chromosome_segmentation.BAD_list[i])
        self.P_initial = S

    def find_optimal_boundaries(self):
        for i in range(self.candidates_count + 1):
            self.score[i] = self.L[0, i]

            last_boundary_index = -1
            current_optimal = self.score[i]

            for k in range(i):

                parameter_penalty = self.get_parameter_penalty(self.best_boundaries_count[k] + 1,
                                                               len(self.chromosome_segmentation.BAD_list))

                likelihood = self.score[k] + self.L[k + 1, i]
                candidate = likelihood + parameter_penalty
                if candidate > current_optimal:
                    current_optimal = candidate
                    self.score[i] = likelihood
                    last_boundary_index = k
            if last_boundary_index != -1:
                self.has_boundary[last_boundary_index] = True
            for j in range(last_boundary_index + 1, i):
                self.has_boundary[j] = False

            self.best_boundaries_count[i] = [int(x) for x in self.has_boundary].count(1)
            assert ([int(x) for x in self.has_boundary].count(1) == self.best_boundaries_count[i])

        self.boundaries_indexes = [self.candidate_numbers[i] for i in range(self.candidates_count) if
                                   self.has_boundary[i]]
        self.boundary_numbers = [-1] + [i for i in range(self.candidates_count) if self.has_boundary[i]] + [
            self.candidates_count]
        cumulative_counts = [0] + [x + 1 for x in self.boundaries_indexes] + [self.last_snp_number + 1]
        self.segments_container.snps_counts = [cumulative_counts[i + 1] - cumulative_counts[i] for i in
                                               range(len(cumulative_counts) - 1)]
        self.segments_container.covers = [
            sum(x[1] + x[2] for x in self.snps_array[cumulative_counts[i]:cumulative_counts[i + 1]]) for i in
            range(len(cumulative_counts) - 1)]

        for i in range(len(self.has_boundary)):
            if self.has_boundary[i]:
                self.segments_container.boundaries_positions.append(
                    (self.snps_positions[self.candidate_numbers[i]] + self.snps_positions[
                        self.candidate_numbers[i] + 1]) / 2)

    def estimate_BAD(self):
        for n in range(len(self.boundary_numbers) - 1):
            first = self.boundary_numbers[n] + 1
            last = self.boundary_numbers[n + 1]
            likelihoods = self.P[:, first, last] - self.L[first, last]
            estimated_BAD_index = np.argmax(likelihoods)
            self.segments_container.BAD_estimations.append(self.chromosome_segmentation.BAD_list[estimated_BAD_index])
            self.segments_container.likelihoods.append(list(likelihoods))

    def estimate_sub_chr(self):
        if self.total_snps == 0:
            return

        self.construct_initial_likelihood_matrices()
        tuples = self.split_into_overlapping_regions(self.candidates_count + 1,
                                                     self.chromosome_segmentation.atomic_region_length,
                                                     self.chromosome_segmentation.overlap)
        boundary_set = set()
        # print("{} segments: {}".format(len(tuples), tuples))
        counter = 0
        for first, last in tuples:
            counter += 1
            if self.gs.verbose:
                print(
                    'Making {} out of {} segments from {} to {} for {} (part {} of {}).'.format(
                        counter, len(tuples), first,
                        last, self.chromosome_segmentation.chromosome,
                        self.index_in_chromosome,
                        len(self.chromosome_segmentation.get_sub_chromosomes_slices())))
            atomic_region_segmentation = AtomicRegionSegmentation(self, first, last)
            atomic_region_segmentation.estimate()
            boundary_set |= set(atomic_region_segmentation.boundaries_indexes)
        self.candidate_numbers = sorted(list(boundary_set))
        self.candidates_count = len(self.candidate_numbers)
        if self.gs.verbose:
            print('SNPs in part: {}'.format(len(self.snps_positions)))
        # print('{} candidates'.format(self.candidates_count))

        self.initialize_boundaries_arrays()

        self.estimate()
        self.estimate_BAD()
        if self.gs.verbose:
            print(
                '\n'.join(map(str, zip(self.segments_container.BAD_estimations, self.segments_container.snps_counts))))


class ChromosomeSegmentation:  # chromosome
    def __init__(self, genome_segmentator, chromosome, length=0):
        self.gs = genome_segmentator

        self.chromosome = chromosome  # name
        self.length = length  # length, bp

        (self.snps_array,  #
         self.total_snps_count,
         self.snps_positions) = self.read_file()  # unpack
        if self.total_snps_count == 0:
            return
        self.unique_snp_positions = len(set(snp[0] for snp in self.snps_array))
        self.total_read_coverage = sum(x[1] + x[2] for x in self.snps_array)
        self.critical_gap_factor = 1 - 10 ** (- 1 / np.sqrt(self.total_snps_count))
        self.effective_length = self.snps_positions[-1] - self.snps_positions[0]

        self.segments_container = BADSegmentsContainer()

    def read_file(self):
        count = 0
        snps = []
        positions = []
        for line in self.gs.file:
            line_tuple = self.unpack_line_or_false(line)
            if not line_tuple:
                continue
            count += 1
            pos, _, _ = line_tuple
            snps.append(line_tuple)
            positions.append(pos)
        self.gs.file.seek(0)
        return snps, count, positions

    def unpack_line_or_false(self, line):
        try:
            # FIXME int pos ref_c alt_c
            chr, pos, ref_c, alt_c = line.strip('\n').split('\t')[:4]
        except ValueError:
            return False
        if chr != self.chromosome:
            return False
        return pos, ref_c, alt_c

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

    #  TODO: (N+1)/N paradox
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
        if not self.total_snps_count or self.total_snps_count < self.gs.snp_per_chr_tr:
            return

        start_t = time.clock()
        self.adjust_critical_gap()

        #  boundary for first snp
        if self.snps_positions[0] <= self.critical_gap_factor * self.effective_length:
            self.segments_container.boundaries_positions.append(self.snps_positions[0])
        else:
            self.segments_container.boundaries_positions.append((1, self.snps_positions[0]))
        if self.gs.verbose:
            print('Distance splits {}'.format(self.get_sub_chromosomes_slices()))

        for part, (st, ed) in enumerate(self.get_sub_chromosomes_slices(), 1):
            # check
            if len(set(self.snps_positions[st: ed])) <= self.gs.min_segment_length:
                self.segments_container += BADSegmentsContainer(
                    boundaries_positions=[],
                    BAD_estimation=[0],
                    likelihoods=[[0] * len(self.gs.BAD_list)],
                    snps_counts=[ed - st],
                    covers=[0],
                )
            else:
                sub_chromosome = SubChromosomeSegmentation(self.gs, self, self.snps_array[st: ed], part)
                sub_chromosome.estimate_sub_chr()

                self.segments_container += sub_chromosome.segments_container
            if ed != self.total_snps_count:
                self.segments_container.boundaries_positions += [(self.snps_positions[ed - 1], self.snps_positions[ed])]

        #  boundary for last snp
        if self.length - self.snps_positions[-1] <= self.critical_gap_factor * self.effective_length:
            self.segments_container.boundaries_positions.append(self.length)
        else:
            self.segments_container.boundaries_positions.append((self.snps_positions[-1] + 1, self.length))
        if self.gs.verbose:
            print('\nTotal SNPs: {},'
                  '\nEstimated BADs: {},'
                  '\nSNP counts {}'
                  '\nCritical gap {:.0f}'
                  '\nBoundaries distances: {}'
                  .format(len(self.snps_positions), self.segments_container.BAD_estimations,
                          self.segments_container.snps_counts,
                          self.critical_gap_factor * self.effective_length,
                          list(map(lambda x: (x, 1) if isinstance(x, (int, float)) else (x[0], x[1] - x[0]),
                                   self.segments_container.boundaries_positions))))
            print('{} time: {} s\n'.format(self.chromosome, time.clock() - start_t))


class GenomeSegmentator:  # gs
    def __init__(self, file, out, segmentation_mode='corrected', extra_states=None, b_penalty='CAIC',
                 prior=None, verbose=False, additional_file_path=False, log_file_path=None):

        self.additional_file = additional_file_path  # path to file with SNPs intersection
        self.log_file_path = log_file_path
        self.verbose = verbose
        self.mode = segmentation_mode  # 'corrected' or 'binomial'
        self.b_penalty = b_penalty  # boundary penalty mode ('CAIC', ')
        if extra_states is None:
            self.BAD_list = [1, 2, 3, 4, 5]
        else:
            self.BAD_list = sorted([1, 2, 3, 4, 5] + extra_states)
        if prior is None:
            self.prior = dict(zip(self.BAD_list, [1] * len(self.BAD_list)))
        else:
            self.prior = prior

        self.file = file  # input file in .vcf format
        self.out = open(out, 'w')  # output file in .bed format

        self.snp_per_chr_tr = 100  # minimal number of snps in chromosome to start segmentation
        self.allele_reads_tr = 5  # "minimal read count on each allele" snp filter
        self.atomic_region_length = 600  # length of an atomic region in snps
        self.overlap = 300  # length of regions overlap in snps
        self.min_segment_length = 2  # minimal segment length in snps
        self.chromosomes = sorted(list(ChromPos.chrs.keys()))  # {'chr1': length_in_bp, ...}

        self.chr_segmentations = []  # list of ChromosomeSegmentation instances

        for chromosome in self.chromosomes:
            chr_segmentation = ChromosomeSegmentation(self, chromosome, ChromPos.chrs[chromosome])
            if self.verbose:
                print('{} total SNP count: {}'.format(chromosome, chr_segmentation.total_snps_count))
            self.chr_segmentations.append(chr_segmentation)

    @staticmethod
    def append_BAD_segments(chr, total_snps_count_tr):
        segments_to_write = []
        cur = None
        counter = 0
        sc = chr.segments_container
        if chr.total_snps_count >= total_snps_count_tr:
            for boundary in chr.boundaries_positions:
                if cur is None:
                    if isinstance(boundary, tuple):
                        cur = boundary[1]
                    else:
                        cur = 1
                elif isinstance(boundary, tuple):
                    segments_to_write.append(
                        [chr.chromosome, cur, boundary[0] + 1, sc.BAD_estimations[counter]] + sc.likelihoods[
                            counter] +
                        [sc.snps_counts[counter], sc.covers[counter]])
                    cur = boundary[0] + 1
                    segments_to_write.append([chr.chromosome, cur, boundary[1], 0] + [0] * len(chr.BAD_list) + [0, 0])
                    cur = boundary[1]
                    counter += 1
                else:
                    segments_to_write.append(
                        [chr.chromosome, cur, math.floor(boundary) + 1, sc.BAD_estimations[counter]] +
                        sc.likelihoods[counter] +
                        [sc.snps_counts[counter], sc.covers[counter]])
                    cur = math.floor(boundary) + 1
                    counter += 1

        return segments_to_write

    def write_BAD_to_file(self, chr):
        segments = self.append_BAD_segments(chr, self.snp_per_chr_tr)

        filtered_segments = self.filter_segments(segments, self.min_segment_length)
        for segment in filtered_segments:
            if segment[3] == 0:  # BAD == 0
                continue
            self.out.write(pack(segment))

    # noinspection PyTypeChecker
    def estimate_BAD(self):
        self.out.write(
            pack(['#chr', 'start', 'end', 'BAD'] + ['Q{:.2f}'.format(BAD) for BAD in self.BAD_list] + ['SNP_count',
                                                                                                       'sum_cover']))
        for j in range(len(self.chr_segmentations)):
            chromosome = self.chr_segmentations[j]
            chromosome.estimate_chr()
            self.write_BAD_to_file(chromosome)
            self.chr_segmentations[j] = None

    def filter_segments(self, segments, snp_number_tr=2):
        is_bad_left = False
        is_bad_segment = False
        for k in range(len(segments)):
            if segments[k][4 + len(self.BAD_list)] <= snp_number_tr and segments[k][3] != 0:  # если k сегмент "плохой"
                if is_bad_segment:  # если k-1 тоже "плохой"
                    is_bad_left = True
                    for j in range(3, 4 + len(self.BAD_list)):
                        segments[k - 1][j] = 0
                else:
                    is_bad_left = False
                is_bad_segment = True  # текущий сегмент плохой, следующий шаг цикла
            else:  # k сегмент хороший
                if is_bad_segment and not is_bad_left and k > 1:  # а k-1 плохой и k-2 хороший
                    # if segments[k][3] < segments[k - 1][3] and segments[k - 2][3] < segments[k - 1][3]:
                    #     # если BAD k-1 сегмента больше BAD k-2 и k сегментов
                    #     if segments[k][3] > segments[k - 2][3]:  # если BAD k сегмента больше BAD k-2
                    #         segments[k - 1][3] = segments[k][3]  # присвоить BAD k сегмента
                    #     else:  # если BAD k-2 сегмента больше BAD k
                    #         segments[k - 1][3] = segments[k - 2][3]  # присвоить BAD k-2 сегмента
                    #
                    #     for j in range(4, 4 + len(self.i_list)):
                    #         segments[k - 1][j] = 0
                    is_bad_left = True
                if is_bad_left and is_bad_segment:
                    for j in range(3, 4 + len(self.BAD_list)):
                        segments[k - 1][j] = 0
                    is_bad_left = True

                is_bad_segment = False  # текущий сегмент хороший, следующий шаг цикла

        return segments


# FIXME
def check_input_file_format():
    # for lines in opened_file:
    pass


def segmentation_start():
    # TODO: Global version (here and in setup.py)
    args = docopt(__doc__, version='BAD segmentation v0.1')
    schema = Schema({
        '<file>': And(
            Const(os.path.exists, error='Input file should exist'),
            Use(open, error='Input file should be readable'),
            Const(check_input_file_format, error='Wrong input file format')
        ),
        '--output': And(
            Const(os.path.exists, error='Output path should exist'),
            Const(lambda x: os.access(x, os.W_OK), error='No write permissions')
        ),
        Optional('--log'): Or(
            None, And(
                Const(os.path.exists, error='Log path should exist'),
                Const(lambda x: os.access(x, os.W_OK), error='No write permissions for log file')
            )
        ),
        str: bool
    })
    try:
        args = schema.validate(args)
    except SchemaError as e:
        print(__doc__)
        exit('Error: {}'.format(e))

    input_file = args['<file>']
    file_name = os.path.splitext(os.path.basename(input_file.index_in_chromosome))[0]

    log_file_path = args['--log']
    if log_file_path and os.path.isdir(log_file_path):
        log_file_path += file_name + '.log'

    output_file_path = args['--output']
    if os.path.isdir(output_file_path):
        output_file_path += file_name + '.bed'

    verbose = args['--verbose']
    mode = 'corrected'
    b_penalty = 'CAIC'
    states = [4 / 3, 1.5, 2.5, 6]
    t = time.clock()
    GS = GenomeSegmentator(file=input_file,
                           out=output_file_path,
                           segmentation_mode=mode,
                           extra_states=states,
                           b_penalty=b_penalty,
                           verbose=verbose,
                           additional_file_path=args['--add'],
                           log_file_path=log_file_path
                           )
    try:
        GS.estimate_BAD()
    except Exception as e:
        raise e
    if verbose:
        print('Total time: {} s'.format(time.clock() - t))
