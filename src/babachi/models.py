import math
import pandas as pd
import numpy as np
from typing import List

from babachi.helpers import pack
from babachi.chrom_wrapper import ChromosomesWrapper, init_wrapper


class BADSegmentsContainer:
    allowed_fields = ['boundaries_positions', 'BAD_estimations', 'likelihoods',
        'snps_counts', 'covers', 'snp_id_counts']

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
                        current_position = 0
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
    def __init__(self, chrom, start, end, BAD, likelihoods, snps_count, total_cover, snp_id_count):
        self.chr = chrom
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


def filter_segments(segments: List[BADSegment], post_seg_filter: int = None):
    for segment in segments:
        if post_seg_filter is None or segment.snp_id_count >= post_seg_filter:
            if segment.BAD != 0:
                yield segment


## SNPs handlers

class ChromosomeSNPsHandler:
    def __init__(self, chromosome, positions: np.ndarray, read_counts: np.ndarray):
        assert positions.shape[0] == read_counts.shape[0]
        if not type(read_counts) is np.ndarray:
            raise ValueError('Not a numpy array provided.')
        if read_counts.shape[1] != 2:
            raise ValueError(f'Wrong data shape {read_counts.shape}')
        self.chromosome = chromosome
        self.read_counts = read_counts
        self.positions = positions


class GenomeSNPsHandler:
    def __init__(self, *chrom_handlers: ChromosomeSNPsHandler):
        self.data = {}
        self.chromosomes_order = []
        for chrom_handler in chrom_handlers:
            self.data[chrom_handler.chromosome] = chrom_handler
            self.chromosomes_order.append(chrom_handler.chromosome)
    
    def __getitem__(self, item) -> ChromosomeSNPsHandler:
        return self.data[item]

    @classmethod
    def from_df(cls, data: pd.DataFrame, chrom_wrapper: ChromosomesWrapper=None):
        chrom_wrapper = init_wrapper(chrom_wrapper)
        gb = data.groupby(['chr'])
        snp_handlers = []
        for chromosome, group_df in [(group, gb.get_group(group)) for group in gb.groups]:
            if chromosome in chrom_wrapper.chromosomes:
                positions = group_df['start'].values.astype(np.uint32)
                read_counts = group_df[['ref_counts', 'alt_counts']].to_numpy().astype(np.float32)
                snp_handlers.append(ChromosomeSNPsHandler(chromosome, positions, read_counts))
        return cls(*snp_handlers)
