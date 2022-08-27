__all__ = ['ChromosomesWrapper', 'GenomeSNPsHandler', 'pack', 'nucleotides', 'init_wrapper']

import pandas as pd
import numpy as np
nucleotides = ['A', 'T', 'G', 'C']

chr_l = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
         145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
         101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
         156040895, 57227415]


class ChromosomesWrapper:
    def __init__(self, chromosomes_df=None):
        if chromosomes_df is None:
            self.sorted_chromosomes = [f'chr{i}' for i in range(1, 23)] + ['chrX', 'chrY']
            self.chromosomes = dict(zip(self.sorted_chromosomes, chr_l))
        else:
            self.chromosomes = pd.Series(chromosomes_df.length.values,
                                         index=chromosomes_df.chromosome).to_dict()
            self.sorted_chromosomes = chromosomes_df.chromosome.tolist()


def init_wrapper(wrapper):
    if wrapper is None:
        return ChromosomesWrapper()
    else:
        return wrapper


def pack(values):
    return '\t'.join(map(str, values)) + '\n'


class GenomeSNPsHandler:
    def __init__(self, data: pd.DataFrame, chrom_wrapper: ChromosomesWrapper):
        self.data = {}
        gb = data.groupby(['chr'])
        self.chromosomes_order = data['chr'].unique()
        for chromosome, group_df in [(group, gb.get_group(group)) for group in gb.groups]:
            if chromosome in chrom_wrapper.chromosomes:
                self.data[chromosome] = ChromosomeSNPsHandler.from_df(chromosome, group_df)


class ChromosomeSNPsHandler:
    def __init__(self, chromosome, data: np.ndarray):
        if not type(data) is np.ndarray:
            raise ValueError('Not a numpy array provided.')
        if data.shape[0] != 3:
            raise ValueError(f'Wrong data shape {data.shape}')
        self.chromosome = chromosome
        self.data = data

    @classmethod
    def from_df(cls, chromosome, counts_df):
        numpy_df = counts_df[['start', 'ref_counts', 'alt_counts']].transpose().to_numpy()
        return cls(chromosome, numpy_df)

