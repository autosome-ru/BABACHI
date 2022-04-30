__all__ = ['ChromosomesWrapper', 'pack', 'nucleotides', 'init_wrapper']

import pandas as pd

nucleotides = ['A', 'T', 'G', 'C']

chr_l = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
         145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
         101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
         156040895, 57227415]


class ChromosomesWrapper:
    def __init__(self, chromosomes_df=None):
        if chromosomes_df is None:
            self.sorted_chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
            self.chromosomes = dict(zip(self.sorted_chromosomes, chr_l))
        else:
            self.chromosomes = pd.Series(chromosomes_df.length.values,
                                         index=chromosomes_df.chromosomes)
            self.sorted_chromosomes = chromosomes_df.chromosomes.tolist()


def init_wrapper(wrapper):
    if wrapper is None:
        return ChromosomesWrapper()
    else:
        return wrapper


def pack(values):
    return '\t'.join(map(str, values)) + '\n'


