import pandas as pd
from babachi.logging import root_logger
from babachi.chrom_wrapper import init_wrapper
from urllib.request import Request, urlopen
import os
import re


nucleotides = ['A', 'T', 'G', 'C']
df_header = ['chr', 'start', 'end', 'ID', 'ref', 'alt', 'ref_counts', 'alt_counts', 'sample_id']


class InputParser:
    def __init__(self, allele_reads_tr=5, snp_strategy='SEP', force_sort=False, to_filter=True,
                 chromosomes_wrapper=None, filter_no_rs=False, logger=None):
        self.allele_reads_tr = allele_reads_tr
        self.to_filter = to_filter
        self.force_sort = force_sort
        if logger is None:
            logger = root_logger
        self.logger = logger
        self.snp_strategy = snp_strategy
        self.filter_no_rs = filter_no_rs
        self.chromosomes_wrapper = init_wrapper(chromosomes_wrapper)
        if force_sort:
            self.chromosomes_order = self.chromosomes_wrapper.sorted_chromosomes
        else:
            self.chromosomes_order = []

    def _filter_record(self, record, line_number, sample_id_list):
        if sample_id_list is not None:
            samples = [record.samples[sample_id] for sample_id in sample_id_list]
        else:
            samples = record.samples
        if record.chrom not in self.chromosomes_wrapper.chromosomes:
            self.logger.warning(f'Chromosome length for {record.chrom} in line #{line_number} not available')
            return
        if self.to_filter:
            if len(record.alts) != 1:
                return
            if record.ref not in nucleotides or record.alts[0] not in nucleotides:
                return
            if self.filter_no_rs and record.id == '.':
                return
        result = []
        ref_read_sum = 0
        alt_read_sum = 0
        filter_out = True
        names = []
        for sample_id, sample in samples.items():
            sample_ref_read_count, sample_alt_read_count = map(int, sample['AD'])
            if self.to_filter:
                if min(sample_ref_read_count, sample_alt_read_count) < self.allele_reads_tr:
                    continue
                if '/'.join(map(str, sample['GT'])) != '0/1':
                    continue
            filter_out = False
            if self.snp_strategy == 'ADD':
                ref_read_sum += sample_ref_read_count
                alt_read_sum += sample_alt_read_count
                names.append(sample.name)
            elif self.snp_strategy == 'SEP':
                result.append(
                    (sample_ref_read_count, sample_alt_read_count, sample.name)
                )
            else:
                raise ValueError
        if filter_out:
            return
        if self.snp_strategy == 'ADD':
            result.append(
                (ref_read_sum, alt_read_sum, ','.join(names))
            )
        if record.chrom not in self.chromosomes_order:
            self.chromosomes_order.append(record.chrom)
        return result

    @staticmethod
    def df_to_counts(df):
        return list(zip(*df[['start', 'ref_counts', 'alt_counts']].transpose().to_numpy()))

    def read_bed(self, file_path) -> pd.DataFrame:
        """
        :param file_path: input bed file path
        :return: pd.DataFrame
        """
        df = pd.read_table(file_path, header=None, comment='#')
        df = df[df.columns[:len(df_header)]]
        df.columns = df_header
        if self.to_filter:
            df = df[df['chr'].isin(self.chromosomes_wrapper.chromosomes)]
            df = df[df[['ref_counts', 'alt_counts']].min(axis=1) >= self.allele_reads_tr]
        return df

    @staticmethod
    def check_if_vcf(file_path):
        result = True
        with open(file_path) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                if re.match(r'^chr(\d+|X|Y)\t\d+\t\d+\t', line):
                    result = False
                f.seek(0)
                break
        return result

    def read_file(self, file_path, samples_list=None) -> pd.DataFrame:
        is_vcf = samples_list is not None or self.check_if_vcf(file_path)
        if is_vcf:
            self.logger.debug('Reading as VCF file')
            return self.read_vcf(file_path, sample_list=samples_list)
        else:
            self.logger.debug('Reading as BED file')
            return self.read_bed(file_path)

    @staticmethod
    def check_record(record, previous_line):
        if previous_line is not None:
            if record.chrom == previous_line.chrom:
                if record.start < previous_line.start:
                    raise ValueError(f'VCF file is not sorted. Please sort input file.')
        return record

    def read_vcf(self, file_path, sample_list=None) -> pd.DataFrame:
        """
        :param sample_list: optional, list of sample names or sample IDs to work with
        :param file_path: input VCF file
        :return: None pd.DataFrame
        """
        self.logger.info('Reading input file...')
        try:
            from pysam import VariantFile
        except ImportError as e:
            print(f'Please install pysam package (https://pysam.readthedocs.io/en/latest/installation.html)')
            raise e
        vcfReader = VariantFile(file_path, 'r')
        if sample_list is None:
            sample_indices = None
        elif all(isinstance(sample, int) for sample in sample_list):
            sample_indices = sample_list
        else:
            sample_indices = []
            for sample in sample_list:
                sample_index = vcfReader.header.samples.index(sample)
                if sample_index == -1:
                    raise ValueError('Error: Sample {} was not found in header'.format(sample))
                sample_indices.append(vcfReader.header.samples[sample_index])
        previous_line = None
        result = []
        df_columns = df_header
        for line_number, record in enumerate(vcfReader.fetch(), 1):
            previous_line = self.check_record(record, previous_line)
            filter_result = self._filter_record(record, line_number, sample_indices)
            if filter_result:
                for counts in filter_result:
                    result.append([record.chrom, record.start,
                                   record.stop, record.id, record.ref, record.alts[0], *counts])
        return pd.DataFrame.from_records(result, columns=df_columns)

def pack(values):
    return '\t'.join(map(str, values)) + '\n'


def craft_prior(states, string, p):
    if string == 'uniform':
        return None
    minimum_ploidy = {
        1: 2,
        4 / 3: 7,
        3 / 2: 5,
        2: 3,
        5 / 2: 7,
        3: 4,
        4: 5,
        5: 6,
        6: 7,
    }
    return {
        state:
            p ** (minimum_ploidy[state] - 1)
        for state in states
    }


# IO functions
def read_url_file(url):
    url_request = Request(url)
    return urlopen(url_request)

def make_file_path_from_dir(out_path, file_name, ext='badmap.bed'):
    if os.path.isdir(out_path):
        return os.path.join(out_path, f'{file_name}.{ext}')
    else:
        return out_path


def read_snps_file(file_path, chrom_sizes=None, snp_strategy='SEP', samples_list=None,
                   allele_reads_tr=5, force_sort=False, to_filter=False, filter_no_rs=False):
    if chrom_sizes is not None:
        chrom_sizes = pd.read_table(
            chrom_sizes, 
            header=None, 
            names=['chromosome', 'length']
        )
    chromosomes_wrapper = init_wrapper(chrom_sizes)
    input_parser = InputParser(
        allele_reads_tr=allele_reads_tr,
        force_sort=force_sort,
        to_filter=to_filter,
        filter_no_rs=filter_no_rs,
        snp_strategy=snp_strategy,
        chromosomes_wrapper=chromosomes_wrapper,
    )
    file = input_parser.read_file(file_path=file_path, samples_list=samples_list)
    return file, chromosomes_wrapper
