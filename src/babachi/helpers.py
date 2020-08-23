__all__ = ['ChromosomePosition', 'pack', 'parse_input_file']

chr_l = [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973,
         145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718,
         101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468,
         156040895, 57227415]


class ChromosomePosition:
    sorted_chromosomes = ['chr' + str(i) for i in range(1, 23)] + ['chrX', 'chrY']
    chromosomes = dict(zip(sorted_chromosomes, chr_l))
    genome_length = sum(chr_l)

    def __init__(self, chr, pos):
        if chr not in self.chromosomes:
            raise ValueError("Not in valid chromosomes {}".format(chr))
        self.chr = chr
        self.pos = int(pos)

    def __lt__(self, other):
        if self.chr == other.chr:
            return self.pos < other.pos
        else:
            return self.chr < other.chr

    def __gt__(self, other):
        if self.chr == other.chr:
            return self.pos > other.pos
        else:
            return self.chr > other.chr

    def __le__(self, other):
        if self.chr == other.chr:
            return self.pos <= other.pos
        else:
            return self.chr <= other.chr

    def __ge__(self, other):
        if self.chr == other.chr:
            return self.pos >= other.pos
        else:
            return self.chr >= other.chr

    def __eq__(self, other):
        return (self.chr, self.pos) == (other.chr, other.pos)

    def __ne__(self, other):
        return (self.chr, self.pos) != (other.chr, other.pos)

    def distance(self, other):
        if self.chr != other.chr:
            return float('inf')
        return abs(self.pos - other.pos)


def pack(values):
    return '\t'.join(map(str, values)) + '\n'


def parse_input_file(opened_file, allele_reads_tr=0, force_sort=False):
    snps_collection = {chromosome: [] for chromosome in ChromosomePosition.chromosomes}
    chromosomes_order = []
    for line_number, line in enumerate(opened_file, 1):
        if line[0] == '#':
            continue
        else:
            parsed_line = line.strip().split('\t')
            if parsed_line[0] not in ChromosomePosition.chromosomes:
                print('Invalid chromosome name: {} in line #{}'.format(parsed_line[0], line_number))
                return False
            try:
                if int(parsed_line[1]) <= 0 or int(parsed_line[5]) < 0 or int(parsed_line[6]) < 0:
                    raise ValueError
            except ValueError:
                print('Position, Reference allele read counts, Alternative allele read counts must be'
                      ' a non-negative integer in line #{}'.format(line_number))
        ref_read_count = int(parsed_line[5])
        alt_read_count = int(parsed_line[6])
        if parsed_line[0] not in chromosomes_order:
            chromosomes_order.append(parsed_line[0])
        if ref_read_count < allele_reads_tr or alt_read_count < allele_reads_tr:
            continue

        snps_collection[parsed_line[0]].append((int(parsed_line[1]), ref_read_count, alt_read_count))
    if force_sort:
        chromosomes_order = ChromosomePosition.sorted_chromosomes
    return snps_collection, chromosomes_order, opened_file.name