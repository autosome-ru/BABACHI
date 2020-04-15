import gzip
import os
import sys
import json

sys.path.insert(1, "/home/abramov/ASB-Project")
from scripts.HELPERS.paths_for_components import ploidy_path, ploidy_dict_path
from scripts.HELPERS.helpers import make_dict_from_vcf, make_list_from_vcf


def merge_vcfs_add_counts(out_file_name, in_files):
    vcf_dict = dict()
    for file in in_files:
        with gzip.open(file, 'rt') as vcf:
            make_dict_from_vcf(vcf, vcf_dict)

    vcf_keys = list(vcf_dict.keys())
    vcf_keys.sort(key=lambda cords: cords[1])
    vcf_keys.sort(key=lambda cords: cords[0])

    with open(out_file_name, 'w') as out:
        for (chr, pos, ID, REF, ALT) in vcf_keys:
            (R, A) = vcf_dict[(chr, pos, ID, REF, ALT)]
            out.write('\t'.join(map(str, [chr, pos, ID, REF, ALT, R, A])) + '\n')


def merge_vcfs_independent_snps(out_file_name, in_files):
    vcf_list = []
    for file in in_files:
        with gzip.open(file, 'rt') as vcf:
            vcf_list += make_list_from_vcf(vcf)

    vcf_list.sort(key=lambda cords: cords[1])
    vcf_list.sort(key=lambda cords: cords[0])

    with open(out_file_name, 'w') as out:
        for (chr, pos, ID, REF, ALT, R, A) in vcf_list:
            out.write('\t'.join(map(str, [chr, pos, ID, REF, ALT, R, A])) + '\n')


if __name__ == '__main__':
    with open(ploidy_dict_path, 'r') as read_file:
        d = json.loads(read_file.readline())
    key = sys.argv[1]
    mode = 'independent'
    print(key)

    arr = []
    for path in d[key]:
        if os.path.isfile(path):
            arr.append(path)
    if not os.path.isdir(ploidy_path + 'merged_vcfs/'):
        try:
            os.mkdir(ploidy_path + 'merged_vcfs/')
        except:
            pass
    out_file = ploidy_path + 'merged_vcfs/' + key + ".tsv"

    if mode == 'independent':
        merge_vcfs_independent_snps(out_file, arr)
    elif mode == 'add':
        merge_vcfs_add_counts(out_file, arr)
    else:
        raise ValueError(mode)
