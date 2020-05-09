# BAD Segmentation

BAD Segmentation is a tool for constant Background Allelic Dosage (BAD) genome regions calling from
non-phased heterozygous SNPs. It is aimed at estimation of BAD on low-coverage sequencing data, where
the precise estimation of allelic copy numbers is not possible.

BAD corresponds to the ratio of Major copy number to Minor copy number.

BAD Segmentation takes in a vcf-like .tsv file with heterozygous SNPs sorted by genome positions (ascending).
The input file must contain the following first 7 columns:
chromosome, position, ID, reference base, alternative base, reference read count, alternative read count
All lines, starting with # are ignored.

The output is a .bed file with BAD annotations.

## Installation

//to be uploaded on PyPI
temporary:
``` 
pip3 install git+ssh://git@github.com/autosome-ru/BAD_segmentation
```

## Requirements
```
python >= 3.6
```

## Usage
```
segmentation <options>...
python -m segmentation <options>...
```
To get full usage description one can execute:
```
segmentation --help
```

## Quick start
To perform a test run:
```
segmentation --test
```


The BAD Segmentation tool is maintained by Sergey Abramov and Alexandr Boytsov.
