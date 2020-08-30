# BABACHI: Backgroud Allelic Dosage Bayesian Checkpoint Identification
[![DOI](https://zenodo.org/badge/255952669.svg)](https://zenodo.org/badge/latestdoi/255952669) <br>
BABACHI is a tool for constant Background Allelic Dosage (BAD) genome regions calling from
non-phased heterozygous SNPs. It is aimed at estimation of BAD on low-coverage sequencing data, where
the precise estimation of allelic copy numbers is not possible.

BAD corresponds to the ratio of Major copy number to Minor copy number.

BABACHI takes in a vcf-like .tsv file with heterozygous SNPs sorted by genome positions (ascending).
The input file must contain the following first 7 columns:
chromosome, position, ID, reference base, alternative base, reference read count, alternative read count
All lines, starting with # are ignored.

The output is a .bed file with BAD annotations.

## Installation

```
pip3 install babachi 
```

## Requirements
```
python >= 3.6
```

## Usage
```
babachi <options>...
```
To get full usage description one can execute:
```
babachi --help
```

## Quick start
To perform a test run:
```
babachi --test
```


The BABACHI tool is maintained by Sergey Abramov and Alexandr Boytsov.
