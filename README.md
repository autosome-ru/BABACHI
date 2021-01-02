# BABACHI: Backgroud Allelic Dosage Bayesian Checkpoint Identification
[![DOI](https://zenodo.org/badge/255952669.svg)](https://zenodo.org/badge/latestdoi/255952669) <br>
BABACHI is a tool for Background Allelic Dosage (BAD) genomic regions calling from
non-phased heterozygous SNVs. It is aimed at estimation of BAD on low-coverage sequencing data, where
the precise estimation of allelic copy numbers is not possible.

BAD corresponds to the ratio of Major copy number to Minor copy number.

BABACHI takes in a vcf-like .tsv file with heterozygous SNVs sorted by genome positions (ascending).
The input file must contain the following first 7 columns:
chromosome, position, ID, reference base, alternative base, reference read count, alternative read count
All lines, starting with # are ignored.

The output is a .bed file with BAD annotations.
## System Requirements
### Hardware requirements
`BABABCHI` package requires only a standard computer with enough RAM to support the in-memory operations.

### Software requirements
#### OS Requirements
The package can be installed on all major platforms (e.g. BSD, GNU/Linux, OS X, Windows) from Python Package Index (PyPI) and GitHub.
The package has been tested on the following systems:
+ Windows: Windows 10 
+ Linux: Ubuntu 18.04
#### Python Dependencies
`BABACHI` mainly depends on the following Python 3 packages:
```
docopt>=0.6.2
numpy>=1.18.0
schema>=0.7.2
contextlib2>=0.5.5
pandas>=1.0.4
matplotlib>=3.2.1
seaborn>=0.10.1
```
## Installation
### Install from PyPi
```
pip3 install babachi 
```
### Install from Github
```
git clone https://github.com/autosome-ru/BABACHI
cd BABACHI
python3 setup.py install
```
- `sudo`, if required
The package should take less than 1 minute to install.

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
This will produce the following message:
```
Usage:
    babachi <file> [-O <path> |--output <path>] [-q | --quiet] [--allele_reads_tr <int>] [--force-sort] [--visualize] [--boundary-penalty <float>] [--states <string>]
    babachi (--test) [-O <path> |--output <path>] [-q | --quiet] [--allele_reads_tr <int>] [--force-sort] [--visualize] [--boundary-penalty <float>]
    babachi visualize <file> (-b <badmap>| --badmap <badmap>) [-q | --quiet] [--allele_reads_tr <int>]
    babachi -h | --help

Arguments:
    <file>            Path to input file in tsv format with columns:
                      chr pos ID ref_base alt_base ref_read_count alt_read_count.
    <badmap>          Path to badmap .bed format file
    <int>             Non negative integer
    <float>           Non negative number
    <states_string>   String of states separated with "," (to provide fraction use "/", e.g. 4/3). Each state must be >= 1


Options:
    -h, --help                      Show help.
    -q, --quiet                     Less log messages during work time.
    -b <badmap>, --badmap <badmap>  Input badmap file
    -O <path>, --output <path>      Output directory or file path. [default: ./]
    --allele_reads_tr <int>         Allelic reads threshold. Input SNPs will be filtered by ref_read_count >= x and
                                    alt_read_count >= x. [default: 5]
    --force-sort                    Do chromosomes need to be sorted
    --visualize                     Perform visualization of SNP-wise AD and BAD for each chromosome.
                                    Will create a directory in output path for the .svg visualizations.
    --boundary-penalty <float>      Boundary penalty coefficient [default: 9]
    --states <states_string>        States string [default: 1,2,3,4,5,6,1.5,2.5]
    --test                          Run segmentation on test file
```

## Demo
To perform a test run:
```
babachi --test
```
The test run takes approximately 2 minutes on a standard computer.
<br>
The result is a file named `test.bed` that will be produced in the root directory of the project (if `-O` option is not used).
The contents of the `test.bed` file should start as follows:
```
#chr	start	end	BAD	Q1.00	Q1.33	Q1.50	Q2.00	Q2.50	Q3.00	Q4.00	Q5.00	Q6.00	SNP_count	sum_cover
chr1	1	125183196	2	-63.47825919524621	-24.598710473939718	-8.145646624117944	-2.000888343900442e-11	-30.773041699645546	-78.80480783186977	-189.88685134708248	-299.82657588596703	-401.6012195141575	1325	17280
```
Each row represents a single segment with a constant estimated BAD. The columns are as follows:
- #chr:  chromosome
- start: segment start position
- end: segment end position
- BAD: estimated BAD
- Q<b>X</b>: the logarithmic likelyhood of the segment to have BAD = <b>X</b>
- SNP_count: number of SNPs of the segment
- sum_cover: the total read coverage of all SNPs of the segment
<br>
The BABACHI tool is maintained by Sergey Abramov and Alexandr Boytsov.
