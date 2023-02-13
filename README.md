# BABACHI: Background Allelic Dosage Bayesian Checkpoint Identification
[![DOI](https://zenodo.org/badge/255952669.svg)](https://zenodo.org/badge/latestdoi/255952669) <br>
BABACHI is a tool for estimation of relative Background Allelic Dosage (BAD) from
non-phased heterozygous SNVs. It estimates BAD directly from enriched sequencing data, where
the precise estimation of allelic copy numbers is not possible.

BAD corresponds to the ratio of Major allele copy number to Minor allele copy number. More details and algorithm description are available [here](https://www.nature.com/articles/s41467-021-23007-0)

## Files format
BABACHI accepts either a BED file with heterozygous SNVs or a standard VCF file sorted by genomic positions (ascending).
- a BED file should begin with the following 8 columns; additional columns are permitted and ignored:
```chromosome, start, end, ID, reference base, alterntive base, reference read count, alternative read count, sample_id```.<br>
Lines starting with # are ignored.

- a VCF file should have GT and AD fields. 

(!) We suggest to use only <b>common SNPs</b> for BAD estimation. User is expected to filter common variants prior to BABACHI usage. See [here](https://www.nature.com/articles/s41467-021-23007-0) for more details.

The output is a BED file with BAD annotations. The file format is described in the [Demo](#demo) section
## System Requirements
### Hardware requirements
`BABACHI` package requires only a standard computer with enough RAM to support in-memory operations.

### Software requirements
#### OS Requirements
The package can be installed on GNU/Linux and OS X platforms from Python Package Index (PyPI) and GitHub.
The package has been tested on the following systems:
+ MacOS Monterey v12.5
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
numba>=0.53.1
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
or
```
pip3 install 'babachi @ git+https://github.com/autosome-ru/BABACHI.git'
```
- `sudo`, if required.
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
    babachi (<file> | --test) [options]
    babachi visualize <file> (-b <badmap>| --badmap <badmap>) [options]
    babachi filter <file> [options]

Arguments:
    <file>            Path to input VCF file. Expected to be sorted by (chr, pos)

    <path>            Path to the file
    <int>             Non negative integer
    <float>           Non negative number
    <states-string>   String of states separated with "," (to provide fraction use "/", e.g. 4/3).
                      Each state must be >= 1
    <samples-string>  Comma-separated sample names or indices
    <prior-string>    Either "uniform" or "geometric"
    <file-or-link>    Path to existing file or link


Arguments:
    -h, --help                              Show help
    --version                               Show version

    -O <path>, --output <path>              Output directory or file path. [default: ./]
    --test                                  Run segmentation on test file

    -v, --verbose                           Write debug messages
    --sample-list <samples-string>          Comma-separated sample names or integer indices to use from input VCF
    --snp-strategy <snp-strategy>           Strategy for the SNPs at the same genomic position (from different samples).
                                            Either add read counts 'ADD' or treat as separate events 'SEP'. [default: SEP]

    -n, --no-filter                         Skip filtering of input file
    --filter-no-rs                          Filter variants without assigned ID in VCF file.
    -f, --force-sort                        Chromosomes in output file will be sorted in numerical order
    -j <int>, --jobs <int>                  Number of jobs to use, parallel by chromosomes [default: 1]
    --chrom-sizes <file-or-link>            File with chromosome sizes (can be a link), default is hg38
    -a <int>, --allele-reads-tr <int>       Allelic reads threshold. Input SNPs will be filtered by ref_read_count >= x and
                                            alt_read_count >= x. Required for correct estimations in underlying statistical model [default: 5]
    -p <string>, --prior <prior-string>     Prior to use. Can be either uniform or geometric [default: uniform]
    -g <float>, --geometric-prior <float>   Coefficient for geometric prior [default: 0.98]
    -s <string>, --states <states-string>   States string [default: 1,2,3,4,5,6]

    -B <float>, --boundary-penalty <float>  Boundary penalty coefficient [default: 4]
    -Z <int>, --min-seg-snps <int>          Only allow segments containing Z or more unique SNPs (IDs/positions) [default: 3]
    -R <int>, --min-seg-bp <int>            Only allow segments containing R or more base pairs [default: 1000]
    -P <int>, --post-segment-filter <int>   Remove segments with less than P unique SNPs (IDs/positions) from output [default: 0]
    -A <int>, --atomic-region-size <int>    Atomic region size in # of SNPs [default: 600]
    -C <int>, --chr-min-snps <int>          Minimum number of SNPs at a chromosome to start segmentation [default: 100]
    -S <int>, --subchr-filter <int>         Exclude subchromosomes with less than C unique SNPs  [default: 3]

Visualization:
    -b <path>, --badmap <path>              BADmap file created with BABACHI
    --visualize                             Perform visualization of SNP-wise AD and BAD for each chromosome.
                                            Will create a directory in output path for the <ext> visualizations
    -z, --zip                               Zip visualizations directory
    -e <ext>, --ext <ext>                   Extension to save visualizations with [default: svg]
```

## Demo
To perform a test run:
```
babachi --test
```
The test run takes approximately 2 minutes on a standard computer.
<br>
The result is a file named `test.bed` that will be created in the working directory (if `-O` option was not provided).
The contents of the `test.bed` file should have the following format:
```
#chr	start	end	BAD	Q1.00	Q1.33	Q1.50	Q2.00	Q2.50	Q3.00	Q4.00	Q5.00	Q6.00	SNP_count	sum_cover
chr1	1	125183196	2	-63.47825919524621	-24.598710473939718	-8.145646624117944	-2.000888343900442e-11	-30.773041699645546	-78.80480783186977	-189.88685134708248	-299.82657588596703	-401.6012195141575	1325	17280
```
Each row represents a single segment with a constant estimated BAD. The most important columns are:
- \#chr:  chromosome
- start: segment start position
- end: segment end position
- BAD: estimated BAD score

Additional columns:
- SNP_count: number of SNPs in the segment
- SNP_ID_count: number of unique SNPs in the segment
- sum_cover: the total read coverage of all SNPs of the segment
- Q<b>X</b>: the logarithmic likelihood of the segment to have BAD = <b>X</b> 

The BABACHI tool is maintained by Sergey Abramov and Alexandr Boytsov.
