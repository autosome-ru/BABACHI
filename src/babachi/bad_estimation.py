"""
Usage:
    babachi (<file> | --test) [options]
    babachi visualize <file> (-b <badmap>| --badmap <badmap>) [options]
    babachi filter <file> [options]

Arguments:
    <file>            Path to input VCF file. Expected to be sorted by (chr, pos)

    <path>            Path to the file
    <int>             Non negative integer
    <float>           Non negative number
    <states-string>   Allowed states separated with "," (to provide fraction use "/", e.g. 4/3).
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
    --snp-strategy <snp-strategy>           Strategy for the SNPs on the same position (from different samples).
                                            Either add read counts 'ADD' or treat as separate events 'SEP'. [default: SEP]

    -n, --no-filter                         Skip filtering of input file
    --filter-no-rs                          Filter variants without assigned ID in VCF file.
    -f, --force-sort                        Chromosomes in output file will be sorted in numerical order
    -j <int>, --jobs <int>                  Number of jobs to use, parallel by chromosomes [default: 1]
    --chrom-sizes <file-or-link>            File with chromosome sizes (can be a link), default is hg38
    -a <int>, --allele-reads-tr <int>       Allelic reads threshold. Input SNPs will be filtered by ref_read_count >= x and
                                            alt_read_count >= x. [default: 5]
    -p <string>, --prior <prior-string>     Prior to use. Can be either uniform or geometric [default: uniform]
    -g <float>, --geometric-prior <float>   Coefficient for geometric prior [default: 0.98]
    -s <string>, --states <states-string>   States string [default: 1,2,3,4,5,6]

    -B <float>, --boundary-penalty <float>  Boundary penalty coefficient [default: 4]
    -Z <int>, --min-seg-snps <int>          Only allow segments containing Z or more unique SNPs (IDs/positions) [default: 3]
    -R <int>, --min-seg-bp <int>            Only allow segments containing R or more base pairs [default: 1000]
    -P <int>, --post-segment-filter <int>   Remove segments with less than P unique SNPs (IDs/positions) from output [default: 0]
    -A <int>, --atomic-region-size <int>    Atomic region size in SNPs [default: 600]
    -C <int>, --chr-min-snps <int>          Minimum number of SNPs on a chromosome to start segmentation [default: 100]
    -S <int>, --subchr-filter <int>         Exclude subchromosomes with less than C unique SNPs  [default: 3]

Visualization:
    -b <path>, --badmap <path>              BADmap file created with BABACHI
    --visualize                             Perform visualization of SNP-wise AD and BAD for each chromosome.
                                            Will create a directory in output path for the <ext> visualizations
    -z, --zip                               Zip visualizations directory
    -e <ext>, --ext <ext>                   Extension to save visualizations with [default: svg]
"""
import logging
import validators
from schema import Schema, And, Use, SchemaError, Const, Or
import time

from docopt import docopt # TODO Use click instead of docopt
import os
from babachi.visualize_segmentation import BabachiVisualizer
from babachi.helpers import read_url_file, read_snps_file, df_header, make_file_path_from_dir, craft_prior
from babachi.validators import check_states, check_samples
from babachi.version import __version__
from babachi.logging import root_logger, set_logger_config
from babachi.segmentation import GenomeSegmentator
from babachi.models import GenomeSNPsHandler


def segmentation_start():
    args = docopt(__doc__, version=__version__)
    if args['--test']:
        args['<file>'] = os.path.join(os.path.dirname(__file__), 'tests', 'test.bed')

    schema = Schema({
        '<file>': And(
            Const(os.path.exists, error='Input file should exist'),
        ),
        '--boundary-penalty': And(
            Use(float),
            Const(lambda x: x >= 0), error='Boundary penalty coefficient should be non negative number'
        ),
        '--min-seg-snps': And(
            Use(int),
            Const(lambda x: x >= 0), error='Min SNPs in segment parameter coefficient should be non negative integer'
        ),
        '--min-seg-bp': And(
            Use(int),
            Const(lambda x: x >= 0), error='Min segment length coefficient should be non negative integer'
        ),
        '--prior': Const(
            lambda x: x in ('uniform', 'geometric'),
            error='Invalid prior string. Must be "uniform" or "geometric"'
        ),
        '--snp-strategy': Const(
            lambda x: x in ('ADD', 'SEP'),
            error='SNP strategy should be either ADD or SEP'
        ),
        '--chrom-sizes': Or(
            Const(lambda x: x is None),
            And(
                Const(validators.url, error='Not a valid URL'),
                Use(read_url_file)
            ),
            Const(os.path.exists, error='File should exist')
        ),
        '--badmap': Or(
            Const(lambda x: x is None and not args['visualize']),
            And(
                Const(os.path.exists, error='Badmap file should exist'),
                Const(lambda x: os.access(x, os.R_OK), error='No read permission for badmap file')
            )),
        '--output': Or(
            And(
                Const(os.path.exists),
                Const(lambda x: os.access(x, os.W_OK), error='No write permissions')
            ),
            And(
                Const(lambda x: not os.path.exists(x)),
                Const(lambda x: os.access(os.path.dirname(x) if os.path.dirname(x) != '' else '.', os.W_OK),
                      error='No write permissions')
            )
        ),
        '--states': Use(
            check_states, error='''Incorrect value for --states.
            Must be "," separated list of numbers or fractions in the form "x/y", each >= 1'''
        ),
        '--sample-list': Or(
            Use(
                check_samples, error='Invalid sample list'
            ),
            Const(lambda x: x is None),
        ),
        '--allele-reads-tr': And(
            Use(int),
            Const(lambda x: x >= 0), error='Allelic reads threshold must be a non negative integer'
        ),
        '--post-segment-filter': And(
            Use(int),
            Const(lambda x: x >= 0), error='Segments length filter (in SNPs) must be a non negative integer'
        ),
        '--subchr-filter': And(
            Use(int),
            Const(lambda x: x >= 0), error='Subchromosome post filter (in SNPs) must be a non negative integer'
        ),
        '--chr-min-snps': And(
            Use(int),
            Const(lambda x: x >= 0), error='Subchromosome pre filter (in SNPs) must be a non negative integer'
        ),
        '--atomic-region-size': And(
            Use(int),
            Const(lambda x: x >= 200), error='Atomic region size (in SNPs) must be a non negative integer'
        ),
        '--jobs': Use(int, error='Number of jobs should be positive integer'),
        '--geometric-prior': And(
            Use(float, error='Geometric prior coefficient should be a number'),
            Const(lambda x: 0 <= x <= 1, error='Coefficient should be between 0 and 1')
        ),
        '--ext': str,
        str: bool
    })
    try:
        args = schema.validate(args)
    except SchemaError as e:
        print(__doc__)
        exit('Error: {}'.format(e))

    verbose = args['--verbose']
    if verbose:
        level = logging.DEBUG
    else:
        level = logging.INFO

    set_logger_config(root_logger, level)
    full_name = args['<file>']
    file_name, ext = os.path.splitext(os.path.basename(full_name))

    try:
        snps, chrom_wrapper = read_snps_file(
            file_path=full_name,
            chrom_sizes=args['--chrom-sizes'],
            snp_strategy=args['--snp-strategy'],
            allele_reads_tr=args['--allele-reads-tr'],
            force_sort=args['--force-sort'],
            to_filter=not args['--no-filter'] or args['filter'],
            filter_no_rs=args['--filter-no-rs']
        )
    except Exception as e:
        raise ValueError("Can not read the input file", *e.args)
    if args['filter']:
        snps = snps[df_header]
        snps['#chr'] = snps['chr']
        snps[['#chr', 'start', 'end', 'ID', 'ref', 'alt', 'ref_counts', 'alt_counts', 'sample_id']].to_csv(
            make_file_path_from_dir(args['--output'], file_name, 'snps.bed'),
            sep='\t', index=False)
        root_logger.info('Succesfully converted to BED')
        exit(0)
        return
    chrom_order = snps['chr'].unique()
    chrom_sizes = {x: chrom_wrapper.chromosomes[x] for x in chrom_order}
    snps_collection = GenomeSNPsHandler.from_df(data=snps, chrom_wrapper=chrom_wrapper)

    if not args['visualize']:
        badmap_file_path = make_file_path_from_dir(args['--output'], file_name)
        mode = 'corrected'
        t = time.perf_counter()
        genome_segmentator = GenomeSegmentator(
            snps_collection=snps_collection,
            chrom_sizes=chrom_sizes,
            segmentation_mode=mode,
            states=args['--states'],
            b_penalty=args['--boundary-penalty'],
            allele_reads_tr=args['--allele-reads-tr'],
            min_seg_snps=args['--min-seg-snps'],
            min_seg_bp=args['--min-seg-bp'],
            post_seg_filter=args['--post-segment-filter'],
            atomic_region_size=args['--atomic-region-size'],
            chr_filter=args['--chr-min-snps'],
            subchr_filter=args['--subchr-filter'],
            prior=craft_prior(
                args['--states'],
                args['--prior'],
                args['--geometric-prior']
            ),
            jobs=args['--jobs'],
            logger=root_logger,
            logger_level=level,  # workaround for mp logging,
        )
        try:
            bad_segments = genome_segmentator.estimate_BAD()
            genome_segmentator.write_BAD(bad_segments, badmap_file_path)
        except Exception as e:
            raise e
        root_logger.debug('Total time: {} s'.format(time.perf_counter() - t))
    else:
        badmap_file_path = args['--badmap']
    if args['--visualize'] or args['visualize']:
        visualizer = BabachiVisualizer(chromosomes_wrapper=chrom_wrapper)
        visualizer.init_from_snps_collection(
            snps_collection=snps_collection,
            to_zip=args['--zip'],
            ext=args['--ext'],
            BAD_file=badmap_file_path
        )
