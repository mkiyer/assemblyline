'''
AssemblyLine: transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer

@author: mkiyer
'''
import os
import sys
import argparse
import logging
import json
import pickle
import collections

from assemblyline.lib2.gtf import sort_gtf, GTF
from assemblyline.lib2.aggregate import add_sample_gtf

EXIT_ERROR = 1
EXIT_SUCCESS = 0


class Sample(object):
    REF_ID = 'Ref'

    def __init__(self, name, gtf_file, bam_file=None):
        self.name = name
        self.gtf_file = gtf_file
        self.bam_file = bam_file

    @staticmethod
    def parse_txt_file(filename, header=True, sep='\t'):
        with open(filename) as f:
            if header:
                f.next()
            # table rows
            for line in f:
                fields = line.strip().split(sep)
                yield Sample(*fields)

    def gtf_valid(self):
        if self.gtf_file is None:
            logging.error('Sample %s GTF file is None' % (self.name))
            return False
        if not os.path.exists(self.gtf_file):
            logging.error('Sample %s GTF file not found' % (self.name))
            return False
        return True

    def bam_valid(self):
        if self.bam_file is None:
            logging.error('Sample %s BAM file is None' % (self.name))
            return False
        if not os.path.exists(self.bam_file):
            logging.error('Sample %s BAM file not found' % (self.name))
            return False
        return True


class Assembler(object):
    VERSION = '0.4.0'

    class Default:
        VERBOSE = False
        GUIDED = False
        GTF_EXPR_ATTR = 'FPKM'
        MIN_FRAG_LENGTH = 200
        MAX_FRAG_LENGTH = 400
        FRAC_MAJOR_ISOFORM = 0.01
        MAX_ISOFORMS = 100
        OUTPUT_DIR = 'assemblyline'

    class Status(object):
        FIELDS = ('create', 'aggregate')

        def __init__(self):
            for f in Assembler.Status.FIELDS:
                setattr(self, f, False)

        def write(self, filename):
            d = dict((f, getattr(self, f)) for f in Assembler.Status.FIELDS)
            json.dump(d, open(filename, 'w'))

        @staticmethod
        def load(filename):
            d = json.load(open(filename))
            self = Assembler.Status()
            for f in Assembler.Status.FIELDS:
                setattr(self, f, d[f])

    class Results(object):
        TMP_DIR = 'tmp'
        STATUS_FILE = 'status.json'
        ARGS_FILE = 'args.pickle'
        SAMPLE_ID_MAP_FILE = 'sample_id_map.txt'
        TRANSFRAGS_GTF_FILE = 'transfrags.gtf'
        TRANSFRAGS_FAIL_GTF_FILE = 'transfrags.fail.gtf'
        AGGREGATE_STATS_FILE = 'aggregate_stats.txt'

        def __init__(self, output_dir):
            self.output_dir = output_dir
            self.tmp_dir = \
                os.path.join(output_dir, Assembler.Results.TMP_DIR)
            self.args_file = \
                os.path.join(output_dir, Assembler.Results.ARGS_FILE)
            self.status_file = \
                os.path.join(output_dir, Assembler.Results.STATUS_FILE)
            self.sample_id_map_file = \
                os.path.join(output_dir, Assembler.Results.SAMPLE_ID_MAP_FILE)
            self.transfrags_gtf_file = \
                os.path.join(output_dir, Assembler.Results.TRANSFRAGS_GTF_FILE)
            self.transfrags_fail_gtf_file = \
                os.path.join(output_dir,
                             Assembler.Results.TRANSFRAGS_FAIL_GTF_FILE)
            self.aggregate_stats_file = \
                os.path.join(output_dir,
                             Assembler.Results.AGGREGATE_STATS_FILE)

    @staticmethod
    def _create_arg_parser():
        parser = argparse.ArgumentParser()
        parser.add_argument('-v', '--verbose', dest='verbose',
                            action="store_true",
                            default=Assembler.Default.VERBOSE,
                            help='Enabled detailed logging (for debugging)')
        parser.add_argument('--gtf-expr-attr',
                            dest='gtf_expr_attr',
                            default=Assembler.Default.GTF_EXPR_ATTR,
                            metavar='ATTR',
                            help='GTF attribute field containing expression '
                            ' [default=%(default)s]')
        parser.add_argument('--min-frag-length', dest='min_frag_length',
                            type=int, metavar='N',
                            default=Assembler.Default.MIN_FRAG_LENGTH,
                            help='Length (bp) of smallest valid fragment '
                            'across all experiments [default=%(default)s]')
        parser.add_argument('--max-frag-length', dest='max_frag_length',
                            type=int, metavar='N',
                            default=Assembler.Default.MAX_FRAG_LENGTH,
                            help='Length (bp) of largest valid fragment '
                            'across all experiments [default=%(default)s]')
        parser.add_argument('--ref-gtf', dest='ref_gtf_file',
                            metavar='<gtf_file>',
                            default=None,
                            help='Option reference GTF file of "true" '
                            'validated transcripts to facilitate guided '
                            'assembly and/or noise filtering')
        parser.add_argument('--guided', dest='guided', action='store_true',
                            default=Assembler.Default.GUIDED,
                            help='Enable guided assembly (requires a '
                            'reference GTF to be specified using --ref-gtf)')
        parser.add_argument('--frac-major-isoform',
                            dest='frac_major_isoform', type=float,
                            metavar='X',
                            default=Assembler.Default.FRAC_MAJOR_ISOFORM,
                            help='Report transcript isoforms with expression '
                            'fraction >=X (0.0-1.0) relative to the major '
                            'isoform [default=%(default)s]')
        parser.add_argument('--max-isoforms',
                            dest='max_isoforms', type=int, metavar='N',
                            default=Assembler.Default.MAX_ISOFORMS,
                            help='Maximum isoforms to report for each gene '
                            '[default=%(default)s]')
        parser.add_argument("-o", "--output-dir", dest="output_dir",
                            default=Assembler.Default.OUTPUT_DIR,
                            help="directory where output files will be stored "
                            "(created if it does not exist) "
                            "[default=%(default)s]")
        parser.add_argument('sample_file')
        return parser

    @staticmethod
    def _parse_args():
        parser = Assembler._create_arg_parser()
        args = parser.parse_args()
        # check for valid parameters
        if args.min_frag_length <= 0:
            parser.error("min_transcript_length <= 0")
        if args.max_frag_length <= 0:
            parser.error("max_frag_length <= 0")
        if (args.frac_major_isoform <= 0) or (args.frac_major_isoform > 1):
            parser.error("frac_major_isoform out of range (0.0-1.0)")
        if (args.max_isoforms < 1):
            parser.error("max_isoforms <= 0")
        # check command line parameters
        if not os.path.exists(args.sample_file):
            parser.error("sample file %s not found" % (args.sample_file))
        if args.ref_gtf_file is not None:
            if not os.path.exists(args.ref_gtf_file):
                parser.error("reference GTF file %s not found" %
                             (args.ref_gtf_file))
        elif args.guided:
            parser.error('Guided assembly mode (--guided) requires reference '
                         'GTF (--ref-gtf)')
        if os.path.exists(args.output_dir):
            parser.error("Output directory '%s' already exists" %
                         args.output_dir)
        return args

    @staticmethod
    def create():
        # parse command line args
        args = Assembler._parse_args()
        # setup logging
        if args.verbose:
            level = logging.DEBUG
        else:
            level = logging.INFO
        logging.basicConfig(level=level,
                            format="%(asctime)s - %(levelname)s - %(message)s")
        # create assembler instance
        self = Assembler()
        self.args = args
        self.results = Assembler.Results(args.output_dir)
        self.status = Assembler.Status()
        # create output directories
        results = self.results
        if not os.path.exists(results.output_dir):
            logging.debug("Creating output directory '%s'" %
                          (results.output_dir))
            os.makedirs(results.output_dir)
        if not os.path.exists(results.tmp_dir):
            logging.debug("Creating tmp directory '%s'" % (results.tmp_dir))
            os.makedirs(results.tmp_dir)
        # write assembler configuration to file
        self.write_args()
        # update status and write to file
        self.status.create = True
        self.status.write(self.results.status_file)
        return self

    def load_args(self):
        self.args = pickle.load(open(self.results.args_file))

    def write_args(self):
        pickle.dump(self.args, open(self.results.args_file, 'wb'))

    def parse_samples(self):
        cur_sample_id = 1
        samples = []
        sample_names = set()
        sample_id_map_fileh = open(self.results.sample_id_map_file, 'w')
        print >>sample_id_map_fileh, 'id\tname'
        for sample in Sample.parse_txt_file(self.args.sample_file):
            if not sample.gtf_valid():
                logging.error("Sample '%s' GTF file not found" % sample.name)
                return EXIT_ERROR
            if sample.name in sample_names:
                logging.error("Sample '%s' name is not unique" % sample.name)
                return EXIT_ERROR
            sample_names.add(sample.name)
            print >>sample_id_map_fileh, '%d\t%s' % \
                (cur_sample_id, sample.name)
            sample._id = cur_sample_id
            samples.append(sample)
            cur_sample_id += 1
        sample_id_map_fileh.close()
        return samples

    def log_config(self, logging_func=logging.info):
        logging.info('AssemblyLine version %s' % (Assembler.VERSION))
        logging.info('-' * 78)
        logging.info('verbose logging:          %s' %
                     str(self.args.verbose))
        logging.info('output directory:         %s' %
                     self.args.output_dir)
        logging.info('reference GTF file:       %s' %
                     self.args.ref_gtf_file)
        logging.info('guided:                   %s' %
                     str(self.args.guided))
        logging.info('sample file:              %s' %
                     self.args.sample_file)
        logging.info('GTF expression attribute: %s' %
                     self.args.gtf_expr_attr)
        logging.info('min fragment length:      %d' %
                     self.args.min_frag_length)
        logging.info('max fragment length:      %d' %
                     self.args.max_frag_length)
        logging.info('fraction major isoform:   %f' %
                     self.args.frac_major_isoform)
        logging.info('max isoforms:             %d' %
                     self.args.max_isoforms)
        logging.info('-' * 78)

    def aggregate(self, samples):
        '''
        Aggregate/merge individual sample GTF files
        '''
        r = self.results
        a = self.args

        # setup output files
        tmp_file = os.path.join(r.tmp_dir, 'transcripts.unsorted.gtf')
        tmp_fileh = open(tmp_file, 'w')
        stats_fileh = open(r.aggregate_stats_file, 'w')
        # stats file has header
        fields = ['sample', 'expr_tot', 'expr_quantiles', 'length_quantiles']
        print >>stats_fileh, '\t'.join(fields)

        # aggregate ref gtf
        if a.ref_gtf_file is not None:
            sample = Sample('ref', a.ref_gtf_file)
            sample._id = Sample.REF_ID
            logging.debug('Sample: %s' % sample._id)
            add_sample_gtf(sample, a.gtf_expr_attr, tmp_fileh, stats_fileh,
                           is_ref=True)
        # aggregate sample gtfs
        for sample in samples:
            logging.debug('Sample: %s' % sample._id)
            add_sample_gtf(sample, a.gtf_expr_attr, tmp_fileh, stats_fileh)
        tmp_fileh.close()
        stats_fileh.close()

        logging.info("Sorting GTF")
        retcode = sort_gtf(tmp_file, r.transfrags_gtf_file,
                           tmp_dir=r.tmp_dir)
        if retcode != EXIT_SUCCESS:
            logging.error("Error sorting GTF")
            if os.path.exists(r.transfrags_gtf_file):
                os.remove(r.transfrags_gtf_file)
            return EXIT_ERROR
        os.remove(tmp_file)

        # update status and write to file
        self.status.aggregate = True
        self.status.write(self.results.status_file)
        return EXIT_SUCCESS

    def assemble(self):
        # parse gtf file
        for lines in GTF.parse_loci(self.results.transfrags_gtf_file):

            t_dict = collections.OrderedDict()
            exon_boundaries = []
            for line in lines:
                f = GTF.Feature.from_str(line)
                t_id = t_dict[f.attrs[GTF.Attr.TRANSCRIPT_ID]]
                if f.feature == 'transcript':
                    if t_id in t_dict:
                        # error duplicate transcript ids
                        pass
                    else:
                        t_dict[t_id] = f
                    continue
                elif f.feature == 'exon':
                    continue

            pass


def main():
    # instantiate assembler from command line
    assembler = Assembler.create()
    # log command line arguments
    assembler.log_config()
    # parse sample file
    logging.info('Parsing sample file')
    samples = assembler.parse_samples()
    logging.debug('Samples: %d' % (len(samples)))
    # parse gtf files
    logging.info('Aggregating GTF files')
    assembler.aggregate(samples)

    return 0


if __name__ == '__main__':
    sys.exit(main())
