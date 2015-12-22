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

from assemblyline.lib2.gtf import sort_gtf, GTF
from assemblyline.lib2.aggregate import add_sample_gtf
from assemblyline.lib2.assemble import Transfrag, Locus

EXIT_ERROR = 1
EXIT_SUCCESS = 0


class AssemblyLineError(Exception):
    pass


class Sample(object):
    REF_ID = 'R'

    def __init__(self, name, gtf_file, bam_file=None, sample_id=None):
        self.name = name
        self.gtf_file = gtf_file
        self.bam_file = bam_file
        self._id = sample_id

    @staticmethod
    def parse_tsv(filename, header=True, sep='\t'):
        cur_sample_id = 1
        samples = []
        sample_names = set()
        with open(filename) as f:
            if header:
                f.next()
            # table rows
            for line in f:
                fields = line.strip().split(sep)
                sample = Sample(*fields)
                if not sample.gtf_valid():
                    m = "Sample '%s' GTF file not found" % sample.name
                    logging.error(m)
                    raise AssemblyLineError(m)
                if sample.name in sample_names:
                    m = "Sample '%s' name is not unique" % sample.name
                    logging.error(m)
                    raise AssemblyLineError(m)
                sample_names.add(sample.name)
                if sample._id is None:
                    sample._id = cur_sample_id
                    cur_sample_id += 1
                samples.append(sample)
        return samples

    @staticmethod
    def write_tsv(samples, filename, header=True, sep='\t'):
        with open(filename, 'w') as f:
            if header:
                print >>f, sep.join(['name', 'gtf', 'bam', 'sample_id'])
            for s in samples:
                print >>f, sep.join([s.name, s.gtf_file, s.bam_file,
                                     str(s._id)])

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


class AssemblyLine(object):
    VERSION = '0.4.0'

    class Args:
        VERBOSE = False
        GUIDED = False
        GTF_EXPR_ATTR = 'FPKM'
        MIN_FRAG_LENGTH = 200
        MAX_FRAG_LENGTH = 400
        FRAC_MAJOR_ISOFORM = 0.01
        MAX_ISOFORMS = 100
        OUTPUT_DIR = 'assemblyline'
        RESUME = False

        @staticmethod
        def create():
            parser = argparse.ArgumentParser()
            parser.add_argument('-v', '--verbose', dest='verbose',
                                action="store_true",
                                default=AssemblyLine.Args.VERBOSE,
                                help='Enabled detailed logging '
                                '(for debugging)')
            parser.add_argument('--gtf-expr-attr',
                                dest='gtf_expr_attr',
                                default=AssemblyLine.Args.GTF_EXPR_ATTR,
                                metavar='ATTR',
                                help='GTF attribute field containing '
                                'expression [default=%(default)s]')
            parser.add_argument('--min-frag-length', dest='min_frag_length',
                                type=int, metavar='N',
                                default=AssemblyLine.Args.MIN_FRAG_LENGTH,
                                help='Length (bp) of smallest valid fragment '
                                'across all experiments [default=%(default)s]')
            parser.add_argument('--max-frag-length', dest='max_frag_length',
                                type=int, metavar='N',
                                default=AssemblyLine.Args.MAX_FRAG_LENGTH,
                                help='Length (bp) of largest valid fragment '
                                'across all experiments [default=%(default)s]')
            parser.add_argument('--ref-gtf', dest='ref_gtf_file',
                                metavar='<gtf_file>',
                                default=None,
                                help='Option reference GTF file of "true" '
                                'validated transcripts to facilitate guided '
                                'assembly and/or noise filtering')
            parser.add_argument('--guided', dest='guided', action='store_true',
                                default=AssemblyLine.Args.GUIDED,
                                help='Enable guided assembly (requires a '
                                'reference GTF to be specified using '
                                '--ref-gtf)')
            parser.add_argument('--frac-major-isoform',
                                dest='frac_major_isoform', type=float,
                                metavar='X',
                                default=AssemblyLine.Args.FRAC_MAJOR_ISOFORM,
                                help='Report transcript isoforms with '
                                'expression fraction >=X (0.0-1.0) relative '
                                'to the major isoform [default=%(default)s]')
            parser.add_argument('--max-isoforms',
                                dest='max_isoforms', type=int, metavar='N',
                                default=AssemblyLine.Args.MAX_ISOFORMS,
                                help='Maximum isoforms to report for each '
                                'gene [default=%(default)s]')
            parser.add_argument("-o", "--output-dir", dest="output_dir",
                                metavar='DIR',
                                default=AssemblyLine.Args.OUTPUT_DIR,
                                help='directory where output files will be '
                                'stored (if already exists then --resume must '
                                'be specified) [default=%(default)s]')
            parser.add_argument('--resume', dest='resume',
                                action='store_true',
                                default=AssemblyLine.Args.RESUME,
                                help='resume a previous run in <dir>')
            parser.add_argument('sample_file', nargs='?')
            return parser

        @staticmethod
        def load(filename):
            return pickle.load(open(filename))

        @staticmethod
        def dump(args, filename):
            pickle.dump(args, open(filename, 'wb'))

        @staticmethod
        def log(args, func=logging.info):
            logging.info('AssemblyLine version %s' % (AssemblyLine.VERSION))
            spacer = '-' * 78
            fmt = '{:<35}{:<35}'
            func(spacer)
            func(fmt.format('verbose logging:', str(args.verbose)))
            func(fmt.format('output directory:', str(args.output_dir)))
            func(fmt.format('reference GTF file:', str(args.ref_gtf_file)))
            func(fmt.format('guided mode:', str(args.guided)))
            func(fmt.format('GTF expression attribute:', args.gtf_expr_attr))
            func(fmt.format('min fragment length:', args.min_frag_length))
            func(fmt.format('max fragment length:', args.max_frag_length))
            func(fmt.format('fraction major isoform:',
                            args.frac_major_isoform))
            func(fmt.format('max isoforms:', args.max_isoforms))
            return

    class Status(object):
        FIELDS = ('create', 'aggregate', 'assemble')

        def __init__(self):
            for f in AssemblyLine.Status.FIELDS:
                setattr(self, f, False)

        def write(self, filename):
            d = dict((f, getattr(self, f)) for f in AssemblyLine.Status.FIELDS)
            json.dump(d, open(filename, 'w'))

        @staticmethod
        def load(filename):
            d = json.load(open(filename))
            self = AssemblyLine.Status()
            for f in AssemblyLine.Status.FIELDS:
                setattr(self, f, d[f])
            return self

    class Results(object):
        TMP_DIR = 'tmp'
        STATUS_FILE = 'status.json'
        ARGS_FILE = 'args.pickle'
        SAMPLE_FILE = 'samples.txt'
        TRANSFRAGS_GTF_FILE = 'transfrags.gtf'
        TRANSFRAGS_FAIL_GTF_FILE = 'transfrags.fail.gtf'
        AGGREGATE_STATS_FILE = 'aggregate_stats.txt'

        def __init__(self, output_dir):
            self.output_dir = output_dir
            self.tmp_dir = \
                os.path.join(output_dir, AssemblyLine.Results.TMP_DIR)
            self.args_file = \
                os.path.join(output_dir, AssemblyLine.Results.ARGS_FILE)
            self.status_file = \
                os.path.join(output_dir, AssemblyLine.Results.STATUS_FILE)
            self.sample_file = \
                os.path.join(output_dir, AssemblyLine.Results.SAMPLE_FILE)
            self.transfrags_gtf_file = \
                os.path.join(output_dir,
                             AssemblyLine.Results.TRANSFRAGS_GTF_FILE)
            self.transfrags_fail_gtf_file = \
                os.path.join(output_dir,
                             AssemblyLine.Results.TRANSFRAGS_FAIL_GTF_FILE)
            self.aggregate_stats_file = \
                os.path.join(output_dir,
                             AssemblyLine.Results.AGGREGATE_STATS_FILE)

    @staticmethod
    def _parse_args():
        parser = AssemblyLine.Args.create()
        args = parser.parse_args()
        # check if we are trying to resume a previous run
        if args.resume:
            if not os.path.exists(args.output_dir):
                parser.error("Output directory '%s' does not exist" %
                             args.output_dir)
            # check that basic files exist in directory
            results = AssemblyLine.Results(args.output_dir)
            can_resume = os.path.exists(results.args_file)
            can_resume = can_resume and os.path.exists(results.status_file)
            can_resume = can_resume and os.path.exists(results.sample_file)
            if not can_resume:
                parser.error("Output directory '%s' not valid" %
                             args.output_dir)
        else:
            if os.path.exists(args.output_dir):
                parser.error("Output directory '%s' already exists" %
                             args.output_dir)
            if not os.path.exists(args.sample_file):
                parser.error("sample file %s not found" % (args.sample_file))
            args.sample_file = os.path.abspath(args.sample_file)

            if args.min_frag_length <= 0:
                parser.error("min_transcript_length <= 0")
            if args.max_frag_length <= 0:
                parser.error("max_frag_length <= 0")
            if (args.frac_major_isoform <= 0) or (args.frac_major_isoform > 1):
                parser.error("frac_major_isoform out of range (0.0-1.0)")
            if (args.max_isoforms < 1):
                parser.error("max_isoforms <= 0")

            if args.ref_gtf_file is not None:
                if not os.path.exists(args.ref_gtf_file):
                    parser.error("reference GTF file %s not found" %
                                 (args.ref_gtf_file))
            elif args.guided:
                parser.error('Guided assembly mode (--guided) requires '
                             'reference GTF (--ref-gtf)')
        return args

    @staticmethod
    def create():
        self = AssemblyLine()
        # parse command line args
        args = AssemblyLine._parse_args()
        if args.resume:
            self.results = AssemblyLine.Results(args.output_dir)
            self.status = AssemblyLine.Status.load(self.results.status_file)
            self.samples = Sample.parse_tsv(self.results.sample_file)
            self.args = AssemblyLine.Args.load(self.results.args_file)
        else:
            self.args = args
            self.results = AssemblyLine.Results(args.output_dir)
            self.status = AssemblyLine.Status()
            self.samples = Sample.parse_tsv(self.args.sample_file)

        # setup logging
        if args.verbose:
            level = logging.DEBUG
        else:
            level = logging.INFO
        logging.basicConfig(level=level,
                            format="%(asctime)s - %(levelname)s - %(message)s")

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
        if not args.resume:
            AssemblyLine.Args.dump(self.args, self.results.args_file)
            Sample.write_tsv(self.samples, self.results.sample_file)
            # update status and write to file
            self.status.create = True
            self.status.write(self.results.status_file)
        return self

    def aggregate(self):
        '''
        Aggregate/merge individual sample GTF files
        '''
        r = self.results
        a = self.args
        samples = self.samples

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
        if retcode != 0:
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
        a = self.args
        r = self.results

        # setup bedgraph output files
        file_prefix = os.path.join(a.output_dir, 'loci.unresolved')
        raw_bgfilehd = Locus.open_bedgraph(file_prefix)
        file_prefix = os.path.join(a.output_dir, 'loci.resolved')
        resolved_bgfilehd = Locus.open_bedgraph(file_prefix)

        # parse gtf file
        ignore_ref = not a.guided
        for interval, gtf_lines in GTF.parse_loci(open(r.transfrags_gtf_file)):
            chrom, start, end = interval
            t_dict = Transfrag.parse_gtf(gtf_lines, ignore_ref)
            locus = Locus.create(t_dict.values())
            logging.debug('Locus %s:%d-%d: '
                          '%d transfrags (+: %d, -: %d, .: %d)' %
                          (chrom, start, end, len(t_dict),
                           len(locus.strand_transfrags['+']),
                           len(locus.strand_transfrags['-']),
                           len(locus.strand_transfrags['.'])))

            # write bedgraph files for expression/recurrence data
            locus.write_bedgraph(raw_bgfilehd)

            # resolve unstranded transcripts
            num_resolved = locus.predict_unknown_strands()
            if num_resolved > 0:
                logging.debug('Locus %s:%d-%d: %d '
                              'resolved (+: %d, -: %d, .: %d)' %
                              (chrom, start, end, num_resolved,
                               len(locus.strand_transfrags['+']),
                               len(locus.strand_transfrags['-']),
                               len(locus.strand_transfrags['.'])))

            # write bedgraph files after strand resolved
            locus.write_bedgraph(resolved_bgfilehd)

        # close bedgraph files
        Locus.close_bedgraph(raw_bgfilehd)
        Locus.close_bedgraph(resolved_bgfilehd)

        # update status and write to file
        # self.status.assemble = True
        # self.status.write(self.results.status_file)
        # return EXIT_SUCCESS

def main():
    # instantiate from command line
    A = AssemblyLine.create()
    # log command line arguments
    AssemblyLine.Args.log(A.args)
    logging.debug('Samples: %d' % (len(A.samples)))
    # parse gtf files
    msg = 'Aggregating GTF files'
    if A.status.aggregate:
        logging.info('[SKIPPING] %s' % msg)
    else:
        logging.info(msg)
        A.aggregate()
    #
    A.assemble()

    return 0


if __name__ == '__main__':
    sys.exit(main())
