'''
AssemblyLine: transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2016 Matthew Iyer

@author: mkiyer
'''
import os
import logging

from assemblyline.lib2.base import AssemblyLineError


class Sample(object):
    REF_ID = 'R'

    def __init__(self, gtf_file, _id):
        self._id = _id
        self.gtf_file = gtf_file

    @staticmethod
    def parse_tsv(filename, header=False, sep='\t'):
        cur_sample_id = 1
        gtf_files = set()
        ids = set()
        samples = []
        with open(filename) as f:
            if header:
                f.next()
            # table rows
            for line in f:
                fields = line.strip().split(sep)
                gtf_file = fields[0]
                if gtf_file in gtf_files:
                    m = "GTF file '%s' is not unique" % gtf_file
                    raise AssemblyLineError(m)
                if not Sample.gtf_valid(gtf_file):
                    m = "GTF file '%s' is not valid" % gtf_file
                    raise AssemblyLineError(m)
                if len(fields) > 1:
                    _id = fields[1]
                    if _id in ids:
                        m = "sample_id '%s' is not unique" % _id
                        raise AssemblyLineError(m)
                else:
                    _id = cur_sample_id
                    cur_sample_id += 1
                samples.append(Sample(gtf_file, _id))
        return samples

    @staticmethod
    def write_tsv(samples, filename, header=True, sep='\t'):
        with open(filename, 'w') as f:
            if header:
                print >>f, sep.join(['gtf', 'sample_id'])
            for s in samples:
                print >>f, sep.join([s.gtf_file, str(s._id)])

    @staticmethod
    def gtf_valid(gtf_file):
        if gtf_file is None:
            logging.error('GTF file %s is None' % (gtf_file))
            return False
        if not os.path.exists(gtf_file):
            logging.error('GTF file %s not found' % (gtf_file))
            return False
        return True
