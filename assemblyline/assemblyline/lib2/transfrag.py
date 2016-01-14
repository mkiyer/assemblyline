'''
AssemblyLine: transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer

@author: mkiyer
'''
import collections
import logging

from base import Exon
from gtf import GTF, GTFError


class Transfrag(object):
    __slots__ = ('chrom', 'start', 'end', 'strand', '_id', 'sample_id',
                 'expr', 'is_ref', 'exons')

    def __init__(self, chrom, start, end, strand, _id, sample_id, expr,
                 is_ref, exons=None):
        self.chrom = chrom
        self.start = start
        self.end = end
        self.strand = strand
        self._id = _id
        self.sample_id = sample_id
        self.expr = expr
        self.is_ref = is_ref
        self.exons = [] if exons is None else exons

    @property
    def length(self):
        return sum((e.end - e.start) for e in self.exons)

    @staticmethod
    def from_gtf(f):
        '''GTF.Feature object to Transfrag'''
        return Transfrag(f.seqid, f.start, f.end, f.strand,
                         f.attrs[GTF.Attr.TRANSCRIPT_ID],
                         f.attrs[GTF.Attr.SAMPLE_ID],
                         float(f.attrs[GTF.Attr.EXPRESSION]),
                         bool(int(f.attrs[GTF.Attr.REF])))

    @staticmethod
    def parse_gtf(gtf_lines, ignore_ref):
        '''
        returns OrderedDict key is transcript_id value is Transfrag
        '''
        t_dict = collections.OrderedDict()
        for gtf_line in gtf_lines:
            f = GTF.Feature.from_str(gtf_line)
            t_id = f.attrs[GTF.Attr.TRANSCRIPT_ID]
            is_ref = bool(int(f.attrs[GTF.Attr.REF]))

            if is_ref and ignore_ref:
                continue

            if f.feature == 'transcript':
                if t_id in t_dict:
                    raise GTFError("Transcript '%s' duplicate detected" % t_id)
                t = Transfrag.from_gtf(f)
                t_dict[t_id] = t
            elif f.feature == 'exon':
                if t_id not in t_dict:
                    logging.error('Feature: "%s"' % str(f))
                    raise GTFError("Transcript '%s' exon feature appeared in "
                                   "gtf file prior to transcript feature" %
                                   t_id)
                t = t_dict[t_id]
                t.exons.append(Exon(f.start, f.end))
        return t_dict
