'''
AssemblyLine: transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2016 Matthew Iyer

@author: mkiyer
'''
import collections


Exon = collections.namedtuple('Exon', ['start', 'end'])


class AssemblyLineError:
    pass


class Strand:
    POS = 0
    NEG = 1
    NA = 2

    FROM_GTF = {'+': POS, '-': NEG, '.': NA}
    TO_GTF = {POS: '+', NEG: '-', NA: '.'}

    @staticmethod
    def from_gtf(s):
        return Strand.FROM_GTF[s]

    @staticmethod
    def to_gtf(s):
        return Strand.TO_GTF(s)
