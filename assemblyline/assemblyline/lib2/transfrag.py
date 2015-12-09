'''
AssemblyLine: transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer

@author: mkiyer
'''
import collections



class Transfrag(object):
    __slots__ = ('chrom', 'start', 'end', 'strand', 'expr', 'exons')

    @staticmethod
    def from_gtf_line(line):
        f = GTF.Feature.from_str(line)
        self = Transfrag()
        self.chrom = f.seqid
        self.start = f.start
        self.end = f.end
        self.strand = f.strand
        self.expr = f.attrs[GTF.Attr.EXPRESSION]


        pass




    def __str__(self):
        return ("<%s(chrom='%s', start='%d', end='%d', strand='%s', "
                "score='%s' exons='%s', attrs='%s'" %
                (self.__class__.__name__, self.chrom, self.start, self.end,
                 strand_int_to_str(self.strand), str(self.score), self.exons,
                 self.attrs))

    @property
    def length(self):
        return sum((e.end - e.start) for e in self.exons)

    def iterintrons(self):
        #e1 = self.exons[0]
        #for e2 in self.exons[1:]:
        #    yield (e1.end,e2.start)
        #    e1 = e2
        e1 = self.exons[0]
        for j in xrange(1, len(self.exons)):
            e2 = self.exons[j]
            yield e1.end, e2.start
            e1 = e2

    def introns(self):
        return list(self.iterintrons())

    def to_bed12(self):
        block_sizes = []
        block_starts = []
        for e0, e1 in self.exons:
            block_starts.append(e0 - self.tx_start)
            block_sizes.append(e1 - e0)
        # write
        s = '\t'.join([self.chrom,
                       str(self.start),
                       str(self.end),
                       str(self.attrs["transcript_id"]),
                       '0',
                       strand_int_to_str(self.strand),
                       str(self.start),
                       str(self.start),
                       '0',
                       str(len(self.exons)),
                       ','.join(map(str,block_sizes)) + ',',
                       ','.join(map(str,block_starts)) + ','])
        return s

    def to_gtf_features(self, source=None, score=1000):
        if source is None:
            source = 'assemblyline'
        # transcript feature
        f = GTFFeature()
        f.seqid = self.chrom
        f.source = source
        f.feature_type = 'transcript'
        f.start = self.start
        f.end = self.end
        f.score = score
        f.strand = strand_int_to_str(self.strand)
        f.phase = '.'
        f.attrs = self.attrs
        features = [f]
        # exon features
        for i,e in enumerate(self.exons):
            f = GTFFeature()
            f.seqid = self.chrom
            f.source = source
            f.feature_type = 'exon'
            f.start = e.start
            f.end = e.end
            f.score = score
            f.strand = strand_int_to_str(self.strand)
            f.phase = '.'
            f.attrs = self.attrs.copy()
            f.attrs["exon_number"] = i
            features.append(f)
        return features
