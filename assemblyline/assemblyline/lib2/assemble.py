'''
AssemblyLine: transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer

@author: mkiyer
'''
import collections
import bisect
import logging

from gtf import GTF, GTFError

Exon = collections.namedtuple('Exon', ['start', 'end'])


class AssemblyError(Exception):
    pass


def find_exon_boundaries(transcripts):
    '''
    input: a list of transcripts (not Node objects, these are transcripts)
    parsed directly from an input file and not added to an isoform graph
    yet.

    output: sorted list of exon boundaries
    '''
    exon_boundaries = set()
    # first add introns to the graph and keep track of
    # all intron boundaries
    for transcript in transcripts:
        # add transcript exon boundaries
        for exon in transcript.exons:
            # keep track of positions where introns can be joined to exons
            exon_boundaries.add(exon.start)
            exon_boundaries.add(exon.end)
    # sort the intron boundary positions and add them to interval trees
    return sorted(exon_boundaries)


def split_exon(exon, boundaries):
    """
    partition the exon given list of node boundaries

    generator yields (start,end) intervals for exon
    """
    if exon.start == exon.end:
        return
    # find the indexes into the intron boundaries list that
    # border the exon.  all the indexes in between these two
    # are overlapping the exon and we must use them to break
    # the exon into pieces
    start_ind = bisect.bisect_right(boundaries, exon.start)
    end_ind = bisect.bisect_left(boundaries, exon.end)
    if start_ind == end_ind:
        yield exon.start, exon.end
    else:
        yield exon.start, boundaries[start_ind]
        for j in xrange(start_ind, end_ind-1):
            yield boundaries[j], boundaries[j+1]
        yield boundaries[end_ind-1], exon.end


def split_exons(t, boundaries):
    # split exons that cross boundaries and to get the
    # nodes in the transcript path
    for exon in t.exons:
        for start, end in split_exon(exon, boundaries):
            yield start, end


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


class Locus(object):

    class NodeData:
        __slots__ = ('strands', 'samples', 'exprs')

        def __init__(self):
            self.strands = {GTF.POS_STRAND: False,
                            GTF.NEG_STRAND: False,
                            GTF.NO_STRAND: False}
            self.samples = {GTF.POS_STRAND: set(),
                            GTF.NEG_STRAND: set(),
                            GTF.NO_STRAND: set()}
            self.exprs = {GTF.POS_STRAND: 0.0,
                          GTF.NEG_STRAND: 0.0,
                          GTF.NO_STRAND: 0.0}

    def __init__(self):
        self.chrom = None
        self.boundaries = None
        self.node_data = None
        self.strand_transfrags = None

    @staticmethod
    def create(transfrags):
        self = Locus()
        self.boundaries = find_exon_boundaries(transfrags)
        self.node_data = collections.defaultdict(Locus.NodeData)
        self.strand_transfrags = {GTF.POS_STRAND: [],
                                  GTF.NEG_STRAND: [],
                                  GTF.NO_STRAND: []}
        self.chrom = None
        for t in transfrags:
            if self.chrom is None:
                self.chrom = t.chrom
            elif self.chrom != t.chrom:
                logging.error('Locus.create: "chrom" mismatch')
                raise AssemblyError('Locus.create: transfrag chromosomes do '
                                    'not match')
            self._add_transfrag(t)
            self.strand_transfrags[t.strand].append(t)
        return self

    def _add_transfrag(self, t):
        for n in split_exons(t, self.boundaries):
            nd = self.node_data[n]
            nd.strands[t.strand] = True
            if not t.is_ref:
                nd.exprs[t.strand] += t.expr
                nd.samples[t.strand].add(t.sample_id)

    def _remove_transfrag(self, t):
        for n in split_exons(t, self.boundaries):
            nd = self.node_data[n]
            nd.strands[t.strand] = False
            nd.samples[t.strand].remove(t.sample_id)
            nd.exprs[t.strand] -= t.expr

    def _predict_strand(self, nodes):
        total_length = sum((n[1]-n[0]) for n in nodes)
        strand_expr = {GTF.POS_STRAND: 0.0, GTF.NEG_STRAND: 0.0}
        strand_length = {GTF.POS_STRAND: 0, GTF.NEG_STRAND: 0}

        for n in nodes:
            nd = self.node_data[n]
            length = n[1] - n[0]
            frac_length = float(length) / total_length
            strand_expr['+'] += nd.exprs['+'] * frac_length
            strand_expr['-'] += nd.exprs['-'] * frac_length
            if nd.strands['+']:
                strand_length['+'] += frac_length
            if nd.strands['-']:
                strand_length['-'] += frac_length

        # if transfrag supported by stranded coverage choose strand with
        # greatest length-normalized expression level
        total_expr = sum(strand_expr.values())
        if total_expr > 0:
            if strand_expr['+'] > strand_expr['-']:
                return '+'
            elif strand_expr['-'] > strand_expr['+']:
                return '-'
        # if transfrag not supported by stranded coverage choose strand
        # with greatest support from reference transcripts
        total_strand_length = sum(strand_length.values())
        if total_strand_length > 0:
            if strand_length['+'] > strand_length['-']:
                return '+'
            elif strand_length['-'] > strand_length['+']:
                return '-'
        return '.'

    def predict_unknown_strands(self):
        # iteratively predict strand until no new transfrags can be
        # predicted
        iterations = 0
        num_resolved = 0
        while(len(self.strand_transfrags['.']) > 0):
            resolved = []
            unresolved = []
            for t in self.strand_transfrags['.']:
                nodes = list(split_exons(t, self.boundaries))
                new_strand = self._predict_strand(nodes)
                if new_strand != '.':
                    resolved.append((t, new_strand))
                    num_resolved += 1
                else:
                    unresolved.append(t)
            # break when no new transfrags could have strand predicted
            if len(resolved) == 0:
                break
            # clear Locus data for transfrags with unknown strand and re-add
            # them to the Locus
            for t, new_strand in resolved:
                self._remove_transfrag(t)
                t.strand = new_strand
                self._add_transfrag(t)
                self.strand_transfrags[t.strand].append(t)
            self.strand_transfrags['.'] = unresolved
            iterations += 1

        if iterations > 0:
            logging.debug('predict_unknown_strands: %d iterations' %
                          iterations)
        return num_resolved

    @staticmethod
    def open_bedgraph(file_prefix):
        attrs = ('expression', 'recurrence')
        strand_names = {'+': 'pos', '-': 'neg', '.': 'none'}
        filehs = {}
        for a in attrs:
            filehs[a] = {}
            for s in ('+', '-', '.'):
                filename = '%s.%s.%s.bedgraph' % (file_prefix, a,
                                                  strand_names[s])
                filehs[a][s] = open(filename, 'w')
        return filehs

    @staticmethod
    def close_bedgraph(filehs):
        for adict in filehs.itervalues():
            for fileh in adict.itervalues():
                fileh.close()

    def get_bedgraph_data(self):
        '''
        Returns node attribute data in tuples
            chrom, start, end, strand, total expression, sample recurrence
        '''
        for n in sorted(self.node_data):
            nd = self.node_data[n]
            for strand in ('+', '-', '.'):
                yield (self.chrom, n[0], n[1], strand,
                       nd.exprs[strand], len(nd.samples[strand]))

    def write_bedgraph(self, bgfiledict):
        '''
        bgfiledict: dictionary structure containing file handles opened
                    for writing obtained using Locus.open_bedgraph()
        '''
        for tup in self.get_bedgraph_data():
            chrom, start, end, strand, expr, recur = tup
            if expr > 0:
                line = '\t'.join(map(str, [chrom, start, end, expr]))
                print >>bgfiledict['expression'][strand], line
            if recur > 0:
                line = '\t'.join(map(str, [chrom, start, end, recur]))
                print >>bgfiledict['recurrence'][strand], line
