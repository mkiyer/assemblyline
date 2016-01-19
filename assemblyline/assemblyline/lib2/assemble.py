'''
AssemblyLine: transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer

@author: mkiyer
'''
import collections
import bisect
import logging
import networkx as nx

from gtf import GTF


class AssemblyError(Exception):
    pass


def find_splice_sites(transfrags):
    '''
    input: a list of transfrags

    output: sorted list of exon boundaries
    '''
    splice_sites = set()
    # first add introns to the graph and keep track of
    # all intron boundaries
    for t in transfrags:
        # add intron boundaries
        for start, end in t.iterintrons():
            splice_sites.add(start)
            splice_sites.add(end)
    # sort the intron boundary positions
    return sorted(splice_sites)


def find_transfrag_nodes(t, splice_sites):
    '''
    assumes transfrag exons are sorted by genomic position
    '''
    if t.start == t.end:
        return

    # subset splice sites list to boundaries of transfrag
    t_start_ind = bisect.bisect_right(splice_sites, t.start)
    t_end_ind = bisect.bisect_left(splice_sites, t.end)
    if t_start_ind == t_end_ind:
        yield t.start, t.end
    splice_sites = splice_sites[t_start_ind:t_end_ind]

    for exon in t.exons:
        # find the indexes into the splice sites list that border the exon.
        start_ind = bisect.bisect_right(splice_sites, exon.start)
        end_ind = bisect.bisect_left(splice_sites, exon.end)
        if start_ind == end_ind:
            yield exon.start, exon.end
        else:
            yield exon.start, splice_sites[start_ind]
            # all the splice sites in between the exon borders must overlap
            for j in xrange(start_ind, end_ind-1):
                yield splice_sites[j], splice_sites[j+1]
            yield splice_sites[end_ind-1], exon.end
        # subset splice sites as we move along the transcript
        splice_sites = splice_sites[end_ind-1:]


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

    def _check_strand_ambiguous(self, nodes):
        '''
        Checks list of nodes for strandedness. If strand is unambiguous,
        then return pos or neg strand. If ambiguous, return unstranded.
        '''
        strands = {GTF.POS_STRAND: False,
                   GTF.NEG_STRAND: False}
        for n in nodes:
            nd = self.node_data[n]
            if nd.strands['+'] or nd.exprs['+'] > 0:
                strands['+'] = True
            if nd.strands['-'] or nd.exprs['-'] > 0:
                strands['-'] = True
        if strands['+'] and strands['-']:
            return '.'
        elif strands['+']:
            return '+'
        elif strands['-']:
            return '-'
        return '.'

    def impute_unknown_strands(self):
        # iteratively predict strand of unstranded transfrags
        # stop when no new transfrag strands can be imputed
        iterations = 0
        num_resolved = 0
        while(len(self.strand_transfrags['.']) > 0):
            resolved = []
            unresolved = []
            for t in self.strand_transfrags['.']:
                nodes = list(split_exons(t, self.boundaries))
                new_strand = self._check_strand_ambiguous(nodes)
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

    def create_directed_graph(self, strand):
        '''
        build strand-specific graph
        '''
        def add_node(G, n, expr):
            """add node to graph"""
            if n not in G:
                G.add_node(n, length=(n.end - n.start), expr=0.0)
            nd = G.node[n]
            nd['expr'] += expr

        # initialize transcript graph
        transfrags = self.strand_transfrags[strand]
        boundaries = find_exon_boundaries(transfrags)
        G = nx.DiGraph()

        # add transcripts
        for t in transfrags:
            # split exons that cross boundaries and get the
            # nodes that made up the transfrag
            nodes = [n for n in split_exons(t, boundaries)]
            if strand == '-':
                nodes.reverse()
            # add nodes/edges to graph
            u = nodes[0]
            add_node(G, u, t.expr)
            for i in xrange(1, len(nodes)):
                v = nodes[i]
                add_node(G, v, t.expr)
                G.add_edge(u, v)
                u = v

        # set graph attributes
        G.graph['boundaries'] = boundaries
        G.graph['strand'] = strand
        return G

    def create_directed_graphs(self):
        for strand, transfrags in self.strand_transfrags.iteritems():
            # create strand-specific directed graph
            G = self.create_directed_graph(strand)
            # collapse consecutive nodes in graph
            H, node_chain_map = collapse_strand_specific_graph(G, introns=True)
            # get connected components of graph which represent independent genes
            # unconnected components are considered different genes
            Gsubs = nx.weakly_connected_component_subgraphs(H)

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
