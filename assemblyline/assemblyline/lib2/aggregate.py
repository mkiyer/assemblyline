'''
AssemblyLine: transcriptome meta-assembly from RNA-Seq
Copyright (C) 2012-2015 Matthew Iyer

@author: mkiyer
'''
import collections
import operator

from stat import scoreatpercentile
from gtf import GTF


def _make_transcript_feature(exon_features):
    f = GTF.Feature()
    f.seqid = exon_features[0].seqid
    f.source = exon_features[0].source
    f.feature = 'transcript'
    f.start = exon_features[0].start
    f.end = exon_features[-1].end
    f.score = exon_features[0].score
    f.strand = exon_features[0].strand
    f.phase = '.'
    f.attrs = exon_features[0].attrs.copy()
    if 'exon_number' in f.attrs:
        del f.attrs['exon_number']
    return f


def _read_transfrags(sample, gtf_expr_attr, is_ref=False):
    '''
    Process individual sample GTF file
      - Reads entire GTF file into memory.
      - Renames "gene_id" and "transcript_id" attributes for
        consistency and to conserve space.
    '''
    t_dict = collections.OrderedDict()
    t_id_map = {}
    t_expr_map = {}
    cur_t_id = 1
    for f in GTF.parse(open(sample.gtf_file)):
        t_id = f.attrs[GTF.Attr.TRANSCRIPT_ID]
        if f.feature == 'transcript':
            # save expression
            expr = f.attrs[gtf_expr_attr]
            t_expr_map[t_id] = expr
            # rename transcript id
            if t_id not in t_id_map:
                new_t_id = "%s.T%d" % (sample._id, cur_t_id)
                t_id_map[t_id] = new_t_id
                cur_t_id += 1
                t_dict[new_t_id] = []    # init t_dict
        elif f.feature == 'exon':
            # lookup expression
            if is_ref:
                expr = 0.0
            else:
                expr = float(t_expr_map[t_id])
            new_t_id = t_id_map[t_id]
            # store exon feature
            attrs = ((GTF.Attr.TRANSCRIPT_ID, new_t_id),
                     (GTF.Attr.SAMPLE_ID, sample._id),
                     (GTF.Attr.REF, str(int(is_ref))),
                     (gtf_expr_attr, expr))
            f.attrs = collections.OrderedDict(attrs)
            t_dict[new_t_id].append(f)
    return t_dict


def add_sample_gtf(sample, gtf_expr_attr, output_fileh, stats_fileh,
                   is_ref=False):
    '''
    Reads and renames transfrags
    Normalizes expression by total filtered expression
    '''
    # read gtf file into dict of transcripts
    t_dict = _read_transfrags(sample, gtf_expr_attr, is_ref)

    # normalize expression by total
    exprs = []
    lengths = []
    for t_id, features in t_dict.iteritems():
        # sort features (exons) by start position
        features.sort(key=operator.attrgetter('start'))
        expr = 0.0 if is_ref else float(features[0].attrs[gtf_expr_attr])
        exprs.append(expr)
        lengths.append(sum((f.end - f.start) for f in features))
        # write transcript
        feature = _make_transcript_feature(features)
        feature.attrs[GTF.Attr.EXPRESSION] = expr
        print >>output_fileh, str(feature)
        # write exons
        for i, feature in enumerate(features):
            feature.attrs[GTF.Attr.EXPRESSION] = expr
            feature.attrs['exon_number'] = '%d' % (i + 1)
            print >>output_fileh, str(feature)

    # compute and write stats
    expr_quantiles = (scoreatpercentile(exprs, q) for q in range(0, 101))
    expr_quantiles = ','.join(map(str, expr_quantiles))
    length_quantiles = (int(round(scoreatpercentile(lengths, q)))
                        for q in range(0, 101))
    length_quantiles = ','.join(map(str, length_quantiles))
    fields = [sample._id, len(t_dict), expr_quantiles,
              length_quantiles]
    print >>stats_fileh, '\t'.join(map(str, fields))
