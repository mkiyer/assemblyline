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
    cur_t_id = 1
    expr_tot = 0.0
    for f in GTF.parse(open(sample.gtf_file)):
        if f.feature != 'exon':
            continue
        # rename transcript id
        t_id = f.attrs[GTF.Attr.TRANSCRIPT_ID]
        expr = f.attrs.get(gtf_expr_attr, '0')
        if t_id not in t_id_map:
            new_t_id = "%s.T%d" % (sample._id, cur_t_id)
            t_id_map[t_id] = new_t_id
            cur_t_id += 1
            expr_tot += float(expr)  # update total expression
        else:
            new_t_id = t_id_map[t_id]
        # update feature attributes
        attrs = ((GTF.Attr.TRANSCRIPT_ID, new_t_id),
                 (GTF.Attr.SAMPLE_ID, sample._id),
                 (GTF.Attr.REF, str(int(is_ref))),
                 (gtf_expr_attr, expr))
        f.attrs = collections.OrderedDict(attrs)
        # store feature
        if new_t_id not in t_dict:
            t_dict[new_t_id] = []
        t_dict[new_t_id].append(f)
    return t_dict, expr_tot


def add_sample_gtf(sample, gtf_expr_attr, output_fileh, stats_fileh,
                   is_ref=False):
    '''
    Reads and renames transfrags
    Normalizes expression by total filtered expression
    '''
    # read gtf file into dict of transcripts
    t_dict, expr_tot = _read_transfrags(sample, gtf_expr_attr, is_ref)

    # normalize expression by total
    exprs = []
    lengths = []
    for t_id, features in t_dict.iteritems():
        # sort features (exons) by start position
        features.sort(key=operator.attrgetter('start'))

        expr = float(features[0].attrs[gtf_expr_attr])
        expr_norm = 0.0 if is_ref else expr * (1.0e6 / expr_tot)
        exprs.append(expr_norm)
        lengths.append(sum((f.end - f.start) for f in features))
        # write transcript
        feature = _make_transcript_feature(features)
        feature.attrs[GTF.Attr.EXPRESSION] = expr_norm
        print >>output_fileh, str(feature)
        # write exons
        for i, feature in enumerate(features):
            feature.attrs[GTF.Attr.EXPRESSION] = expr_norm
            feature.attrs['exon_number'] = '%d' % (i + 1)
            print >>output_fileh, str(feature)

    # compute and write stats
    expr_quantiles = (scoreatpercentile(exprs, q) for q in range(0, 101))
    expr_quantiles = ','.join(map(str, expr_quantiles))
    length_quantiles = (int(round(scoreatpercentile(lengths, q)))
                        for q in range(0, 101))
    length_quantiles = ','.join(map(str, length_quantiles))
    fields = [sample._id, len(t_dict), expr_tot, expr_quantiles,
              length_quantiles]
    print >>stats_fileh, '\t'.join(map(str, fields))
