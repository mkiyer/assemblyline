'''
Created on Dec 20, 2013

@author: mkiyer
'''
'''
Created on Dec 18, 2013

@author: mkiyer
'''
import os
import sys
import logging
import argparse
import collections

from assemblyline.lib.gtf import GTFFeature
from assemblyline.lib.base import Category
from assemblyline.lib.transcript import parse_gtf

# Translation of our transcript categories to GENCODE gene types was done 
# after careful examination of the GENCODE gene type description web page:
# http://www.gencodegenes.org/gencode_biotypes.html
PROTEIN_CATEGORY_MAP = {Category.OPP_STRAND: 'antisense',
                        Category.ENCOMPASSING_SAME_STRAND: 'sense_overlapping',
                        Category.ENCOMPASSING_OPP_STRAND: 'antisense',
                        Category.INTRONIC_SAME_STRAND: 'sense_intronic',
                        Category.INTRONIC_OPP_STRAND: 'antisense',
                        Category.INTRONIC_SAME_STRAND: 'sense_intronic',
                        Category.INTERLEAVING_SAME_STRAND: 'sense_intronic',
                        Category.INTERLEAVING_OPP_STRAND: 'antisense',
                        Category.INTERGENIC: 'lincRNA'}

# Map GENCODE gene type to categories
GENCODE_CATEGORY_MAP = {'IG_C_gene': 'protein_coding',
                     'IG_D_gene': 'protein_coding',
                     'IG_J_gene': 'protein_coding',
                     'IG_V_gene': 'protein_coding',
                     'TR_C_gene': 'protein_coding',
                     'TR_J_gene': 'protein_coding',
                     'TR_V_gene': 'protein_coding',
                     'TR_D_gene': 'protein_coding',
                     'TEC': 'protein_coding',
                     'nonsense_mediated_decay': 'protein_coding',
                     'non_stop_decay': 'protein_coding',                     
                     'retained_intron': 'protein_coding', 
                     'protein_coding': 'protein_coding', 
                     'ambiguous_orf': 'protein_coding', 
                     'Mt_rRNA': 'ncRNA',
                     'Mt_tRNA': 'ncRNA',
                     'miRNA': 'ncRNA',
                     'misc_RNA': 'ncRNA',
                     'rRNA': 'ncRNA',
                     'snRNA': 'ncRNA',
                     'snoRNA': 'ncRNA',
                     '3prime_overlapping_ncrna': 'ncRNA',
                     'processed_transcript': 'lncRNA', 
                     'lincRNA': 'lncRNA',
                     'sense_intronic': 'lncRNA', 
                     'sense_overlapping': 'lncRNA',   
                     'antisense': 'lncRNA', 
                     'IG_C_pseudogene': 'pseudogene',
                     'IG_J_pseudogene': 'pseudogene',
                     'IG_V_pseudogene': 'pseudogene',
                     'TR_V_pseudogene': 'pseudogene',
                     'TR_J_pseudogene': 'pseudogene',
                     'Mt_tRNA_pseudogene': 'pseudogene',
                     'tRNA_pseudogene': 'pseudogene',
                     'snoRNA_pseudogene': 'pseudogene',
                     'snRNA_pseudogene': 'pseudogene',
                     'scRNA_pseudogene': 'pseudogene',
                     'rRNA_pseudogene': 'pseudogene',
                     'misc_RNA_pseudogene': 'pseudogene',
                     'miRNA_pseudogene': 'pseudogene',
                     'pseudogene': 'pseudogene',
                     'processed_pseudogene': 'pseudogene', 
                     'polymorphic_pseudogene': 'pseudogene', 
                     'retrotransposed': 'pseudogene', 
                     'transcribed_processed_pseudogene': 'pseudogene',
                     'transcribed_unprocessed_pseudogene': 'pseudogene',
                     'unitary_pseudogene': 'pseudogene', 
                     'unprocessed_pseudogene': 'pseudogene'} 

def main():
    logging.basicConfig(level=logging.DEBUG,
                        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s")
    parser = argparse.ArgumentParser()
    parser.add_argument('--rename', dest='rename', action='store_true')
    parser.add_argument('gtf_file')
    args = parser.parse_args()
    gtf_file = args.gtf_file
    rename = args.rename
    if not os.path.exists(gtf_file):
        parser.error("GTF file '%s' not found" % (gtf_file))
    # parse transcripts
    num_transcripts = 0
    # keep track of redundant gene/transcript counts
    gene_count = collections.defaultdict(lambda: 1)
    transcript_count = collections.defaultdict(lambda: 1)
    t_id_map = {}
    g_id_map = {}
        
    for transcripts in parse_gtf(open(gtf_file)):
        for t in transcripts:
            catstr = t.attrs['category']
            catint = Category.to_int(catstr)
            gene_type = t.attrs.get('gene_type', None)
            ref_gene_type = t.attrs['ref_gene_type']
            if catint == Category.SAME_STRAND:
                # impute gene type
                new_gene_type = ref_gene_type
            else:
                if gene_type == 'protein_coding':
                    # don't change protein coding genes
                    new_gene_type = gene_type
                elif t.length < 250:
                    # categorize small RNA separately
                    new_gene_type = 'misc_RNA'
                else:
                    if ref_gene_type == 'protein_coding':
                        # categorize based on overlap with reference
                        new_gene_type = PROTEIN_CATEGORY_MAP[catint]
                    else:
                        # reference is also non-coding
                        new_gene_type = 'lincRNA'
            # get gene category
            gene_category = GENCODE_CATEGORY_MAP[new_gene_type]
            new_gene_name = None
            # build new gene name
            if rename:
                ref_gene_name = t.attrs['ref_gene_name']
                if ref_gene_name == 'None':
                    new_gene_name = str(t.attrs['source'])
                elif catint == Category.SAME_STRAND:
                    new_gene_name = str(ref_gene_name)
                else:
                    new_gene_name = '%s.%s' % (ref_gene_name, catstr)
                gene_id = t.attrs['gene_id']
                t_id = t.attrs['transcript_id']
                
    t_id_map = {}
    g_id_map = {}
    tss_id_map = {}
    for feature in GTFFeature.parse(open(args.gtf_file)):
        t_id = feature.attrs['transcript_id']
        g_id = feature.attrs['gene_id']
        tss_id = feature.attrs['tss_id']
        if t_id not in t_id_map:
            new_t_id = "T%06d" % (cur_t_id)
            t_id_map[t_id] = new_t_id
            cur_t_id += 1
        else:
            new_t_id = t_id_map[t_id]
        if g_id not in g_id_map:
            new_g_id = "G%06d" % (cur_g_id)
            g_id_map[g_id] = new_g_id
            cur_g_id += 1
        else:                
                
                gene_num = gene_count[new_gene_name]
                gene_count[new_gene_name] += 1
                new_gene_name = '%s.%d' % (new_gene_name, gene_num)
                transcript_num = transcript_count[new_gene_name]
                transcript_count[new_gene_name] += 1
                new_gene_name = '%s.%d' % (new_gene_name, transcript_num)
            # write new attributes
            for f in t.to_gtf_features(source='assemblyline', score=1000):
                f.attrs['gene_type'] = new_gene_type
                f.attrs['gene_category'] = gene_category 
                if rename:
                    if 'gene_name' in f.attrs:
                        f.attrs['orig_gene_name'] = f.attrs['gene_name']
                    f.attrs['gene_name'] = new_gene_name        
                print str(f)
            num_transcripts += 1
    return 0

if __name__ == '__main__':
    sys.exit(main())