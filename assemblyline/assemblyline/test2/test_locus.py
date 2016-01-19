import os

from assemblyline.lib2.gtf import GTF
from assemblyline.lib2.transfrag import Transfrag
from assemblyline.lib2.assemble import find_splice_sites, find_transfrag_nodes

INPUT_FILE_DIR = "input_files"


def get_gtf_path(filename):
    return os.path.join(os.path.dirname(__file__), INPUT_FILE_DIR, filename)


def read_gtf(filename):
    return list(GTF.parse_loci(open(get_gtf_path(filename))))


def test_parse_loci():
    loci = read_gtf('parse_loci.gtf')
    assert len(loci) == 3
    assert loci[0][0] == ('chrTest1', 10, 50)
    assert loci[1][0] == ('chrTest1', 50, 200)
    assert loci[2][0] == ('chrTest2', 100, 200)


def test_find_splice_sites():
    loci = read_gtf('splice_sites.gtf')
    assert len(loci) == 1
    interval, gtf_lines = loci[0]
    assert interval == ('chr1', 0, 500)
    transfrag_dict = Transfrag.parse_gtf(gtf_lines)
    splice_sites = find_splice_sites(transfrag_dict.values())
    assert splice_sites == [100, 200, 250, 300, 400]


def test_find_transfrag_nodes():
    loci = read_gtf('splice_sites.gtf')
    interval, gtf_lines = loci[0]
    t_dict = Transfrag.parse_gtf(gtf_lines)
    splice_sites = find_splice_sites(t_dict.values())
    # A
    t = t_dict['A']
    nodes = list(find_transfrag_nodes(t, splice_sites))
    assert nodes == 1
