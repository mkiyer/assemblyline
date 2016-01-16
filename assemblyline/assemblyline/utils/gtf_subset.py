import argparse
import logging

from assemblyline.lib2.gtf import GTF


def main():
    logging.basicConfig(level=logging.DEBUG)
    parser = argparse.ArgumentParser()
    parser.add_argument("gtf_file")
    parser.add_argument("gtf_attr")
    parser.add_argument("values_file")
    args = parser.parse_args()

    values = set()
    with open(args.values_file) as fileh:
        for line in fileh:
            line = line.strip().split('\t')
            values.add(line[0])

    with open(args.gtf_file) as fileh:
        for f in GTF.parse(fileh):
            val = f.attrs[args.gtf_attr]
            if val in values:
                print str(f)

if __name__ == '__main__':
    main()
