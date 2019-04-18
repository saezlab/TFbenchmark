#!/usr/bin/env python2
# coding: utf-8


import sys, getopt, argparse
from libs.fimo_parser import AnnotateTFBS



def main():
    parser = argparse.ArgumentParser(description="Script to parse TFBS predictions from fimo results and add ensemble regulatory and GERP conservation annotations")
    parser.add_argument('-outfile', help='Path to outfile containing the annotated TFBS', required = True)
    parser.add_argument('-infile', help='Path to fimo.txt file', required = True)
    args = parser.parse_args()

    AnnotateTFBS(args.infile, args.outfile)


if __name__ == '__main__':
    main()

