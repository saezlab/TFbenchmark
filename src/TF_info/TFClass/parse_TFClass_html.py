#!/usr/bin/env python2
# coding: utf-8


import sys, getopt, argparse
from libs.parseTFClass import ParseTFClassHTMLurl



def main():
    parser = argparse.ArgumentParser(description="Script to parse TFClass html")
    parser.add_argument('-outfile', help='Path to outfile', required=True)
    args = parser.parse_args()
    
    ParseTFClassHTMLurl(args.outfile, URL='http://www.edgar-wingender.de/huTF_classification.html')

if __name__ == '__main__':
    main()

