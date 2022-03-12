#! /usr/bin/env python3
from snapy.harp.radio_model import write_observation
import argparse, glob

parser = argparse.ArgumentParser()
parser.add_argument('-i',
        help='select input file')
parser.add_argument('-d',
        help='select nc output')
parser.add_argument('-o',
        default = 'none',
        help='output file name')
args = vars(parser.parse_args())

write_observation(args['i'], args['d'], args['o'])
