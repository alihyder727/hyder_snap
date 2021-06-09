#! /usr/bin/env python3
import argparse, glob
from chemistry_parser import *

chem_choices = glob.glob('*.md')
chem_choices = [x[:-4] for x in chem_choices]

parser = argparse.ArgumentParser()

parser.add_argument('--chem',
  default='kessler94',
  choices=chem_choices,
  help='select chemistry')

# Parse command-line inputs
args = vars(parser.parse_args())

vlist = read_variables('%s.md' % args['chem'])
coeffs = read_coefficients('%s.md' % args['chem'])
verb_str = read_verbatim('%s.md' % args['chem'])

header_str = write_header('%s.md' % args['chem'])
footer_str = write_footer()
defs_str = write_definitions(vlist)
coeffs_str = write_coefficients(coeffs)
react_str = write_all_reactions('%s.md' % args['chem'], vlist, ['simple', 'custom'])

result = '\n'.join([header_str,
  add_prefix(defs_str, '  '),
  add_prefix(verb_str, '  '),
  add_prefix(coeffs_str, '  '),
  add_prefix(react_str, '  '),
  footer_str])

print(result)
