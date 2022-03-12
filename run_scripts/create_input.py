#! /usr/bin/env python3
from snapy.harp.radio_model import create_input
from snapy.planet_gravity import jupiter_gravity
from snapy.harp.radio_arguments import *

# create model input file from template file
if '.' not in args['input']:
    args['input'] += '.tmp'
clat = float(args['clat'])
args['grav'] = str(jupiter_gravity(clat))
inpfile = create_input(args['input'], args)
