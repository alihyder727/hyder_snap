#! /usr/bin/env python3
from snapy.harp.radio_model import create_inputs, run_forward, write_observation
from snapy.planet_gravity import jupiter_gravity
from snapy.harp.radio_arguments import *

if __name__ == '__main__':
# create model input file from template file
    if '.' not in args['input']:
        args['input'] += '.tmp'
    lat = float(args['lat'])
    args['grav'] = str(jupiter_gravity(lat))
    inpfile = create_inputs(args['input'], args)

# run hard foward
    pid = run_forward('harp_radio_js.ex', inpfile)

# write results
    outfile = write_observation(inpfile, pid + '-main.nc')
