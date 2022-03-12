#! /usr/bin/env python3
from snapy.harp.radio_model import create_input, run_forward, write_observation
from snapy.planet_gravity import jupiter_gravity
from snapy.harp.radio_arguments import *
import os

if __name__ == '__main__':
# create model input file from template file
    if '.' not in args['input']:
        args['input'] += '.tmp'
    clat = float(args['clat'])
    args['grav'] = str(jupiter_gravity(clat))
    inpfile = create_input(args['input'], args)

# run hard foward
    pid = run_forward(args['exe'], inpfile)

# write results
    if os.path.exists(pid + '-mcmc.nc'):
      outfile = write_observation(inpfile, pid + '-mcmc.nc')
    elif os.path.exists(pid + '-main.nc'):
      outfile = write_observation(inpfile, pid + '-main.nc')
    else:
      raise ValueError("neither main nor mcmc file exists. Abort.")
