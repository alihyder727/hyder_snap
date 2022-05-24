from netCDF4 import Dataset
from astropy.io import fits
#from multiprocessing import Process
from tqdm import tqdm
from subprocess import check_call
from numpy import zeros
from multiprocessing import Pool
import shutil, os

def main_to_mcmc(fname, field = 'main'):
  #fname = 'vla_ideal_saturn-n1000'
  #os.remove('%s-mcmc.nc' % fname)
  #shutil.copy('%s-main.nc' % fname, '%s-tmp.nc' % fname)

  if field == 'main':
    data = Dataset('%s-main.nc' % fname, 'r+')
  else:
    data = Dataset('%s.%s.nc' % (fname, field), 'r+')
  msk = fits.open('%s.fits' % fname)[3].data
  nstep, nwalker = msk.shape

  data['time'][:] = range(nstep)
  data['time'].long_name = "steps"
  data['time'].units = "1"
  data['x2'][:] = range(len(data['x2'][:]))
  data['x2'].long_name = "model"
  data['x2'].units = "1"
  data['x3'][:] = range(len(data['x3'][:]))
  data['x3'].long_name = "walker"
  data['x3'].units = "1"

  for key in data.variables.keys():
    print('processing %s ...' % key)
    dimensions = data[key].dimensions

    # copy new state variable
    # i1 records the most recent new state
    i1 = 0
    if dimensions == ('time', 'x1', 'x2', 'x3'):
      new_data = data[key][:,:,:,:]
      for j in tqdm(range(nwalker)):
        for i2 in range(1,nstep):
          if msk[i2,j] == 1:
            new_data[i1:i2,:,:,j] = new_data[i1,:,:,j]
            i1 = i2
        new_data[i1:,:,:,j] = new_data[i1,:,:,j]
      # save data
      data[key][:] = new_data
    elif dimensions == ('time', 'x2', 'x3'):
      new_data = data[key][:,:,:]
      for j in tqdm(range(nwalker)):
        for i2 in range(1,nstep):
          if msk[i2,j] == 1:
            new_data[i1:i2,:,:,j] = new_data[i1,:,:,j]
            i1 = i2
        new_data[i1:,:,j] = new_data[i1,:,j]
      # save data
      data[key][:] = new_data
    elif dimensions == ('time', 'x3'):
      new_data = data[key][:,:]
      for j in tqdm(range(nwalker)):
        for i2 in range(1,nstep):
          if msk[i2,j] == 1:
            new_data[i1:i2,:,:,j] = new_data[i1,:,:,j]
            i1 = i2
        new_data[i1:,j] = new_data[i1,j]
      # save data
      data[key][:] = new_data

  data.close()

  if field == 'main':
    check_call('ncks -d time,0,%d %s-main.nc %s-mcmc.nc' %
      (nstep-1, fname, fname), shell = True)
    os.remove('%s-main.nc' % fname)
  else:
    check_call('ncks -d time,0,%d %s.%s.nc %s.%s-mcmc.nc' %
      (nstep-1, fname, field, fname, field), shell = True)
    os.remove('%s.%s.nc' % (fname, field))
