from ..athena.athena_read import athinput
from numpy import logspace, log10, zeros, arccos, array, pi
from netCDF4 import Dataset
from .utils import get_rt_bands, get_ray_out
import re, subprocess

def create_input(tmpfile, args):
    with open(tmpfile, 'r') as file:
      tmpinp = file.read()

    if ':' in args['plevel']:
      pmax, pmin, np = tuple(map(float, args['plevel'].split(':')))
      np = int(np)
      plevel = logspace(log10(pmax), log10(pmin), np)
      plevel = ['%10.2f'%x for x in plevel]
    elif ',' in args['plevel']:
      plevel = args['plevel'].split(',')
      np = len(plevel)

    if args['tem'] == '0':
      Tp = ['0.']*np
    else:
      Tp = args['tem'].split(',')
    assert(len(Tp) == np)

    if args['nh3'] == '0':
      NH3p  = ['0.']*np
    else:
      NH3p = args['nh3'].split(',')
    assert(len(NH3p) == np)

  # adjust minimum number of walkers
    #args['nwalker'] = str(max(int(args['nwalker']), 2*np))

    var = [x for x in args['var'].split()]

    name = tmpfile.split('/')[-1].split('.')[0]
    if args['output'] != '':
      name += '-' + args['output']

    inpfile = re.sub('\[problem_id\]', name, tmpinp)
    inpfile = re.sub('\[logname\]', name, inpfile)
    inpfile = re.sub('\[obsname\]', args['obs'], inpfile)
    inpfile = re.sub('\[grav\]', '-' + args['grav'], inpfile)
    inpfile = re.sub('\[plevel\]', ' '.join(plevel), inpfile)
    inpfile = re.sub('\[variables\]', ' '.join(var), inpfile)
    inpfile = re.sub('\[Tp\]', ' '.join(Tp), inpfile)
    inpfile = re.sub('\[NH3p\]', ' '.join(NH3p), inpfile)
    inpfile = re.sub('\[Tstd\]', args['sT'], inpfile)
    inpfile = re.sub('\[Xstd\]', args['sNH3'], inpfile)
    inpfile = re.sub('\[Tlen\]', args['zT'], inpfile)
    inpfile = re.sub('\[Xlen\]', args['zNH3'], inpfile)
    inpfile = re.sub('\[lwalker\]', str(int(args['nwalker'])//int(args['nodes'])), inpfile)
    if args['M']:
        inpfile = re.sub('\[M\]', 'true', inpfile)
    else:
        inpfile = re.sub('\[M\]', 'false', inpfile)
    if args['d']:
        inpfile = re.sub('\[diff\]', 'true', inpfile)
    else:
        inpfile = re.sub('\[diff\]', 'false', inpfile)

    for key in args.keys():
      try :
        inpfile = re.sub('\[%s\]' % key, args[key], inpfile)
      except TypeError :
        pass

    with open(name + '.inp', 'w') as file:
        file.write(inpfile)
    print('Input file written to %s.inp' % name)
    return name + '.inp'

def run_forward(exefile, inpfile):
    script = [exefile, '-i', inpfile]
    process = subprocess.Popen(script,
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE)
    while True:
        output = process.stdout.readline()
        #err = process.stderr.readline()
        #if (err != ''):
        #  raise Exception(err.decode('UTF-8'))
        if output == b'' and process.poll() is not None:
            break
        if output:
            print(output.decode('UTF-8'), end = '')
    process.poll()
    #print(out.decode('UTF-8'), end = '\r')
    #print(err.decode('UTF-8'))

    out, err = subprocess.Popen('./combine.py',
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE).communicate()
    print(out.decode('UTF-8'), end = '')
    print(err.decode('UTF-8'), end = '')

    inp = athinput(inpfile)
    return inp['job']['problem_id']

def write_observation(inpfile, datafile, output = 'none'):
    freq = get_rt_bands(inpfile)[:,0]
    num_bands = len(freq)

# read radiation toa
    data = Dataset(datafile, 'r')
    amu = arccos(data['mu_out'][:])/pi*180.
    amu = amu.reshape((num_bands, -1))
    num_dirs = amu.shape[1]
    tb = data['radiance'][0,:,:,0]
    tb = tb.reshape((num_bands, num_dirs, -1))

# write to file
    if output == 'none':
      outfile = '.'.join(inpfile.split('.')[:-1]) + '.out'
    else:
      outfile = output
    with open(outfile, 'w') as file:
      for k in range(tb.shape[2]):
        file.write('# Brightness temperatures of input model %s - model %d\n' % (inpfile, k))
        file.write('%12s' % '# Freq (GHz)')
        for i in range(num_dirs):
          file.write('%10.2f' % amu[0,i])
        file.write('\n')
        for i in range(num_bands):
          file.write('%12.2f' % freq[i])
          for j in range(num_dirs):
            file.write('%10.2f' % tb[i,j,k])
          file.write('\n')
    print('Brightness temperatures written to %s' % outfile)
    return outfile
