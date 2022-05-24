#! /usr/bin/env python3
import yaml, os, re, argparse, glob, shutil
from itertools import product

parser = argparse.ArgumentParser()
parser.add_argument('-i', '--input',
    choices = glob.glob('*.yml'),
    required = True,
    help = 'yaml instruction file'
    )
parser.add_argument('-hpc',
    action='store_true',
    default=False,
    help = 'enable using hpc'
    )
args = vars(parser.parse_args())

def replace_environ_var(envs, old_str):
  new_str, count = old_str, 0
  while '$' in new_str:
    # replace environment variables
    m = re.search('\${{\s*([^\s]*)\s*}}', new_str)
    if m:
      src = m.group(0)
      dst = m.group(1)
      new_str = re.sub('\${{\s*([^\s]*)\s*}}', os.environ[dst], new_str)
    # replace user defined variables
    for key,var in envs.items(): 
      new_str = re.sub('\$' + key, str(var), new_str)
    count += 1
    if count > 10:
      raise RuntimeError('environment variable(s) not found in "%s"' %  new_str)
  return new_str

def replace_vector_var(vecs, old_str):
  new_str, count = old_str, 0
  for key,var in vecs.items():
    new_str = re.sub('\$' + key, str(var), new_str)

def write_job_steps(file, job, detail, major, minor, hpc = False):
  print('Processing %s.%s-%s' % (major, job, minor))
  # write header
  file.write('# %d.%d. execute job %s' % (major, minor, job))
  for key, value in detail['strategy']['vector'].items():
    file.write(', %s = %s' % (key, value[minor]))
  file.write('\n')

  # write submission file
  use_hpc = False
  if ('NODES' in detail['env']) and hpc: use_hpc = True

  # count matrix jobs
  sub_jobs_matrix = []
  for key, value in detail['strategy']['matrix'].items():
    sub_jobs_matrix.append([(key,x) for x in value])

  sub_jobs_list = list(map(dict, product(*sub_jobs_matrix)))
  for it, mdict in enumerate(sub_jobs_list):
    if use_hpc:
      fsub_name = '%s-%s-%s' % (job, minor, it)
      with open('submit-gl.tmp', 'r') as fsub:
        tmpsub = fsub.read()
      subfile = re.sub('<name>', fsub_name, tmpsub)
      subfile = re.sub('<nodes>', str(detail['env']['NODES']), subfile)
      fsub = open(fsub_name + '.sub', 'w')
      fsub.write(subfile)
      original_file = file
      file = fsub

    file.write('# %d.%d.%d. ' % (major, minor, it))

    # set environment variable from matrix
    for key, value in mdict.items():
      detail['env'][key] = mdict[key]
      file.write('%s = %s ' % (key, value))
    file.write('\n')

    # write jobs stesp
    for step in detail['steps']:
      # write step name
      step_name = replace_environ_var(detail['env'], step['name'])
      file.write('## %s\n' % step_name)

      # write step command
      step_run = replace_environ_var(detail['env'], step['run'])
      if use_hpc:
        step_run = step_run.split('&')[0]
        # replace MPI commands
        step_run = re.sub('mpiexec', 'srun', step_run)
        # replace mv commands
        m = re.match('\s*mv\s+([^\s]+)\s+([^\s]+)', step_run)
        if m:
          dst = m.group(2)
          #step_run = re.sub(dst, '../../' + dst, step_run)
      if step['run'][-1] == '\n':
        file.write('%s' % step_run)
      else:
        file.write('%s\n' % step_run)
      file.write('\n')

    if use_hpc:
      file.close()
      # create a new case folder
      os.system('mkdir -p _work/%s' % fsub_name)

      # move submission file to the case folder
      if os.path.isfile('_work/%s/%s.sub' % (fsub_name, fsub_name)):
        os.remove('_work/%s/%s.sub' % (fsub_name, fsub_name))
      shutil.move(fsub_name + '.sub', '_work/%s/' % fsub_name)

      # link scripts
      links = glob.glob('*.py') + glob.glob('*.ex') + glob.glob('*.inp') + glob.glob('*.tmp')
      for fname in links:
        os.system('ln -sf %s _work/%s/' % (os.path.realpath(fname), fsub_name))

      # write submission script
      file = original_file
      file.write('## submit to hpc\n')
      file.write('cd _work/%s\n' % fsub_name)
      file.write('sbatch %s.sub\n' % fsub_name)
      file.write('cd ../../\n\n')

  return len(sub_jobs_list)

def write_all_jobs(file, job, detail, major, hpc = False):
  keys = list(detail['strategy']['vector'].keys())
  if len(keys) == 0:
    nsub_jobs_vector = 0
  else:
    for key in keys:
      if detail['strategy']['vector'][key][-4:] == ".txt":
        var = os.popen('cat %s' % detail['strategy']['vector'][key]).read()
        detail['strategy']['vector'][key] = list(map(str, var.split()))[::4]
    nsub_jobs_vector = len(detail['strategy']['vector'][keys[0]])

  njobs = 0
  # loop over vector jobs
  if nsub_jobs_vector > 0:
    for minor in range(nsub_jobs_vector):
      # set environment variable from vector
      for key, value in detail['strategy']['vector'].items():
        detail['env'][key] = value[minor]
      njobs += write_job_steps(file, job, detail, major, minor, hpc = hpc)
  else:
    njobs += write_job_steps(file, job, detail, major, 0, hpc = hpc)

  return njobs

if __name__ == '__main__':
  yml_file = args['input']
  wrk_file = '%s.sh' % yml_file[:-4]

  with open(yml_file, 'r') as file:
    try:
      work = yaml.safe_load(file)
    except yaml.YAMLError as exc:
      print(exc)

  with open(wrk_file, 'w') as file:
    file.write('#! /usr/bin/bash\n\n')
    file.write('# This file is automatically generated by %s ' % os.path.basename(__file__))
    file.write('based on %s\n\n' %  yml_file)
    finished_jobs, count, major, njobs = [], 0, 0, 0
    while len(finished_jobs) < len(work['jobs']):
      for job, detail in work['jobs'].items():
        # job has been finished
        if job in finished_jobs:
          continue

        # job was skipped
        if ('skip' in detail) and detail['skip']:
          finished_jobs.append(job)
          continue

        # job has dependency
        if 'needs' in detail:
          dependence_clear = True
          if not isinstance(detail['needs'], list):
            detail['needs'] = [detail['needs']]
          for dep in detail['needs']:
            if dep not in finished_jobs:
              dependence_clear = False
          if not dependence_clear:
            continue

        # dependency clear
        major += 1
        if 'env' not in detail:
          detail['env'] = []

        # write job vector and matrix
        if 'strategy' not in detail:
          detail['strategy'] = {}
          detail['strategy']['vector'] = {}
          detail['strategy']['matrix'] = {}
        else :
          if 'vector' not in detail['strategy']:
            detail['strategy']['vector'] = {}
          if 'matrix' not in detail['strategy']:
            detail['strategy']['matrix'] = {}

        # write matrix grid
        for key in detail['strategy']['matrix']:
          with open('%s-%s-grid.txt' % (job, key), 'w') as fgrid:
            for x in detail['strategy']['matrix'][key]:
              fgrid.write('%s ' % x)

        njobs += write_all_jobs(file, job, detail, major, hpc = args['hpc'])
        finished_jobs.append(job)

      count += 1
      if count > 10:
        raise RuntimeError("Loops count exceeds 1000")
    file.write('# Total number of jobs = %d\n' % njobs)
  print('workflow written to %s' % wrk_file)

  os.system('chmod +x %s' % wrk_file)
