#! /usr/bin/env python3
import argparse, glob, os, yaml
from chemistry_parser import *

chemdir = 'reactions'
chem_choices = glob.glob('%s/*.yml' % chemdir)
chem_choices = [os.path.basename(x)[:-4] for x in chem_choices]

template_choices = glob.glob('templates/*.tmp')
template_choices = [os.path.basename(x)[:-4] for x in template_choices] + ['none']

parser = argparse.ArgumentParser()

parser.add_argument('-c', '--chemistry',
  choices = chem_choices,
  required = True,
  help = 'select chemistry')

parser.add_argument('-t', '--template',
  choices = template_choices,
  default = 'none',
  help = 'select template')

# parse command-line inputs
args = vars(parser.parse_args())

def merge_list(list1, list2, fname):
  for v in list2:
    if v not in list1:
      list1.append(v)

def merge_dictionary(dict1, dict2, fname):
  for key, value in dict2.items():
    if key in dict1:
      if isinstance(dict1[key], list):
        if isinstance(value, list):
          merge_list(dict1[key], value, fname+':'+key)
        else:
          msg = 'Conflict key "%s" in file %\n' % (key, fname)
          msg += '%s was a list but is not now' % key
          raise ValueError(msg)
      elif isinstance(dict1[key], dict):
        if isinstance(value, dict):
          merge_dictionary(dict1[key], value, fname+':'+key)
        else:
          msg = 'Conflict key "%s" in file %\n' % (key, fname)
          msg += '%s was a dictionary but is not now' % key
          raise ValueError(msg)
      else:
        if dict1[key] != value:
          msg = 'Conflict key "%s" in file %s\n' % (key, fname)
          msg += 'Value of "%s" was %s\n' % (key, dict1[key])
          msg += 'But is %s now\n' % value
          raise ValueError(msg)
    else:
      dict1[key] = value

def load_yaml_dict(fname, data):
  with open(fname, 'r') as file:
    try:
      data1 = yaml.safe_load(file)
    except yaml.YAMLError as exc:
      print(exc)

  if 'using' in data1:
    for fname in data1['using']:
      load_yaml_dict('%s/%s.yml' % (chemdir, fname), data)

  merge_dictionary(data, data1, fname)
  return data

def dump_yaml_dict(fname, data):
  with open(fname, 'w') as file:
    yaml.dump(data, file, default_flow_style = False)

if __name__ == '__main__':
  data = {}
  load_yaml_dict('%s/%s.yml' % (chemdir, args['chemistry']), data)

  freactions = args['chemistry'].replace('-','_')
  reaction_calls = write_all_reactions(freactions, data)
  definitions = write_definitions(data['symbols'])
  index_comment = write_index_comments(data['symbols'])

  ftemplate = args['template']
  if ftemplate != 'none':
    with open('templates/%s.tmp' % ftemplate, 'r') as file:
      tmp = file.read()
    tmp = re.sub('<header>', 'reactions_src/%s.hpp' % freactions, tmp)
    tmp = re.sub('<common>', add_prefix(data['common'], '  '), tmp)
    tmp = re.sub('<index_comments>', index_comment, tmp)
    tmp = re.sub('<definitions>', definitions, tmp)
    tmp = re.sub('<reaction_calls>', reaction_calls, tmp)

    with open('%s.cpp' % ftemplate, 'w') as file:
      file.write(tmp)
    print('%s chemistry written to %s.cpp' % (args['chemistry'], ftemplate))
