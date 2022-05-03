#! /usr/bin/env python3
import re, sympy, io

def get_derivative(rstr, v):
  if '?' in rstr:  # has condition
    m = re.search('([^\?]+)\?(.+):(.+)', rstr)
    condition = m[1].strip()
    branch1 = m[2].strip()
    branch2 = m[3].strip()
    dstr1 = get_derivative(branch1, v)
    dstr2 = get_derivative(branch2, v)
    if dstr1 == '0' and dstr2 == '0':
      return '0'
    else:
      return '%s ? %s : %s' % (condition, dstr1, dstr2)
  else:
    return str(sympy.diff(rstr, v))

def split_stoichiometry(var):
  result = re.match('^[1-9][0-9]*', var)
  if (result):
    c = int(result.group())
    var = re.sub('^[1-9][0-9]*', '', var)
  else:
    c = 1
  return var, int(c)

def parse_reaction_formula(formula):
  formula = formula.split('->')
  inp = formula[0].split('+')
  inp = [x.strip() for x in inp]
  out = formula[1].split('+')
  out = [x.strip() for x in out]
  return inp, out

def parse_full_reaction(text):
  rmap = {}
  uncondition = text.split('|')
  if len(uncondition) > 1:
    rmap['condition'] = uncondition[1].strip()
  else:
    rmap['condition'] = None

  uncondition = uncondition[0]
  fields = uncondition.split(';')
  rmap['formula'] = fields[0].strip()
  if len(fields) > 1:
    rmap['rate'] = fields[1].strip()
  #if len(fields) > 2:
  #  rmap['enthalpy'] = fields[2].strip()
  if len(fields) > 2:
    raise("ERROR")

  return rmap

def add_prefix(lines, s):
  lines = io.StringIO(lines)
  newlines = ''
  for line in lines:
    newlines += s + line
  return newlines

def add_surfix(lines, s):
  lines = io.StringIO(lines)
  newlines = ''
  for line in lines:
    newlines += line[:-1] + s + line[-1]
  return newlines

def read_variables(fname):
  variables, inblock = [], False
  with open(fname, 'r') as file:
    lines = file.readlines()
    for line in lines:
      line = line.strip()
      if line == '# variables':
        inblock = True
        continue
      if inblock:
        if len(line) == 0: break
        if line[0] == '-':
          variables += [line.split('-')[1].strip()]
  return variables

def read_relations(fname):
  relations, inblock = {}, False
  with open(fname, 'r') as file:
    lines = file.readlines()
    for line in lines:
      line = line.strip()
      if line == '# relations':
        inblock = True
        continue
      if inblock:
        if len(line) == 0: break
        if line[0] == '-':
          fields = line[1:].split('->')
          name, replace = fields[0].strip(), fields[1].strip()
          name = re.sub(' ', '', name)
          relations[name] = replace
  return relations

def read_reactions(fname, rtype):
  reactions, inblock = [], False
  with open(fname, 'r') as file:
    lines = file.readlines()
    for line in lines:
      line = line.strip()
      if line == '# %s reactions' % rtype:
        inblock = True
        continue
      if inblock:
        if len(line) == 0: break
        if line[0] == '-':
          reactions += [line[1:].strip()]
  return reactions

def read_coefficients(fname):
  coeffs, inblock = [], False
  with open(fname, 'r') as file:
    lines = file.readlines()
    for line in lines:
      line = line.strip()
      if line == '# coefficients':
        inblock = True
        continue
      if inblock:
        if len(line) == 0: break
        if line[0] == '-':
          coeffs += [line[1:].strip()]
  return coeffs 

def read_verbatim(fname):
  texts, inblock, inverb = [], False, False
  with open(fname, 'r') as file:
    lines = file.readlines()
    for line in lines:
      line = line.strip()
      if line == '# verbatim':
        inblock = True
        continue
      if inblock:
        if line[:3] == '~~~':
          inverb = not inverb;
          if inverb: continue
          else: break
        texts += [line]
  return '\n'.join(texts) + '\n'

def write_definitions(symbols):
  output = ''
  for key, value in symbols.items():
    if key != 'T':
      output += 'Real q%s = %s;\n' % (get_legal_name(key), value)
    else:
      output += 'Real T = gdata.q[IDN];\n'
  return add_prefix(output, '  ')

def write_index_comments(symbols):
  slist = list(symbols.keys())
  output = '// 0 : T\n'
  for i,key in enumerate(slist):
    if key != 'T':
      output += '// %d : %s\n' % (i, key)
  return output

def write_coefficients(coeffs):
  return add_surfix(add_prefix('\n'.join(coeffs)+'\n', 'Real '), ';')

def write_reaction_string(rstr, istoi, ostoi, slist):
  output = ''
  if '?' in rstr: # has condition
    m = re.search('([^\?]+)\?(.+):(.+)', rstr)
    condition = m[1].strip()
    branch1 = m[2].strip()
    branch2 = m[3].strip()
    output += 'if (%s) {\n' % condition
    output += write_reaction_string(branch1, istoi, ostoi, slist)
    if branch2 == '0':
      output += '}\n'
    else:
      output += '} else {\n'
      output += write_reaction_string(branch2, istoi, ostoi, slist)
      output += '}\n'
    return add_prefix(output, '  ')
  else: # no condition
    vlist = ['T'] + ['q%d' % (j+1) for j in range(len(istoi)+len(ostoi))]

    # reactants
    output += '// reactants\n'
    enthalpy, i = [], 1
    for name,c in istoi:
      ix = slist.index(name)
      if c == 0: continue
      if c == 1:
        rate_line = 'rate[%d] -= %s;\n' % (ix, rstr)
        enthalpy += ['H%d' % i]
      else:
        rate_line = 'rate[%d] -= %d*(%s);\n' % (ix, c, rstr)
        enthalpy += ['%d*H%d' % (c,i)]
      output += rate_line

      # jacobian
      for v in vlist: 
        jstr = str(sympy.diff(rstr, v))
        if jstr == '0': continue
        output += 'jac[%d][%d] -= %s;\n' % (ix, vlist.index(v), jstr)
      i += 1

    enthalpy_change = '+'.join(enthalpy)
    output += '\n'

    # resultants
    output += '// resultants\n'
    enthalpy = []
    for name,c in ostoi:
      ix = slist.index(name)
      if c == 0: continue
      if c == 1:
        rate_line = 'rate[%d] += %s;\n' % (ix, rstr)
        enthalpy += ['H%d' % i]
      else:
        rate_line = 'rate[%d] += %d*(%s);\n' % (ix, c, rstr)
        enthalpy += ['%d*H%d' % (c,i)]
      output += rate_line

      # jacobian
      for v in vlist: 
        jstr = str(sympy.diff(rstr, v))
        if jstr == '0': continue
        output += 'jac[%d][%d] += %s;\n' % (ix, vlist.index(v), jstr)
      i += 1

    enthalpy_change += ' - ' + ' - '.join(enthalpy)
    output += '\n'

    # enthalpy change
    output += '// enthalpy change\n'
    output += 'rate[%d] += (%s)*(%s);\n' % (vlist.index('T'), rstr, enthalpy_change)
    for v in vlist: 
      jstr = str(sympy.diff(rstr, v))
      if jstr == '0': continue
      output += 'jac[%d][%d] += (%s)*(%s);\n' %  \
        (vlist.index('T'), vlist.index(v), jstr, enthalpy_change)

    return add_prefix(output, '  ')

def get_legal_name(name):
  name = name.replace('(','')
  name = name.replace(')','')
  name = name.replace(',','')
  name = name.replace('-','_')
  return name

def write_reaction_function(fname, reaction, symbols, functions, count):
  space = '  '
  # formal reaction string
  formula = reaction['formula'].strip()
  reaction_split = reaction['function'].split(':')
  short_name = reaction_split[0]
  if len(reaction_split) > 1:
    coefficients = reaction_split[1].split(',')
  else:
    coefficients = []
  slist = list(symbols.keys())
  print('working on reaction %s' % formula)

  # check that all parameters are defined
  myfunc = functions[short_name]
  for x in myfunc['parameters']:
    if '=' in x:
      error = True
      for coeff in coefficients:
        if x in coeff:
          error = False
      if error:
        raise ValueError('Parameter %s undefiend in reaction function %s' % (x[:-1], short_name))

  inps, outs = parse_reaction_formula(formula)
  # istoi/ostoi is a tuple of (name, count)
  istoi = [list(split_stoichiometry(x)) for x in inps]
  ostoi = [list(split_stoichiometry(x)) for x in outs]
  for x in istoi:
    for y in ostoi:
      if x[0] == y[0]:
        v = min(x[1], y[1])
        x[1] -= v
        y[1] -= v

  # reaction long name
  long_name = short_name + '_'
  long_name += '_'.join([x[0] for x in istoi])
  if '<->' in formula:
    long_name += '_eq_'
  else:
    long_name += '_to_'
  long_name += '_'.join([x[0] for x in ostoi])
  long_name = get_legal_name(long_name)

  # reaction function
  output = '// %s: %s\n' % (short_name, formula)
  arguments = ', '.join(['double q%d' % (j+1) for j in range(len(istoi)+len(ostoi))]) + ', '
  if len([x for x in myfunc['parameters'] if '=' not in x]) > 0:
    arguments += ', '.join(['double H%d' % (j+1) for j in range(len(istoi)+len(ostoi))]) + ', '
    arguments += ', '.join(['double %s' % x for x in myfunc['parameters'] if '=' not in x])
  else:
    arguments += ', '.join(['double H%d' % (j+1) for j in range(len(istoi)+len(ostoi))])
  output += 'inline void %s(double *rate, double **jac,\n%s) {\n' % (long_name, space + arguments)
  for coeff in coefficients:
    output += space + 'double %s;\n' % coeff
  output = output.replace('=', ' = ')

  # reaction rate string
  rstr = myfunc['rate'].strip()
  #print(istoi, ostoi)
  output += write_reaction_string(rstr, istoi, ostoi, slist)

  # substitution 
  if 'substitution' in myfunc:
    for rule in myfunc['substitution']:
      x, y = rule.split('->')
      output = output.replace(x.strip(), y.strip())

  # end reaction
  output += '}\n\n'

  with open('reactions_src/%s.hpp' % fname, 'a') as file:
    file.write(output)
    #print('reaction file written to reactions/%s.hpp' % long_name)

  # write function call
  arguments = arguments.replace('double ', '')
  vmap, hmap, n = {}, {}, 1
  for name,c in istoi:
    vmap['q%d' % n] = 'q' + get_legal_name(name)
    hmap['H%d' % n] = 'u_energy[%d]' % slist.index(name)
    n += 1
  for name,c in ostoi:
    vmap['q%d' % n] = 'q' + get_legal_name(name)
    hmap['H%d' % n] = 'u_energy[%d]' % slist.index(name)
    n += 1
  for i in range(1, n):
    arguments = arguments.replace('H%d,' % i, hmap['H%d' % i] + ',')
    arguments = arguments.replace(', H%d' % i, ', ' + hmap['H%d' % i])
    arguments = arguments.replace('q%d,' % i, vmap['q%d' % i] + ',')
    arguments = arguments.replace(', q%d' % i, ', ' + vmap['q%d' % i])
  output = space + '// R%d,%s: %s\n' % (count, short_name, formula)
  output += space + '%s(rate, jac, %s);\n' % (long_name, arguments)
  output += '\n'

  return output

def write_all_reactions(fname, data):
  nreactions = len(data['reactions'])
  output = ''
  # write header
  with open('reactions_src/%s.hpp' % fname, 'w') as file:
    file.write('/*! \\file %s.hpp\n' % fname)
    file.write(' *  \\brief reaction functions for %s system\n' % fname.replace('_', '-'))
    file.write(' */\n')
    file.write('#ifndef REACTIONS_%s\n' % fname.upper())
    file.write('#define REACTIONS_%s\n\n' % fname.upper())

  # write individual reactions
  for i in range(nreactions):
    output += write_reaction_function(fname, 
      data['reactions'][i], data['symbols'], data['functions'], i+1)
  print('reactions written to reactions_src/%s.hpp' % fname)

  # write footer
  with open('reactions_src/%s.hpp' % fname, 'a') as file:
    file.write('#endif\n')
  return output
