#! /usr/bin/env python3
import re, sympy, io

def get_simple_rate(rate, inps):
  if len(inps) == 1:
    return rate+'*'+inps[0]
  else:
    return rate+'*'+'*'.join(inps)

def get_derivative(rstr, s, vlist):
  dstr = str(sympy.diff(rstr, s))
  result = re.search('Derivative\(.+\)', dstr)

  if (result):
    x = result.group()
    x = re.sub(' ', '', x)
    dstr = re.sub('Derivative\(.+\)', '<sub>', dstr)
  #for v in vlist:
  #  dstr = re.sub(v, 'q[i'+v+']', dstr)
  if (result):
    dstr = re.sub('<sub>', x, dstr)
  return dstr

def split_stoichiometry(var):
  result = re.match('[1-9][0-9]*', var)
  if (result):
    c = result.group()
    var = re.sub(c, '', var)
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

def write_header(fname):
  with open(fname, 'r') as file:
    line = file.readline()
  result = re.search('# \*\*(.+)\*\*', line)
  name = result.group(1)
  header = '''/** @file %s_impl.hpp
 * @brief Implement %s Chemistry
 *
 * This file is automatically generated by make_chemistry.py
 *
 * @author Cheng Li
 * @bug No know bugs.
 */

template<typename D1, typename D2>
void %s::AssembleReactionMatrix(Eigen::DenseBase<D1>& rate,
  Eigen::DenseBase<D2>& jac, Real const q[], Real cv, Real time)
{
  Thermodynamics *pthermo = pmy_chem->pmy_block->pthermo;
'''
  header = header % (name.lower(), name, name)
  return header

def write_footer():
  footer = '''}'''
  return footer

def write_definitions(vlist):
  output = ''
  for i,v in enumerate(vlist):
    output += 'int i%s = index_[%d];\n' % (v,i)
    output += 'Real %s = q[i%s];\n' % (v,v)
  return output

def write_coefficients(coeffs):
  return add_surfix(add_prefix('\n'.join(coeffs)+'\n', 'Real '), ';')

def write_reaction(rmap, vlist, rtype, relations = None):
  formula = rmap['formula']
  rate = rmap['rate']
  inps, outs = parse_reaction_formula(formula)
  istoi = [list(split_stoichiometry(x)) for x in inps]
  ostoi = [list(split_stoichiometry(x)) for x in outs]
  for x in istoi:
    for y in ostoi:
      if x[0] == y[0]:
        v = min(x[1], y[1])
        x[1] -= v
        y[1] -= v

  if rtype == 'simple':
    rstr = get_simple_rate(rate, inps)
  else:
    rstr = rate
  syms = sympy.symbols(vlist)

  #rstr2 = re.sub('([a-zA-z]+\w*)\([^()]+\)', '\g<1>', rstr)
  #for v in vlist:
  #  rstr2 = re.sub(v, 'q[i'+v+']', rstr2)

  output = '// %s; %s\n' % (formula, rate)

  # reactants
  enthalpy = []
  for name,c in istoi:
    if c == 0: continue
    if c == 1:
      rate_line = 'rate(%d) -= %s;\n' % (vlist.index(name), rstr)
      enthalpy += ['deltaU_[i%s]' % name]
    else:
      rate_line = 'rate(%d) -= %d*%s;\n' % (vlist.index(name), c, rstr)
      enthalpy += ['%d*deltaU_[i%s]' % (c,name)]
    output += rate_line
    for s in syms:
      jstr = get_derivative(rstr, s, vlist)
      if jstr == '0': continue
      output += 'jac(%d,%d) -= %s;\n' % (vlist.index(name), vlist.index(str(s)), jstr)

  enthalpy_change = '+'.join(enthalpy)
  output += '\n'

  # resultants
  enthalpy = []
  for name,c in ostoi:
    if c == 0: continue
    if c == 1:
      rate_line = 'rate(%d) += %s;\n' % (vlist.index(name), rstr)
      enthalpy += ['deltaU_[i%s]' % name]
    else:
      rate_line = 'rate(%d) += %d*%s;\n' % (vlist.index(name), c, rstr)
      enthalpy += ['%d*deltaU_[i%s]' % (c,name)]
    output += rate_line
    for s in syms:
      jstr = get_derivative(rstr, s, vlist)
      if jstr == '0': continue
      output += 'jac(%d,%d) += %s;\n' % (vlist.index(name), vlist.index(str(s)), jstr)

  enthalpy_change += ' - ' + ' - '.join(enthalpy)
  output += '\n'

  # enthalpy change
  output += '// enthalpy change\n'
  output += 'rate(%d) += (%s)*(%s)/cv;\n' % (vlist.index('T'), rstr, enthalpy_change)
  for s in syms:
    jstr = get_derivative(rstr, s, vlist)
    if jstr == '0': continue
    output += 'jac(%d,%d) += (%s)*(%s)/cv;\n' %  \
      (vlist.index('T'), vlist.index(str(s)), jstr, enthalpy_change)

  # replace relations
  if relations != None:
    for key in relations:
      output = output.replace(key, relations[key])

  return output

def write_all_reactions(fname, vlist, rtypes):
  output, cond, cond_str = '', {}, {}
  relations = read_relations(fname)

  for rtype in rtypes:
    reactions = read_reactions(fname, rtype)
    for r in reactions:
      rmap = parse_full_reaction(r)
      if rmap['condition'] != None:
        #output = 'if (%s) {\n' % rmap['condition']
        key = rmap['condition'].replace(' ', '')
        if key in cond.keys():
          cond[key] += '\n' + write_reaction(rmap, vlist, rtype, relations = relations)
        else: 
          cond[key] = write_reaction(rmap, vlist, rtype, relations = relations)
          cond_str[key] = rmap['condition']
        #output += '}\n'
      else:
        rstr = write_reaction(rmap, vlist, rtype, relations = relations)
        output += rstr + '\n'

  # print conditioned reactions
  for key in cond.keys():
    output += 'if (%s) {\n' % cond_str[key]
    output += add_prefix(cond[key], '  ')
    output += '}\n\n'

  return output[:-1]
