#! /usr/bin/env python3
import re
from io import StringIO

# match file name
xstart = re.compile(r'^@@ ([^ ]+) -> (.+$)')

# match the action header
xheader = re.compile(
  r'\[([adr]):([+-]?\d),(first\b|second\b|third\b|fourth\b|last3\b|last2\b|last\b|all\b) *~ (.+$)'
  )

# match the replacement
xreplace = re.compile(r'([^ ]+) -> (.+$)')

# match end of the action
xfinish = re.compile(r'\];')

def ErrorMsg(msg):
  print("ERROR: " + msg)
  exit(1)

def LogMsg(msg, em = ''):
  if em != '':
    print(em + ' ' + msg + ' ' + em)
  else:
    print(msg)

def TakeAction(action, steps, src_lines, idx):
  if action == 'r':
    code = src_lines[idx]
    for step in steps:
      match = xreplace.search(step)
      if not match:
        ErrorMsg('Failed to parse action: %s' % step)
      old_str = match.group(1)
      new_str = match.group(2)
      old_code = code
      code = code.replace(old_str, new_str)
      if code != old_code:
        LogMsg('Replacement done at line #%d' % (idx+1,), em = '====')
        LogMsg('Old code:\n%s' % old_code[:-1])
        LogMsg('New code:\n%s' % code)
      else:
        ErrorMsg('Replacement failed.\nCannot find "%s" in "%s"' % (old_str, old_code[:-1]))
    src_lines[idx] = code
  elif action == 'a':
    old_str = src_lines[idx-1]
    new_str = ''
    for step in steps:
      src_lines[idx-1] += step
      new_str += step
    LogMsg('Added the following lines of code at line #%d' % idx, em = '====')
    LogMsg('After "%s":' % old_str.lstrip().rstrip())
    LogMsg(new_str)

def UpgradeFile(src_fname, dst_fname, rules):
  src_lines = open(src_fname, 'r').readlines()
  rules = StringIO(rules)
  rule = rules.readline()

  while rule:
    # action to take
    find_header = xheader.search(rule)
    while rule and (not find_header):
      rule = rules.readline()
      find_header = xheader.search(rule)
    if not rule: break

    action = find_header.group(1)
    offset = int(find_header.group(2))
    kind = find_header.group(3)
    position = find_header.group(4).strip()

    if kind == 'first':
      ncycle = 1
    elif kind == 'second':
      ncycle = 2
    elif kind == 'thrid':
      ncycle = 3
    elif kind == 'fourth':
      ncycle = 4
    elif kind == 'last3':
      ncycle = -3
    elif kind == 'last2':
      ncycle = -2
    elif kind == 'last':
      ncycle = -1 
    else: # all
      ncycle = 0

    # read all steps
    steps, max_steps, count = [], 1000, 0
    step = rules.readline()
    while not xfinish.search(step) and count < max_steps:
      steps.append(step)
      step = rules.readline()
      count += 1
    if count >= max_steps: 
      ErrorMsg('Missing end-of-action token "];"')

    # search for position in file
    cycle = 1
    ids = []
    for i,line in enumerate(src_lines):
      if position in line:
        ids.append(i)
        if cycle == ncycle: break
        cycle += 1

    if len(ids) == 0:
      ErrorMsg('Cannot find %s in %s' % (position, src_fname))

    # take action
    if ncycle == 0:
      for i in ids:
        TakeAction(action, steps, src_lines, i)
    elif ncycle > 0:
      TakeAction(action, steps, src_lines, ids[-1] + offset)
    else: # ncycle < 0
      TakeAction(action, steps, src_lines, ids[ncycle] + offset)

    # next rule
    rule = rules.readline()

  # write to destination file
  with open(dst_fname, 'w') as dst_file:
    for line in src_lines:
      dst_file.write(line)

if __name__ == '__main__':
  # start upgrading a file
  with open('upgrade.rule', 'r') as file:
    find_file = xstart.match(file.readline())
    while find_file:
      src_fname = find_file.group(1)
      dst_fname = find_file.group(2)
      rules = ''
      line = file.readline()
      while line and (not xstart.match(line)):
        rules += line
        line = file.readline()
      LogMsg('... Upgrading File %s ...' % dst_fname)
      UpgradeFile(src_fname, dst_fname, rules)
      find_file = xstart.match(line)
