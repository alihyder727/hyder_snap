#! /usr/bin/env python3
import fnmatch, sys
import os, re
import subprocess

# recursively obtain all source files in a directory
def SourceFiles(directory, pattern):
  matches = []
  for root, dirnames, filenames in os.walk(directory):
    for filename in fnmatch.filter(filenames, pattern):
      if not 'z.junk' in root:
        matches.append(os.path.join(root, filename))
  return matches

# generate patch file
def PatchStructure(patch):
  # athena source files
  files = SourceFiles('src', '*.?pp')
  athena = [f[4:] for f in files]

  replacement, addition = [], []
  files = SourceFiles(patch, '*.?pp')
  for fname in files:
    source = fname[len(patch)+1:]
    if source in athena:
      replacement.append(source)
    else:
      addition.append(source)
  return replacement, addition

# remove all softlinks
def CleanSoftlinks(directory):
  os.system('find -L %s -xtype l -delete' % directory)
  files = SourceFiles('src', '*.?pp.old')
  for f in files:
    if not os.path.isfile(f[:-4]):
      os.rename(f, f[:-4])
    else:
      raise SystemExit('File %s exists' % f[:-4])

# write patch file
def WritePatchFile(fname, patches):
  CleanSoftlinks('src')
  if os.path.isfile(fname):
    os.remove(fname)
  for patch in patches:
    replacement, addition = PatchStructure(patch)
    with open(fname, 'a') as file:
      file.write('# Replacement files from package %s:\n' % patch)
      for f in replacement:
        file.write('%s/%s -> src/%s\n' % (patch, f, f))
      file.write('\n# Additional files from package %s:\n' % patch)
      for f in addition:
        file.write('%s/%s\n' % (patch, f))
      file.write('\n')

nargs = len(sys.argv)
argv = ['drum']
if nargs > 1:
    argv.extend(sys.argv[1:])
WritePatchFile('patch_files', argv)
print('Patch rules written to "patch_files" using package(s): %s.' % ', '.join(argv))
