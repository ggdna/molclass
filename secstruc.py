import math
import subprocess as sub
import os

def setup_args(**kargs):
  """ Returns the list of tokens specifying the command to be run in the pipe. """
  args = [kargs['exec_name'],
          '-material', kargs['material'],   '-sodium', kargs['sodium'],
          '-magnesium', kargs['magnesium'], '-dangles', kargs['dangles'], '-T', kargs['T']]
  if kargs['multi']: args += ['-multi']
  if kargs['pseudo']: args += ['-pseudo']
  return args

def setup_cmd_input(multi, sequences, ordering, structure = ''):
  """ Returns the command-line input string to be given to NUPACK. """
  if not multi:
    cmd_input = '+'.join(sequences) + '\n' + structure
  else:
    n_seqs = len(sequences)
    if ordering == None:
      seq_order = ' '.join([str(i) for i in range(1, n_seqs+1)])
    else:
      seq_order = ' '.join([str(i) for i in ordering])
    cmd_input = str(n_seqs) + '\n' + ('\n'.join(sequences)) + '\n' + seq_order + '\n' + structure
  return cmd_input.strip()

def setup_nupack_input(**kargs):
  """ Returns the list of tokens specifying the command to be run in the pipe, and
  the command-line input to be given to NUPACK.
  Note that individual functions below may modify args or cmd_input depending on their
  specific usage specification. """
  # Set up terms of command-line executable call
  args = setup_args(**kargs)

  # Set up command-line input to NUPACK
  cmd_input = setup_cmd_input(kargs['multi'], kargs['sequences'], kargs['ordering'],
                              kargs.get('structure', ''))

  return (args, cmd_input)

def call_with_file(args, cmd_input, outsuffix):
  """ Performs a NUPACK call, returning the lines of the output in a temporary
  output file. The output file is assumed to have the suffix 'outsuffix'.
  outsuffix includes the period (.) delimiter.
    Ex:
      call_with_file(args, input, '.sample')
  """

  import tempfile

  ## Preliminaries
  # Set up temporary output file
  outfile = tempfile.NamedTemporaryFile(delete=False, suffix=outsuffix)
  outprefix = outfile.name[:-len(outsuffix)]

  # Close the output file so sample can open/write to it.
  # Will reopen it later to get the output.
  outfile.close()

  ## Perform executable call, ignoring pipe output
  args = [str(s) for s in args] # all argument elements must be strings
  cmd_input = outprefix + '\n' + cmd_input # prepend the output file prefix to the input for NUPACK
  p = sub.Popen(args, stdin=sub.PIPE, stdout=sub.PIPE, stderr=sub.STDOUT)
  p.communicate(cmd_input.encode('utf-8'))

  ## Process and return output
  # Read output file and clean it up
  # Note that it was created by us, so it won't be cleaned up automatically
  out = open(outfile.name, "rt")
  output_lines = out.readlines()
  out.close()
  os.remove(outfile.name)

  return output_lines

def pairs(sequences, ordering = None, material = 'rna',
          dangles = 'some', T = 37, multi = True, pseudo = False,
          sodium = 1.0, magnesium = 0.0, cutoff = 0.001):
  """Calls NUPACK's pairs executable on a complex consisting of the unique strands in sequences.
     Returns the probabilities of pairs of bases being bound, only including those pairs
     with probability greater than cutoff.
       sequences is a list of the strand sequences
       See NUPACK User Manual for information on other arguments.
  """

  ## Set up command-line arguments and input
  args, cmd_input = \
    setup_nupack_input(exec_name = 'pairs', sequences = sequences, ordering = ordering,
                       material = material, sodium = sodium, magnesium = magnesium,
                       dangles = dangles, T = T, multi = multi, pseudo = pseudo)
  if multi:
    suffix = '.epairs'
  else:
    suffix = '.ppairs'

  ## Perform call
  output = call_with_file(args, cmd_input, suffix)

  ## Parse and return output
  pair_probs = []
  for l in filter(lambda x: x[0].isdigit(), output):
    if len(l.split()) > 1:
      i,j,p = l.split()
      pair_probs.append(tuple( (int(i),int(j),float(p)) ))

  return pair_probs
