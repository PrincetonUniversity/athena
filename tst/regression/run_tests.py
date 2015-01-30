"""
Regression test script.

Usage: From this directory, call this script with python:
      python run_tests.py

Notes:
  - Requires Python 2.7+.
  - This file should not be modified when adding new scripts.
  - To add a new script, create a new .py file in a scripts/tests/ subdirectory.
  - See scripts/tests/example.py for an example.
  - For more information, check online regression test documentation.
"""

# Python modules
import argparse
import os
import pkgutil
import sys
import traceback

# Prevent generation of .pyc files
sys.dont_write_bytecode = True

# Athena++ modules
import scripts.utils.athena as athena

# Main function
def main(**kwargs):

  # Save existing files
  athena.save_files()

  # Make list of tests to run
  test_names = []
  if kwargs['single'] is not None:
    test_names.extend([name.replace('/','.') for name in kwargs['single']])
  if kwargs['collection'] is not None:
    for directory in kwargs['collection']:
      if directory[-1] == '/':
        directory = directory[:-1]
      test_names.extend([name for _,name,_ in
          pkgutil.iter_modules(path=['scripts/tests/'+directory],
          prefix=directory+'.')])
  if kwargs['single'] is None and kwargs['collection'] is None:
    for _,directory,ispkg in pkgutil.iter_modules(path=['scripts/tests']):
      if ispkg:
        test_names.extend([name for _,name,_ in
            pkgutil.iter_modules(path=['scripts/tests/'+directory],
            prefix=directory+'.')])
  test_names = list(set(test_names))

  # Run tests
  current_dir = os.getcwd()
  test_results = []
  test_errors = []
  try:
    for name in test_names:
      try:
        name_full = 'scripts.tests.' + name
        module = __import__(name_full, globals(), locals(),
            fromlist=['prepare','run','analyze'])
        os.system('rm -rf {0}/bin'.format(current_dir))
        try:
          module.prepare()
        except Exception:
          traceback.print_exc()
          test_errors.append('prepare()')
          raise TestError(name_full.replace('.','/')+'.py')
        try:
          module.run()
        except Exception:
          traceback.print_exc()
          test_errors.append('run()')
          raise TestError(name_full.replace('.','/')+'.py')
        try:
          result = module.analyze()
        except Exception:
          traceback.print_exc()
          test_errors.append('analyze()')
          raise TestError(name_full.replace('.','/')+'.py')
      except TestError as err:
        test_results.append(False)
        print('---> Error in '+str(err))
      else:
        test_results.append(result)
        test_errors.append(None)
      finally:
        os.system('rm -rf {0}/bin'.format(current_dir))

  # Restore any previously-existing files once runs are complete
  finally:
    athena.restore_files()

  # Report test results
  print('')
  print('Results:')
  for name,result,error in zip(test_names,test_results,test_errors):
    result_string = 'passed' if result else 'failed'
    error_string = ' -- unexpected failure in {0} stage'.format(error) \
        if error is not None else ''
    print('    {0}: {1}{2}'.format(name,result_string,error_string))
  print('')
  num_tests = len(test_results)
  num_passed = test_results.count(True)
  test_string = 'test' if num_tests == 1 else 'tests'
  print('Summary: {0} out of {1} {2} passed'.format(num_passed,num_tests,test_string))
  print('')

# Exception for unexpected behavior by individual tests
class TestError(RuntimeError):
  pass

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-s', '--single',
      type=str,
      default=None,
      nargs='+',
      help='name of individual test to run, relative to scripts/tests/, excluding .py')
  parser.add_argument('-c', '--collection',
      type=str,
      default=None,
      nargs='+',
      help='name of collection (directory) of tests to run, relative to scripts/tests/')
  args = parser.parse_args()
  main(**vars(args))
