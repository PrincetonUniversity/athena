# Regression test script

# Python modules
import argparse
import pkgutil
import sys

# Prevent generation of .pyc files
sys.dont_write_bytecode = True

# Athena++ modules
import scripts.utils.athena as athena

# Main function
def main(**kwargs):

  # Save existing files
  athena.save_files()

  # Make list of tests to run
  if kwargs['single'] is None:
    test_names = [name for _,name,_ in pkgutil.iter_modules(['scripts/tests'])]
  else:
    test_names = [kwargs['single']]

  # Run tests
  test_results = []
  try:
    for name in test_names:
      name_full = 'scripts.tests.' + name
      module = __import__(name_full, globals(), locals(), fromlist=['run_test'])
      test_results.append(module.run_test())

  # Restore any previously-existing files once runs are complete
  finally:
    athena.restore_files()

  # Report test results
  print('')
  print('Results:')
  for name,result in zip(test_names,test_results):
    result_string = 'passed' if result else 'failed'
    print('    {0}: {1}'.format(name,result_string))
  print('')
  num_tests = len(test_results)
  num_passed = test_results.count(True)
  test_string = 'test' if num_tests == 1 else 'tests'
  print('Summary: {0} out of {1} {2} passed'.format(num_passed,num_tests,test_string))
  print('')

# Execute main function
if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-s', '--single',
      type=str,
      default=None,
      help='name of individual test to run')
  args = parser.parse_args()
  main(**vars(args))
