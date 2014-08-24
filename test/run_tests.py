#!/usr/bin/python

# Runs Google Test suites and summarizes results

# Modules
import argparse
import glob
import os
import xml.etree.ElementTree as et

# Main function
def main(**kwargs):

  # Delete old logs
  print('Deleting old test logs...')
  try:
    os.system('rm -f logs/*.xml')
  except OSError as err:
    print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
    exit()

  # Run tests
  print('Running tests...')
  executables = glob.glob('bin/run_tests_*')
  for executable in executables:
    command = executable + ' --gtest_output=xml:logs/'
    if kwargs['quiet']:
      command += ' > /dev/null'
    try:
      os.system(command)
    except OSError as err:
      print('OS Error ({0}): {1}'.format(err.errno, err.strerror))
      exit()

  # Parse results
  print('Parsing test results...')
  test_count = 0
  failure_count = 0
  logs = glob.glob('logs/run_tests_*.xml')
  for log in logs:
    tree = et.parse(log)
    testsuites = tree.getroot()
    for testsuite in testsuites:
      for testcase in testsuite:
        test_count += 1
        if testcase.find('failure') is not None:
          if kwargs['quiet']:
            print('\033[1;31mfailure:\033[0;0m {0}.{1}'.format(testsuite.attrib['name'], testcase.attrib['name']))
          failure_count += 1

  # Summarize results
  print('Testing completed:')
  print('  {0} test{1} run'.format(test_count, 's' if test_count != 1 else ''))
  if failure_count == 0:
    print('  \033[1;32mAll tests passed\033[0;0m')
  else:
    print('  \033[1;31m{0} test{1} failed\033[0;0m'.format(failure_count, 's' if failure_count > 1 else ''))

if __name__ == '__main__':
  parser = argparse.ArgumentParser()
  parser.add_argument('-q', '--quiet', action='store_true', default=False,
      help='flag indicating live individual test output is to be suppressed')
  args = parser.parse_args()
  main(**vars(args))
