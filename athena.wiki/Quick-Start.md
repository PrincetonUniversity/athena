### Download the Code

You can download the latest public version on [[the Github repository|https://github.com/PrincetonUniversity/athena-public-version/]]. You can download the code from the repository using any Git client. To download the code using the standard git clone, type:

    > git clone https://github.com/PrincetonUniversity/athena-public-version

To retrieve the latest update, move into the code directory and type:

    > git pull

For more information, please consult the documentation of the Git client you use.



### Testing that the code runs

A quick way of testing that the code works is to run the regression test suite. From within the `athena/` directory, go to the regression directory and run the tests:

    > cd tst/regression
    > python run_tests.py

This may take several minutes, at the end of which a summary will be printed.

### Configuring and compiling

Athena++ is compiled from an automatically generated `Makefile`. From the top-level directory, the code can be configured and compiled as follows:

    > python configure.py <options>
    > make clean
    > make

This places the executable `athena` in the `bin/` directory. Available configuration options are covered in [[Configuring]].

### Running

The code is easily run by specifying an input file:

    > cd bin
    > ./athena -i <path/filename>

Additional command line options are covered in [[Running the Code]].