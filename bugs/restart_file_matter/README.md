Restart file bug

Bug found by: FZ

Explanation: Segmentation fault occurs when one tries to generate restart file
using the input file here. It occurs when, at the same time, num_threads = 1,
restart file is enabled and athena is run as follows:

./athena -i athinput.twopunctures time/nlim=0

with the purpose of creating a restart file for initial conditions. Never tried 
to see if it happens regardless of nlim=0.

This error occurs also in the old codebase. The problem is exactly located in 
src/outputs/restart.cpp
at line 219, when trying to copy "mat" data in z4c.

The reason is not clear. The error DOES NOT occur if, i.e., num_threads>=2.
