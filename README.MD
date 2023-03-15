# Parallel code using MPI to perform relaxation technique 

Before running parallel or sequential program, first run createRandomFile.c
It will generate 'randNumbers' which are randomly generated integers written as binary data to the file
It's currently set to allow runs up to n=10,000

mpicc & mpirun to run parallel 
gcc & ./ to run sequential 
