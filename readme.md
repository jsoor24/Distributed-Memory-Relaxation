# Parallel code using MPI to perform relaxation technique 

This project solves a system of equations in a matrix through the 'relaxation technique'. Each value is updated by taking the average of the 4 values to the north/east/south/west of the target value. This can be done in parallel. 

## Running the code 

Before running parallel or sequential program, first run createRandomFile.c
It will generate 'randNumbers' which are randomly generated integers written as binary data to the file
It's currently set to allow runs up to n=10,000

mpicc & mpirun to run parallel 
gcc & ./ to run sequential 


## Approach 

To solve a n x n matrix with p processes:
1. Split the matrix into p sections where each section has n / p rows (integer division).
If mod = n % p != 0 then add 1 (row) to each of the first ’mod’ rows.
11 / 4 = 2.
11 % 4 = 3.
So each process gets 2 rows but the first 3 process get +1 row
2. Each process calculates their respective offset and reads in from the file
as required
3. Each process calls MPI_Irecv to receive the kill signal from ’above’ or
’below’ when they are done
{LOOP BEGIN}
4. Each process will use MPI_Isend and MPI_Irecv to send and receive
’shared’ rows that are outside of their section
5. Data is copied into ’copy’ array so that data is not overwritten whilst still
calculating
6. Each process computes their ’middle’ rows (averaging north, east, west,
south values)
If a section has 5 rows, the process will compute rows [1, 2, 3]
Rows [0, 4] require the information from ’above’ or ’below’
7. MPI_Wait is called
8. Once it has the data from above/below, a process can compute it’s top/bot-
tom row
9. Wait to ensure all data has been sent to other rows before starting next
iteration
10. If all the computations for a section are below the given precision, the
process is ’done’
{LOOP ENDS}
11. Once a process has finished, it sends a kill signal to ’above’ and ’below’ to
indicate that is has finished
12. The process must also send it’s final values for top/bottom rows to above/-
below so other processes have correct final values


The reasoning for this approach was to limit the amount of data sent over
message passing as this would be an overhead that could slow the program.
Additionally, MPI_Isend/MPI_Irecv are used to reduce the amount of time a
process is blocked waiting for communication. The process will continue to
do independent computations such as copying data or computing middle rows
instead of idly waiting.
