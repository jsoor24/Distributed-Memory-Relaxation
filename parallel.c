#include <errno.h>
#include <math.h>
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

void get_time(struct timespec *time);
double elapsed_time(struct timespec start, struct timespec stop);
void open_file(MPI_File *fr, int rank, char *filename);
void check_file(MPI_File fr, int n);
int get(int i, int j, int n);
void print_data(int rank, int world_size, int rows, int n, double *data);
double average(int i, int j, int n, int rows, double *data, double *above,
               double *below);
int compute_row(int i, int rows, int n, double *data, double *copy,
                double *above, double *below, double precision);
void relaxation(double *data, int rank, int world_size, int rows, int n,
                double precision);

void get_time(struct timespec *time) {
    if (clock_gettime(CLOCK_REALTIME, time) == -1) {
        printf("Error getting clock time");
        exit(6);
    }
}

double elapsed_time(struct timespec start, struct timespec stop) {
    return (double)(stop.tv_sec - start.tv_sec) +
           (double)(stop.tv_nsec - start.tv_nsec) / (double)1000000000L;
}

void open_file(MPI_File *fr, int rank, char *filename) {
    int err = MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY,
                            MPI_INFO_NULL, fr);
    if (err) {
        if (rank == 0) {
            printf("Error opening file for fr\n");
        }
        exit(1);
    }
}

// Checks the value of N for the file and exits if it's too small
void check_file(MPI_File fr, int n) {
    int f_n;
    int err = MPI_File_read(fr, &f_n, 1, MPI_INT, MPI_STATUS_IGNORE);
    if (err) {
        printf("Error reading file line 1\n");
        exit(2);
    }

    if (n > f_n) {
        printf("Provided n: %d too large. File supports max n: %d\n", n, f_n);
        exit(3);
    }
}

// Makes indexing the 2D array easier
int get(int i, int j, int n) { return i * n + j; }

void print_data(int rank, int world_size, int rows, int n, double *data) {
    if (rank == 0) {
        printf("\n");
    }
    for (int r = 0; r < world_size; r++) {
        MPI_Barrier(MPI_COMM_WORLD);
        if (r == rank) {
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < n; j++) {
                    printf("%f\t", data[get(i, j, n)]);
                }
                printf("\n");
            }
        }
    }
    if (rank == 0) {
        printf("\n");
    }
    MPI_Barrier(MPI_COMM_WORLD);
}

double average(int i, int j, int n, int rows, double *data, double *above,
               double *below) {
    double north, east, south, west;

    // Check to see if we need to read from 'borrowed' rows
    north = i == 0 ? above[j] : data[get(i - 1, j, n)];
    south = i == rows - 1 ? below[j] : data[get(i + 1, j, n)];

    // No need to check if outside bounds for east/west because average is
    // never called on first/last column values
    east = data[get(i, j + 1, n)];
    west = data[get(i, j - 1, n)];

    // Average the numbers and return the result
    return (north + east + south + west) / 4.0;
}

// Averages all the numbers in a row
// Takes data and copy (of data)
// Data is treated as mutable and copy is treated as immutable
// relaxation function is responsible for synchronizing the two
// Returns 1 if this row is 'done'
int compute_row(int i, int rows, int n, double *data, double *copy,
                double *above, double *below, double precision) {
    double avg;
    double diff;
    int done = 1;

    for (int j = 1; j < n - 1; j++) {
        avg = average(i, j, n, rows, copy, above, below);
        diff = fabs(copy[get(i, j, n)] - avg);
        done = done * (precision < diff ? 0 : 1);
        data[get(i, j, n)] = avg;
    }

    return done;
}

void relaxation(double *data, int rank, int world_size, int rows, int n,
                double precision) {
    // Any s_ is a send request
    // Any r_ is a recv request
    MPI_Request s_data_above, s_data_below;
    MPI_Request r_data_above, r_data_below;
    MPI_Request s_done_above, s_done_below;
    MPI_Request r_done_above, r_done_below;

    MPI_Request *send_requests;
    int send_count = 0;

    int kill_signal = 1;
    int done = 0;
    int done_above = 0;
    int done_below = 0;
    int received_final_above = 0;
    int received_final_below = 0;

    double *data_copy;
    double *temp_send_above;
    double *temp_send_below;

    double above[n];
    double below[n];

    // Malloc required space for variables
    // and check the result of the malloc
    send_requests = malloc((size_t)2 * sizeof(MPI_Request));
    if (send_requests == NULL) {
        printf("Error allocating memory for 'double *send_requests' from rank "
               "%d'\n",
               rank);
        exit(4);
    }
    data_copy = malloc((size_t)rows * (size_t)n * sizeof(double));
    if (data_copy == NULL) {
        printf(
            "Error allocating memory for 'double *data_copy' from rank %d'\n",
            rank);
        exit(4);
    }
    temp_send_above = malloc((size_t)n * sizeof(double));
    if (temp_send_above == NULL) {
        printf("Error allocating memory for 'double *temp_send_above' from "
               "rank %d'\n",
               rank);
        exit(4);
    }
    temp_send_below = malloc((size_t)n * sizeof(double));
    if (temp_send_below == NULL) {
        printf("Error allocating memory for 'double *temp_send_below' from "
               "rank %d'\n",
               rank);
        exit(4);
    }

    // Open recv handles for kill signals
    if (rank != 0) {
        MPI_Irecv(&done_above, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD,
                  &r_done_above);
    }

    if (rank != world_size - 1) {
        MPI_Irecv(&done_below, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD,
                  &r_done_below);
    }

    // Loop until done
    while (done == 0) {
        done = 1;
        send_count = 0;

        // From MPI_Isend documentation:
        // The sender should not modify any part of the send buffer after a
        // nonblocking send operation is called, until the send completes
        // So need to make copies of the data we need to send
        // No need to make copies of data that isn't going to be sent
        if (rank != 0 && rank != world_size - 1) {
            for (int j = 0; j < n; j++) {
                temp_send_above[j] = data[get(0, j, n)];
                temp_send_below[j] = data[get(rows - 1, j, n)];
            }
        } else if (rank != 0) {
            for (int j = 0; j < n; j++) {
                temp_send_above[j] = data[get(0, j, n)];
            }
        } else if (rank != world_size - 1) {
            for (int j = 0; j < n; j++) {
                temp_send_below[j] = data[get(rows - 1, j, n)];
            }
        }

        // It is necessary to have separate 'done' flags
        // One for when above is finished
        // And one for when we have received the final data from above
        // There is a subtle race condition that can occur and if we have
        // opened a recv handle, later on the code will deadlock
        // waiting for data that will never be sent (or has already been
        // received)

        // If we have NOT received the final data from above, open recv handle
        // If above is NOT done, send the data to above
        if (rank != 0) {
            if (received_final_above == 0) {
                MPI_Irecv(above, n, MPI_DOUBLE, rank - 1, MPI_TAG_UB,
                          MPI_COMM_WORLD, &r_data_above);
            }

            if (done_above == 0) {
                MPI_Isend(temp_send_above, n, MPI_DOUBLE, rank - 1, MPI_TAG_UB,
                          MPI_COMM_WORLD, &s_data_above);
                send_requests[send_count++] = s_data_above;
            }
        }

        // If we have NOT received the final data from below, open recv handle
        // If below is NOT done, send the data to below
        if (rank != world_size - 1) {
            if (received_final_below == 0) {
                MPI_Irecv(below, n, MPI_DOUBLE, rank + 1, MPI_TAG_UB,
                          MPI_COMM_WORLD, &r_data_below);
            }

            if (done_below == 0) {
                MPI_Isend(temp_send_below, n, MPI_DOUBLE, rank + 1, MPI_TAG_UB,
                          MPI_COMM_WORLD, &s_data_below);
                send_requests[send_count++] = s_data_below;
            }
        }

        // Make a copy of the data so that it computes the answer correctly
        // Could have memory optimized this and just swapped the pointers
        // rather than copy each iteration
        for (int i = 0; i < rows * n; i++) {
            data_copy[i] = data[i];
        }

        // Compute 'middle' rows first while waiting for data from above/below
        for (int i = 1; i < rows - 1; i++) {
            done *=
                compute_row(i, rows, n, data, data_copy, NULL, NULL, precision);
            ;
        }

        // Wait on data being sent from above
        // Only wait if we have NOT received the final communications from above
        // Compute top row
        if (rank != 0) {
            if (received_final_above == 0) {
                MPI_Wait(&r_data_above, MPI_STATUS_IGNORE);
            }
            done *= compute_row(0, rows, n, data, data_copy, above, below,
                                precision);
            received_final_above = done_above;
        }

        // Wait on data being sent from below
        // Only wait if we have NOT received the final communications from below
        // Compute bottom row
        if (rank != world_size - 1) {
            if (received_final_below == 0) {
                MPI_Wait(&r_data_below, MPI_STATUS_IGNORE);
            }
            done *= compute_row(rows - 1, rows, n, data, data_copy, above,
                                below, precision);
            received_final_below = done_below;
        }

        // Wait on r_sends
        // Make sure the data has been sent before starting the next iteration
        if (done == 0) {
            MPI_Waitall(send_count, send_requests, MPI_STATUS_IGNORE);
        }
    }

    // printf("---- Rank %d managed to escape. Took %d iterations ----\n", rank,
    // iteration);

    send_count = 0;

    // Send kill signal to above/below
    if (rank != 0) {
        MPI_Isend(&kill_signal, 1, MPI_INT, rank - 1, 1, MPI_COMM_WORLD,
                  &s_done_above);
        send_requests[send_count++] = s_done_above;
    }

    if (rank != world_size - 1) {
        MPI_Isend(&kill_signal, 1, MPI_INT, rank + 1, 1, MPI_COMM_WORLD,
                  &s_done_below);
        send_requests[send_count++] = s_done_below;
    }

    // Wait until kill signals have been received
    MPI_Waitall(send_count, send_requests, MPI_STATUS_IGNORE);

    send_count = 0;

    // Send final data to above/below
    if (rank != 0 && done_above == 0) {
        MPI_Isend(&data[0], n, MPI_DOUBLE, rank - 1, MPI_TAG_UB, MPI_COMM_WORLD,
                  &s_data_above);
        send_requests[send_count++] = s_data_above;
    }

    if (rank != world_size - 1 && done_below == 0) {
        MPI_Isend(&data[get(rows - 1, 0, n)], n, MPI_DOUBLE, rank + 1,
                  MPI_TAG_UB, MPI_COMM_WORLD, &s_data_below);
        send_requests[send_count++] = s_data_below;
    }

    // Don't want to free the data before it has been read out of
    MPI_Waitall(send_count, send_requests, MPI_STATUS_IGNORE);

    // Free variables
    free(send_requests);
    free(data_copy);
    free(temp_send_above);
    free(temp_send_below);

    // Don't return from the function until all processes are finished
    MPI_Barrier(MPI_COMM_WORLD);
}

int main(int argc, char **argv) {
    MPI_File fr;
    char *filename = "randNumbers";
    int err, f_start, count;
    int matrix_n = 8;
    int rows, remainder;
    int world_size, rank;
    double precision = 0.00001;
    struct timespec start, stop;

    // Read in command line arguments
    if (argc == 3) {
        char *p1, *p2;
        long conv1 = strtol(argv[1], &p1, 10);
        double conv2 = strtod(argv[2], &p2);
        if (errno != 0 || *p1 != '\0' || *p2 != '\0') {
            printf("Error converting command line arguments");
            return 1;
        }

        matrix_n = (int)conv1;
        precision = conv2;
    }

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (rank == 0) {
        get_time(&start);
    }

    if (matrix_n < world_size) {
        printf("Error, cannot have more processes than n: %d", matrix_n);
        exit(5);
    }

    // Divide up the rows equally between the processes
    // Take n = 95, world size = 16
    // 95 / 16 = 5
    // 95 % 16 = 15
    // Give each process 5 rows
    // Add 1 row to the first 15 processes
    rows = matrix_n / world_size;
    remainder = matrix_n % world_size;

    // Added this to simplify relaxation code, it would be possible to run
    // with 1 row however it complicates the code for little benefit
    if (rows == 1) {
        if (rank == 0) {
            printf("Error, n: %d too small. Cannot have processes with 1 row\n",
                   matrix_n);
            printf("This error occurred because %d/%d=%d\n", matrix_n,
                   world_size, rows);
            exit(5);
        }
    }

    // Adds 1 to the relevant processes
    // Also calculates the file offset position for each rank
    if (rank <= remainder - 1) {
        rows++;
        f_start = rank * rows * matrix_n;
    } else {
        f_start = (rows + 1) * (remainder * matrix_n) +
                  ((rank - remainder) * rows * matrix_n);
    }

    f_start++;
    count = matrix_n * rows;

    // Open the file
    open_file(&fr, rank, filename);

    // Read first value and check against matrix_n
    // Checks to make sure there are enough values in the file
    if (rank == 0) {
        check_file(fr, matrix_n);
    }

    if (rank == 0) {
        printf("%d processes solving %dx%d matrix (%d items) to %f precision. "
               "Numbers taken from %s file\n",
               world_size, matrix_n, matrix_n, matrix_n * matrix_n, precision,
               filename);
        printf("----------------------------\n\n");
    }

    // Don't want to start reading file if there aren't enough entries
    MPI_Barrier(MPI_COMM_WORLD);

    int *temp = malloc((size_t)count * sizeof(int));
    if (temp == NULL) {
        printf(
            "Error allocating memory for 'int *temp' variable from rank %d'\n",
            rank);
        exit(4);
    }

    // Read from the file at varying offsets (dependent on rank)
    // File is binary integer data so must be read into integer datatype
    err = MPI_File_read_at_all(fr, (MPI_Offset)((size_t)f_start * sizeof(int)),
                               temp, count, MPI_INT, MPI_STATUS_IGNORE);
    if (err) {
        if (rank == 0) {
            printf("Error reading file\n");
        }
        exit(2);
    }
    MPI_File_close(&fr);
    double data[count];

    // TODO not the best way?
    // Could read in doubles rather than read in ints and converting?
    for (int i = 0; i < count; i++) {
        data[i] = (double)temp[i];
    }
    free(temp);

    // print_data(rank, world_size, rows, matrix_n, data);
    // MPI_Barrier(MPI_COMM_WORLD);

    // Call relaxation method
    relaxation(data, rank, world_size, rows, matrix_n, precision);

    if (rank <= 15) {
        print_data(rank, world_size, rows, matrix_n, data);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Print times
    if (rank == 0) {
        get_time(&stop);
        printf("\n\nTotal time: %lf\n", elapsed_time(start, stop));
    }

    MPI_Finalize();
    return 0;
}
