#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <errno.h>

double average(double n, double e, double s, double w);
void print_arr(double **arr, int n);
void copy_arr(double **source, double **target, int n);

// Prints the given array
void print_arr(double **arr, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%f\t", arr[i][j]);
        }
        printf("\n");
    }
}

// Performs a deep copy on the given arrays
void copy_arr(double **source, double **target, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            target[i][j] = source[i][j];
        }
    }
}

double average(double n, double e, double s, double w) {
    return (n + e + s + w) / 4.0;
}

int main(int argc, char *argv[]) {
    struct timespec start, stop;

    // Start timer 
    if (clock_gettime(CLOCK_REALTIME, &start) == -1) {
        printf("Error getting clock time");
        return 1;
    }

    int n = 10;
    int file_n, a;
    double precision = 0.00001;

    // Read in command line arguments 
    if (argc == 3) {
        char *p1, *p2;
        long conv1 = strtol(argv[1], &p1, 10);
        double conv2 = strtod(argv[2], &p2);
        if (errno != 0 || *p1 != '\0') {
            printf("Error converting command line arguments");
            return 1;
        }

        n = (int)conv1;
        precision = conv2;
    }

    printf(
        "Running program with %dx%d matrix and %d thread to %f precision\n\n",
        n, n, 1, precision);
    
    // Create 2D arrays
    double **arr1 = malloc((size_t)n * sizeof(double *));
    if (arr1 == NULL) {
        printf("Malloc failed");
        return 1;
    }
    double *buf = malloc((size_t)n * (size_t)n * sizeof(double));
    if (buf == NULL) {
        printf("Malloc failed");
        return 1;
    }
    double **arr2 = malloc((size_t)n * sizeof(double *));
    if (arr2 == NULL) {
        printf("Malloc failed");
        return 1;
    }
    double *buf2 = malloc((size_t)n * (size_t)n * sizeof(double));
    if (buf2 == NULL) {
        printf("Malloc failed");
        return 1;
    }

    for (int i = 0; i < n; i++) {
        arr1[i] = buf + i * n;
        arr2[i] = buf2 + i * n;
    }

    // Read in from file
    FILE *fr = fopen("randNumbers", "r");
    if (fr == NULL) {
        printf("fopen failed");
        return 1;
    }
    if (fread(&file_n, sizeof(int), 1, fr) == 0) {
        printf("Error reading from file\n");
        return 1;
    }
    if (n > file_n) {
        printf("Provided n: %d. ", n);
        printf("Matrix %dx%d is larger than file size (%d)", n, n, file_n);
        return 1;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            if (fread(&a, sizeof(int), 1, fr) == 0) {
                printf("Error reading from file");
                return 1;
            }
            arr1[i][j] = (double)a;
            arr2[i][j] = (double)a;
        }
    }
    fclose(fr);

    // Loop over computing values
    double **temp;
    int b_break = 0;
    while (b_break == 0) {
        b_break = 1;
        for (int i = 1; i < n - 1; i++) {
            for (int j = 1; j < n - 1; j++) {
                double n = arr1[i - 1][j];
                double e = arr1[i][j + 1];
                double s = arr1[i + 1][j];
                double w = arr1[i][j - 1];
                double result = average(n, e, s, w);

                arr2[i][j] = result;

                // Check if the difference is more or less than precision 
                if (fabs(result - arr1[i][j]) > precision) {
                    b_break = 0;
                }
            }
        }
        temp = arr1;
        arr1 = arr2;
        arr2 = temp;
    }

    print_arr(arr1, n);

    // Tidy up 
    free(arr1);
    free(buf);
    free(arr2);
    free(buf2);

    // Compute total runtime 
    if (clock_gettime(CLOCK_REALTIME, &stop) == -1) {
        printf("Error getting clock time");
        return 1;
    }

    double total_time =
        (double)(stop.tv_sec - start.tv_sec) +
        (double)(stop.tv_nsec - start.tv_nsec) / (double)1000000000L;

    printf("\n\nTotal time: %lf\n", total_time);
}