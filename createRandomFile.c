#include <stdio.h>
#include <stdlib.h>
#include <time.h>

int main() {
    // Creates a file with randomly generated integers up to n 
    int n = 10000;
    srand(time(0));

    FILE *fw = fopen("randNumbers", "wb");
    if (fw == NULL) {
        return 1;
    }

    if (fwrite(&n, sizeof(int), 1, fw) == 0) {
        return 1;
    }

    for (int i = 0; i < n * n; i++) {
        int number = rand() % 1000;
        if (fwrite(&number, sizeof(int), 1, fw) == 0) {
            return 1;
        }
    }

    fclose(fw);
    return 0;
}
