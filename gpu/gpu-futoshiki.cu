#include<iostream>
#include<string>
#include<sstream>
#include<fstream>

#define SIZE 5 //Matrix size
#define INPUTSIZE 2306451

#define gpuErrchk(ans) { gpuAssert((ans), __FILE__, __LINE__); }

// Ege Alpay 19551

inline void gpuAssert(cudaError_t code, const char *file, int line, bool abort = true) {
    if (code != cudaSuccess) {
        fprintf(stderr, "GPUassert: %s %s %d\n", cudaGetErrorString(code), file, line);
        if (abort) exit(code);
    }
}

__device__ bool canUseInRow(int *matrix, int currentRowIndex, int value) {
    int startIndex = currentRowIndex * 5;
    for (int i = startIndex; i < startIndex + 5; i++) {
        if (matrix[i] == value)
            return false;
    }

    return true;
}

__device__ bool canUseInColumn(int *matrix, int currentColIndex, int value) {
    for (int i = 0; i < 5; i++) {
        if (matrix[currentColIndex + 5 * i] == value)
            return false;
    }

    return true;
}

__device__ bool
isSatisfyConstraints(int *matrix, int currentRow, int currentCol, int value, int *constraints, int numConstraints) {
    for (int i = 0; i < numConstraints; i++) {
        // -1 Since constraint index start from 1, not 0
        int startRow = constraints[4 * i] - 1;
        int startCol = constraints[4 * i + 1] - 1;
        int endRow = constraints[4 * i + 2] - 1;
        int endCol = constraints[4 * i + 3] - 1;


        if (currentRow == startRow && currentCol == startCol) {
            if (value < matrix[5 * endRow + endCol] && matrix[5 * endRow + endCol] != -1)
                return false;
        } else if (currentRow == endRow && currentCol == endCol) {
            if (value > matrix[5 * startRow + startCol] && matrix[5 * startRow + startCol] != -1)
                return false;
        }
    }

    return true;
}

__device__ bool isCurrentCellEmpty(int *matrix, int currentRow, int currentCol) {
    return matrix[5 * currentRow + currentCol] == -1;
}

__device__ bool
doesSatisfyRules(int *matrix, int currentRow, int currentCol, int value, int *constraints, int numConstraints) {
    return canUseInRow(matrix, currentRow, value) &&
           canUseInColumn(matrix, currentCol, value) &&
           isSatisfyConstraints(matrix, currentRow, currentCol, value, constraints, numConstraints) &&
           isCurrentCellEmpty(matrix, currentRow, currentCol);
}

__device__ bool canUseInRowCheck(int *matrix, int currentRowIndex, int currentColIndex, int value) {
    int startIndex = currentRowIndex * 5;
    for (int i = startIndex; i < startIndex + 5; i++) {
        if (5 * currentRowIndex + currentColIndex != i && matrix[i] == value)
            return false;
    }

    return true;
}

__device__ bool canUseInColumnCheck(int *matrix, int currentRowIndex, int currentColIndex, int value) {
    for (int row = 0; row < 5; row++) {
        if (currentRowIndex != row && matrix[currentColIndex + 5 * row] == value)
            return false;
    }

    return true;
}

__device__ bool isSolutionCorrect(int *matrix, int *constraints, int numConstraints) {
    for (int row = 0; row < 5; row++) {
        for (int col = 0; col < 5; col++) {
            int value = matrix[5 * row + col];
            // printf("Checking Row: %d, Col: %d, Value: %d \n", row, col, value);

            if (value < 1 || value > 5) {
                // printf("Invalid value: %d \n", value);
                return false;
            }

            if (!canUseInRowCheck(matrix, row, col, value)) {
                // printf("Already in use in ROW: %d \n", value);
                return false;
            }

            if (!canUseInColumnCheck(matrix, row, col, value)) {
                // printf("Already in use in COL: %d \n", value);
                return false;
            }

            if (!isSatisfyConstraints(matrix, row, col, value, constraints, numConstraints)) {
                // printf("Constraints Fails: %d \n", value);
                return false;
            }
        }
    }

    return true;
}

__device__ int *findEmptyCell(int *matrix) {
    int emptyCell[2]; // row, col

    for (int row = 0; row < 5; row++) {
        for (int col = 0; col < 5; col++) {
            if (matrix[row * 5 + col] == -1) {
                emptyCell[0] = row;
                emptyCell[1] = col;
                return emptyCell;
            }
        }
    }

    emptyCell[0] = -1;
    emptyCell[1] = -1;
    return emptyCell;
}

// Stack like implementation was used since recursion on GPU is not recommended
__device__ void solveSingleThread(int *matrix, int *constraints, int numConstraints) {
    bool solved = false;
    bool unroll = false;
    int stack[50];
    int stackIndex = 0;
    int values[25]; // 5*row + col

    // Every cells last taken value will be stored in this array
    for (int i = 0; i < 25; i++)
        values[i] = 1;

    while (!solved) {
        while (unroll) {
            int col = stack[--stackIndex];
            int row = stack[--stackIndex];

            matrix[row * 5 + col] = -1;

            if (values[5 * row + col] != 5) {
                unroll = false;
                values[5 * row + col] += 1;
            } else {
                values[5 * row + col] = 1;
            }
        }

        int *emptyCell = findEmptyCell(matrix);
        int row = emptyCell[0];
        int col = emptyCell[1];

        if (row == -1 && col == -1) {
            solved = true;
            break;
        }

        for (int value = values[5 * row + col]; value <= 5; value++) {
            if (doesSatisfyRules(matrix, row, col, value, constraints, numConstraints)) {
                matrix[row * 5 + col] = value;

                stack[stackIndex++] = row;
                stack[stackIndex++] = col;

                values[5 * row + col] = value;

                break;
            }

            if (value == 5)
                unroll = true;

            if (value == 5 && matrix[row * 5 + col] == -1) {
                values[5 * row + col] = 1;
            }
        }
    }
}

__global__ void solve(int *grids, int *constraints, int *constraint_sizes, int gridCount) {
    int futoshiki[25]; // Each thread will have a local copy
    int constraintSizeForPuzzle;

    int globalId = blockIdx.x * blockDim.x + threadIdx.x;

    if (globalId < gridCount) {
        constraintSizeForPuzzle = constraint_sizes[globalId];

        int localConstraints[60];

        // Create local copy for puzzle (OK)
        for (int i = 0; i < 25; i++) {
            futoshiki[i] = grids[globalId * 25 + i];
        }

        int constraintStartIndex = 0;
        for (int i = 0; i < globalId; i++) {
            constraintStartIndex += constraint_sizes[i] * 4;
        }

        // Create local copy for constraints (OK)
        for (int i = 0; i < constraintSizeForPuzzle; i++) {
            localConstraints[4 * i] = constraints[constraintStartIndex + 4 * i];
            localConstraints[4 * i + 1] = constraints[constraintStartIndex + 4 * i + 1];
            localConstraints[4 * i + 2] = constraints[constraintStartIndex + 4 * i + 2];
            localConstraints[4 * i + 3] = constraints[constraintStartIndex + 4 * i + 3];
        }

        solveSingleThread(futoshiki, localConstraints, constraint_sizes[globalId]);

        /*if (!isSolutionCorrect(futoshiki, localConstraints, constraint_sizes[globalId]))
            printf("Solution is not correct! Thread ID: %d \n", globalId);*/

        // Write back to solutions vector
        for (int i = 0; i < 25; i++) {
            grids[globalId * 25 + i] = futoshiki[i];
        }

        /* if (globalId == 0) {
            printf("********************* Solution ********************* \n");
            for (int i = 0; i < 25; i += 5) {
                printf(" %d ", grids[25 * globalId + i]);
                printf(" %d ", grids[25 * globalId + i + 1]);
                printf(" %d ", grids[25 * globalId + i + 2]);
                printf(" %d ", grids[25 * globalId + i + 3]);
                printf(" %d ", grids[25 * globalId + i + 4]);
                printf("\n");
            }
            printf("********************* Solution ********************* \n");
        } */
    }
}



//You can change any part in order to optimize.
//Just don't forget to measure times

int main(int argc, char **argv) {

    std::string filename(argv[1]);
    std::ifstream file(filename.c_str());
    std::ifstream scout(filename.c_str());

    int no_grids;
    file >> no_grids;

    int dummy;
    scout >> dummy;

    int ***grids = new int **[no_grids];
    int **constraints = new int *[no_grids];

    for (int i = 0; i < no_grids; i++) {
        grids[i] = new int *[SIZE];
        for (int j = 0; j < SIZE; j++) {
            grids[i][j] = new int[SIZE];
        }
    }

    int elem0, elem1, elem2, elem3, elem4;
    int pre_cursor = 0;
    int cursor = 0;
    int csize = 0;

    std::string file_line;
    std::string scout_line;

    int *constraint_sizes = new int[no_grids];

    std::getline(scout, scout_line);//These are for spare lines
    std::getline(scout, scout_line);
    for (int i = 0; i < INPUTSIZE; i++) {
        std::getline(scout, scout_line);
        if (scout_line == "-------") {
            csize = i - pre_cursor - 5;
            constraint_sizes[cursor] = csize;
            cursor++;
            pre_cursor = i + 1;
        }
    }

    for (int i = 0; i < no_grids; i++) {
        constraints[i] = new int[constraint_sizes[i] * 4];
    }

    int NUM_CONSTRAINTS_ELEMENTS = 0;
    std::getline(file, file_line);
    for (int i = 0; i < no_grids; i++) {
        std::getline(file, file_line);
        for (int j = 0; j < SIZE; j++) {
            std::getline(file, file_line);
            //std::cout << "i: " << i << " file_line: " << file_line << std::endl;
            std::istringstream iss(file_line);
            iss >> elem0 >> elem1 >> elem2 >> elem3 >> elem4;
            grids[i][j][0] = elem0;
            grids[i][j][1] = elem1;
            grids[i][j][2] = elem2;
            grids[i][j][3] = elem3;
            grids[i][j][4] = elem4;
        }
        for (int c = 0; c < constraint_sizes[i]; c++) {
            std::getline(file, file_line);
            //std::cout << "i: " << i << "c line: " << file_line << std::endl;
            std::istringstream iss(file_line);
            iss >> elem0 >> elem1 >> elem2 >> elem3;
            constraints[i][4 * c] = elem0;
            constraints[i][4 * c + 1] = elem1;
            constraints[i][4 * c + 2] = elem2;
            constraints[i][4 * c + 3] = elem3;

            NUM_CONSTRAINTS_ELEMENTS += 4;
        }
    }

    // Flatten 3D Grids
    int *grids_1d = new int[no_grids * SIZE * SIZE];
    for (int numGrid = 0; numGrid < no_grids; numGrid++) {
        for (int row = 0; row < SIZE; row++) {
            for (int col = 0; col < SIZE; col++) {
                grids_1d[SIZE * SIZE * numGrid + SIZE * row + col] = grids[numGrid][row][col];
            }
        }
    }

    // Flatten 2D Constraints
    int counter = 0;
    int *constraints_1d = new int[NUM_CONSTRAINTS_ELEMENTS];
    for (int numGrid = 0; numGrid < no_grids; numGrid++) {
        int totalIterationCount = constraint_sizes[numGrid] * 4;
        for (int c = 0; c < totalIterationCount; c++) {
            constraints_1d[counter] = constraints[numGrid][c];
            counter++;
        }
    }

    /* You can access input and constraints like this
    for(int in = 0; in < 25; in++){
      std::cout << "in: " << in << std::endl;
      for(int i = 0; i < SIZE; i++){
        for(int j = 0; j < SIZE; j++){
      std::cout << grids[in][i][j] << " ";
        }
        std::cout << std::endl;
      }

      for(int i = 0; i < constraint_sizes[in]; i++){
        std::cout << constraints[in][4*i] << " " << constraints[in][4*i+1] << " " << constraints[in][4*i+2] << " " << constraints[in][4*i+3] << std::endl;
      }

    }
    */

    float time;
    cudaEvent_t start, stop;

    cudaEventCreate(&start);
    cudaEventCreate(&stop);

    //YOUR MEMORY OPERATIONS//Time accordingly
    int *d_grids;
    int *d_constraints;
    int *d_constraint_sizes;

    cudaEventRecord(start, 0);
    std::cout << "GPU Memory Allocation" << std::endl;
    cudaMalloc((void **) &d_grids, sizeof(int) * no_grids * SIZE * SIZE);
    cudaMalloc((void **) &d_constraints, sizeof(int) * NUM_CONSTRAINTS_ELEMENTS);
    cudaMalloc((void **) &d_constraint_sizes, sizeof(int) * no_grids);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("GPU Memory allocation duration: %f ms \n", time);


    cudaEventRecord(start, 0);
    std::cout << "CPU to GPU" << std::endl;
    cudaMemcpy(d_grids, grids_1d, sizeof(int) * no_grids * SIZE * SIZE, cudaMemcpyHostToDevice);
    cudaMemcpy(d_constraints, constraints_1d, sizeof(int) * NUM_CONSTRAINTS_ELEMENTS, cudaMemcpyHostToDevice);
    cudaMemcpy(d_constraint_sizes, constraint_sizes, sizeof(int) * no_grids, cudaMemcpyHostToDevice);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("CPU to GPU Data Transfer Duration: %f ms \n", time);

//YOUR MEMORY OPERATIONS//


//KERNEL CALL//Time accordingly

    int threadPerBlock = 64;
    int blockCount = (no_grids / threadPerBlock) + 1;

//KERNEL CALL//
    cudaEventRecord(start, 0);
    solve << < blockCount, threadPerBlock >> > (d_grids, d_constraints, d_constraint_sizes, no_grids);
    gpuErrchk(cudaPeekAtLastError());
    gpuErrchk(cudaDeviceSynchronize());
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("Kernel Duration: %f ms \n", time);

    //YOUR MEMORY OPERARIONS//Time accordingly
    cudaEventRecord(start, 0);
    std::cout << "GPU to CPU" << std::endl;
    cudaMemcpy(grids_1d, d_grids, sizeof(int) * no_grids * SIZE * SIZE, cudaMemcpyDeviceToHost);
    cudaDeviceSynchronize();
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("GPU to CPU Data Transfer Duration: %f ms \n", time);

    /*int gridNum = 122310;
    std::cout << "THESE ARE RESULTS: *********************" << std::endl;
    for (int i = 25 * gridNum; i < 25 * gridNum + 25; i += 5) {
        std::cout << grids_1d[i] << " ";
        std::cout << grids_1d[i + 1] << " ";
        std::cout << grids_1d[i + 2] << " ";
        std::cout << grids_1d[i + 3] << " ";
        std::cout << grids_1d[i + 4] << " ";
        std::cout << std::endl;
    }
    std::cout << "RESULTS ENDED *********************" << std::endl;*/


//YOUR MEMORY OPERARIONS//
    cudaEventRecord(start, 0);
    std::cout << "Memory Deallocation" << std::endl;
    cudaFree(d_grids);
    cudaFree(d_constraints);
    cudaFree(d_constraint_sizes);
    cudaEventRecord(stop, 0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&time, start, stop);
    printf("Memory Deallocation Duration: %f ms \n", time);


    // Create and open a text file
    std::ofstream MyFile("_solution.txt");

    // Write to file
    MyFile << no_grids << std::endl;
    for (int i = 0; i < no_grids; i++) {
        MyFile << "-------" << std::endl;

        for (int j = 0; j < 25; j += 5) {
            MyFile << grids_1d[25 * i + j];
            MyFile << " " << grids_1d[25 * i + j + 1];
            MyFile << " " << grids_1d[25 * i + j + 2];
            MyFile << " " << grids_1d[25 * i + j + 3];
            MyFile << " " << grids_1d[25 * i + j + 4];
            MyFile << std::endl;
        }
    }
    MyFile << "-------" << std::endl;

    // Close the file
    MyFile.close();


    //Deallocate
    for (int i = 0; i < no_grids; i++) {
        for (int j = 0; j < SIZE; j++) {
            delete[] grids[i][j];
        }
        delete[] grids[i];
    }
    delete[] grids;

    for (int i = 0; i < no_grids; i++) {
        delete[] constraints[i];
    }
    delete[] constraints;

    delete[] constraint_sizes;

    delete[] grids_1d;
}