#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include "omp.h"

bool solutionFound = false;
double start, end;
int THRESHOLD;

//DON'T CHANGE THIS FUNCTION
void
get_constraints(std::vector<std::pair<std::pair<int, int>, std::pair<int, int> > > &constraints, std::ifstream &file) {
    int ctr = 1;
    std::string constraint;
    while (!file.eof()) {
        std::getline(file, constraint);
        if (constraint != "") {
            std::cout << "Constraint " << ctr++ << ": " << constraint << std::endl;
            std::stringstream ss(constraint);
            int val1, val2, val3, val4; // Coordinate (val1, val2) > (val3, val4)
            ss >> val1 >> val2 >> val3 >> val4;
            constraints.push_back(
                    std::pair<std::pair<int, int>, std::pair<int,
                            int> >(std::pair<int, int>(val1 - 1, val2 - 1),
                                   std::pair<int, int>(val3 - 1, val4 - 1)));
        }
    }
}


//DON'T CHANGE THIS FUNCTION
void read_matrix(int **&matrix, std::ifstream &file, int & size) {
    matrix = new int *[size];

    for (int i = 0; i < size; i++) {
        matrix[i] = new int[size];
    }

    int val;
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            file >> val;
            matrix[i][j] = val;
        }
    }
}

//USE THIS FUNCTION WHILE CHECKING FINAL RESULT
//YOU MAY FURTHER OPTIMIZE THIS FUNCTION TO USE FOR CALCULATION
bool solved(int **&matrix, std::vector<std::pair<std::pair<int, int>, std::pair<int, int> > > & constraints, int & size) {
    for (int i = 0; i < constraints.size(); i++) {
        if (matrix[constraints[i].first.first][constraints[i].first.second] <
            matrix[constraints[i].second.first][constraints[i].second.second])
            return false;
    }

    std::vector<int> rows;
    std::vector<int> cols;
    for (int rc = 0; rc < size; rc++) {
        for (int s = 0; s < size; s++) {
            rows.push_back(matrix[rc][s]);
            cols.push_back(matrix[s][rc]);
        }

        std::sort(rows.begin(), rows.end());
        std::sort(cols.begin(), cols.end());

        if ((rows[0] == -1) || (cols[0] == -1))
            return false;

        for (int i = 0; i < size - 1; i++) {
            if ((rows[i] == rows[i + 1]) || (cols[i] == cols[i + 1])) {
                return false;
            }
        }

        rows.clear();

        cols.clear();
    }

    return true;

}

bool canUseInRow(int **&matrix, int & currentRowIndex, int & value, int & size) {
    int *currentRow = matrix[currentRowIndex];
    for (int i = 0; i < size; i++) {
        if (currentRow[i] == value)
            return false;
    }

    return true;
}

bool canUseInColumn(int **&matrix, int & currentColIndex, int & value, int & size) {
    for (int i = 0; i < size; i++) {
        if (matrix[i][currentColIndex] == value)
            return false;
    }

    return true;
}

bool isSatisfyConstraints(int **&matrix, int & currentRow, int & currentCol, int & value,
                          std::vector<std::pair<std::pair<int, int>, std::pair<int, int> >

                          > &constraints) {

    int numConstraints = constraints.size();
    for (
            int i = 0;
            i < numConstraints;
            i++) {
        std::pair<std::pair<int, int>, std::pair<int, int> > constraint = constraints[i];
        std::pair<int, int> start = constraint.first;
        std::pair<int, int> end = constraint.second;

        if (currentRow == start.
                first && currentCol
                         == start.second) {
            if (
                    value < matrix[end.first][end.second] && matrix[end.first][end.second]
                                                             != -1)
                return false;
        } else if (currentRow == end.
                first && currentCol
                         == end.second) {
            if (value > matrix[start.first][start.second] && matrix[start.first][start.second] != -1)
                return false;
        }
    }

    return true;
}

bool isCurrentCellEmpty(int **&matrix, int & currentRow, int & currentCol) {
    return matrix[currentRow][currentCol] == -1;
}

bool doesSatisfyRules(int **&matrix, int & currentRow, int & currentCol, int & value,
                      std::vector<std::pair<std::pair<int, int>, std::pair<int, int> >

                      > &constraints,
                      int size
) {
    return
            canUseInRow(matrix, currentRow, value, size
            ) &&
            canUseInColumn(matrix, currentCol, value, size
            ) &&
            isSatisfyConstraints(matrix, currentRow, currentCol, value, constraints
            ) &&
            isCurrentCellEmpty(matrix, currentRow, currentCol
            );
}

std::pair<int, int> findEmptyCell(int **&matrix, int & size) {
    for (int i = 0; i < size; i++) { // row
        for (int j = 0; j < size; j++) { // column
            if (matrix[i][j] == -1)
                return std::pair<int, int>(i, j);
        }
    }

    return std::pair<int, int>(-1, -1);
}

void solve(int **&matrix, std::vector<std::pair<std::pair<int, int>, std::pair<int, int> > > &constraints, int & size) {
// Find a empty cell
    std::pair<int, int> emptyCell = findEmptyCell(matrix, size);
    int currentRow = emptyCell.first;
    int currentCol = emptyCell.second;

// If there is no empty cell, then we are done
    if (currentRow == -1 && currentCol == -1) {
        solutionFound = true;
        std::cout << "Can solve: " << solved(matrix, constraints, size) << std::endl;
        return;
    }

    for (int value = 1; value <= size && !solutionFound; value++) {
        if (doesSatisfyRules(matrix, currentRow, currentCol, value, constraints, size)) {
            matrix[currentRow][currentCol] = value;

            solve(matrix, constraints, size);

            // Invalidate failed
            if (!solutionFound)
                matrix[currentRow][currentCol] = -1;
        }
    }

    return;
}

void solveDesign2(int **&matrix, std::vector<std::pair<std::pair<int, int>, std::pair<int, int> > > &constraints, int & size) {
// Find a empty cell
    std::pair<int, int> emptyCell = findEmptyCell(matrix, size);
    int currentRow = emptyCell.first;
    int currentCol = emptyCell.second;

// If there is no empty cell, then we are done
    if (currentRow == -1 && currentCol == -1) {
        solutionFound = true;
        std::cout << "Can solve: " << solved(matrix, constraints, size) << std::endl;
        return;
    }

    for (int value = 1; value <= size && !solutionFound; value++) {
        if (doesSatisfyRules(matrix, currentRow, currentCol, value, constraints, size)) {
            matrix[currentRow][currentCol] = value;

            if (currentRow < THRESHOLD) {
                solveDesign2(matrix, constraints, size);
            } else {
                // create copy
                int **matrixNew;
                matrixNew = new int *[size];
                for (int i = 0; i < size; i++) {
                    matrixNew[i] = new int[size];
                }

                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < size; j++) {
                        matrixNew[i][j] = matrix[i][j];
                    }
                }

#pragma omp task
                solve(matrixNew, constraints, size);
            }

            // Invalidate failed
            if (!solutionFound)
                matrix[currentRow][currentCol] = -1;
        }
    }

    return;
}

int main(int argc, char **argv) {
    std::string filename(argv[1]);

    std::ifstream file;
    file.open(filename.c_str());
    int size;

    file >> size;
    std::cout << "Size: " << size << std::endl;

    THRESHOLD = size / 2;

    int **matrix;

    read_matrix(matrix, file, size);

    std::vector<std::pair<std::pair<int, int>, std::pair<int, int> > > constraints;
    get_constraints(constraints, file);

    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            std::cout << matrix[i][j] << " ";
        }
        std::cout << std::endl;
    }

    start = omp_get_wtime();
#pragma omp parallel
    {
#pragma omp single
        solveDesign2(matrix, constraints, size);
    };
    //Measure only solution function's execution time
    end = omp_get_wtime();
    std::cout << "Duration in seconds: " << end - start << std::endl;

    /*if (!solutionFound)
        std::cout << "Can solve: " << solved(matrix, constraints, size) << std::endl;*/

    //DELETE//
    for (int i = 0; i < size; i++) {
        delete matrix[i];
    }

    delete[] matrix;
    //DELETE//
}
