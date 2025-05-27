#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

char arr[20];
int num = 0;
int j;

using Matrix = vector<vector<double>>;

void printMatrix(const Matrix& table, const vector<double>& ratios = {}) {
    cout << "\n print  this matrix \n";
    for (size_t i = 0; i < table.size(); i++) {
        for (double val : table[i]) {
            cout << setw(10) << val << " ";
        }
        cout << "\n";
    }
    cout << "----------------------------------\n";
}

void printtranspose(const Matrix& table, const vector<double>& ratios = {}) {
    cout << "\n print transpose of this matrix \n";
    for (size_t i = 0; i < table.size(); i++) {
        for (size_t val = 0; val < table.size(); val++) {
            cout << setw(10) << table[val][i] << " ";
        }
        cout << "\n";
    }
    cout << "----------------------------------\n";
}

void printTable(const Matrix& table, const vector<double>& ratios = {}) {
    int p = 1;
    for (int i = 1; i < (num); i++) {
        int m = (i >= j) ? p++ : i;
        cout << setw(10) << arr[i] << m;
    }
    cout << fixed << setprecision(2);

    for (size_t i = 0; i < table.size(); i++) {
        for (double val : table[i]) {
            cout << setw(10) << val << " ";
        }
        if (i < ratios.size() && i != table.size() - 1) {
            cout << " | Ratio: " << setw(10) << ratios[i];
        }
        cout << "\n";
    }
    cout << "----------------------------------\n";
}

Matrix transposeMatrix(const Matrix& original) {
    int rows = original.size();
    int cols = original[0].size();
    Matrix transposed(cols, vector<double>(rows));

    for (int i = 0; i < rows; ++i)
        for (int j = 0; j < cols; ++j)
            transposed[j][i] = original[i][j];

    return transposed;
}

Matrix reconstructFullTable(const Matrix& transposedRaw, int Vars, int Constraints) {
    int rows = transposedRaw[0].size();
    int cols = Vars + Constraints + 1;
    Matrix fullTable(rows, vector<double>(cols, 0.0));

    for (int i = 0; i < Constraints; ++i) {
        for (int j = 0; j < Vars; ++j) {
            fullTable[i][j] = transposedRaw[i][j];
        }
        fullTable[i][Vars + i] = 1;
        fullTable[i][cols - 1] = transposedRaw[i][Vars];
    }

    for (int j = 0; j < Vars; ++j)
        fullTable[rows - 1][j] = -transposedRaw[rows - 1][j];

    return fullTable;
}

vector<double> calculateRatios(const Matrix& table, int pivotCol) {
    vector<double> ratios;
    for (size_t i = 0; i < table.size() - 1; i++) {
        if (table[i][pivotCol] > 0) {
            double ratio = table[i].back() / table[i][pivotCol];
            ratios.push_back(ratio);
        }
        else {
            ratios.push_back(-1.0);
        }
    }
    return ratios;
}

int findPivotColumn(const Matrix& table) {
    int lastRow = table.size() - 1;
    int pivotCol = -1;
    double minValue = 0;
    for (int j = 0; j < table[0].size() - 1; j++) {
        if (table[lastRow][j] < minValue) {
            minValue = table[lastRow][j];
            pivotCol = j;
        }
    }
    return pivotCol;
}

int findPivotRow(const Matrix& table, int pivotCol, vector<double>& ratios) {
    int pivotRow = -1;
    double minRatio = 1e9;
    ratios = calculateRatios(table, pivotCol);

    for (size_t i = 0; i < ratios.size(); i++) {
        if (ratios[i] > 0 && ratios[i] < minRatio) {
            minRatio = ratios[i];
            pivotRow = i;
        }
    }
    return pivotRow;
}

void doPivot(Matrix& table, int pivotRow, int pivotCol) {
    double pivotElement = table[pivotRow][pivotCol];
    int cols = table[0].size();

    for (int j = 0; j < cols; ++j)
        table[pivotRow][j] /= pivotElement;

    for (int i = 0; i < table.size(); ++i) {
        if (i != pivotRow) {
            double factor = table[i][pivotCol];
            for (int j = 0; j < cols; ++j) {
                table[i][j] -= factor * table[pivotRow][j];
            }
        }
    }
}

void checkMultipleSolutions(const Matrix& table, int Vars) {
    vector<int> candidateCols;
    for (int j = 0; j < Vars; ++j) {
        if (table.back()[j] == 0) {
            bool isBasic = false;
            for (int i = 0; i < table.size() - 1; ++i) {
                if (table[i][j] == 1) {
                    bool allZero = true;
                    for (int k = 0; k < table.size() - 1; ++k) {
                        if (k != i && table[k][j] != 0) {
                            allZero = false;
                            break;
                        }
                    }
                    if (allZero) {
                        isBasic = true;
                        break;
                    }
                }
            }
            if (!isBasic) {
                candidateCols.push_back(j);
            }
        }
    }
    if (!candidateCols.empty()) {
        cout << "\n Multiple optimal solutions may exist!\n";
        cout << "You can pivot on:";
        for (int col : candidateCols)
            cout << " x" << col + 1;
        cout << "\n";
    }
}

void printBasicVariables(const Matrix& table, int Vars) {
    int numRows = table.size();
    cout << "\nValues of basic variables:\n";
    for (int j = 0; j < Vars; ++j) {
        bool isBasic = false;
        double value = 0.0;

        for (int i = 0; i < numRows - 1; ++i) {
            if (table[i][j] == 1) {
                bool allZero = true;
                for (int k = 0; k < numRows - 1; ++k) {
                    if (k != i && table[k][j] != 0) {
                        allZero = false;
                        break;
                    }
                }
                if (allZero) {
                    isBasic = true;
                    value = table[i].back();
                    break;
                }
            }
        }

        cout << "x" << j + 1 << " = " << (isBasic ? value : 0.0) << "\n";
    }
}

void simplex(Matrix& table, int Vars) {
    int numRows = table.size();
    int numCols = table[0].size();
    bool degeneracyDetected = false;

    while (true) {
        int pivotCol = findPivotColumn(table);
        if (pivotCol == -1) {
            cout << "\n Optimal solution found.\n";
            cout << "Minimum Z = " << table[numRows - 1][numCols - 1] << "\n";
            printBasicVariables(table, Vars);
            checkMultipleSolutions(table, Vars);
            return;
        }

        vector<double> ratios;
        int pivotRow = findPivotRow(table, pivotCol, ratios);
        if (pivotRow == -1) {
            cout << "\n Error! Solution is unbounded.\n";
            return;
        }

        double minRatio = ratios[pivotRow];
        int count = 0;
        for (double r : ratios) {
            if (r == minRatio)
                count++;
        }

        if (count > 1 && !degeneracyDetected) {
            degeneracyDetected = true;
            cout << "\n Degeneracy detected! Equal ratios found.\n";
        }

        cout << "\n Table before pivot (Row " << pivotRow + 1 << ", Col " << pivotCol + 1 << "):\n";
        printTable(table, ratios);
        doPivot(table, pivotRow, pivotCol);
        cout << "\n Table after pivot (Row " << pivotRow + 1 << ", Col " << pivotCol + 1 << "):\n";
        printTable(table);

        if (degeneracyDetected) {
            cout << "\n degeneracy solutions  exist!\n";
            cout << "Current Z = " << table[numRows - 1][numCols - 1] << "\n";
            printBasicVariables(table, Vars);
            return;
        }
    }
}

int main() {
    int Vars, Constraints;
    cout << "Enter number of variables: ";
    cin >> Vars;
    cout << "Enter number of constraints: ";
    cin >> Constraints;

    if (Vars != Constraints) {
        cout << "You can't handle this problem by using Dual method.\n";
        return 0;
    }

    int ROWS = Constraints + 1;
    int COLS = Vars + Constraints + 1;

    Matrix rawTable(ROWS, vector<double>(Vars + 1, 0.0));
    cout << "Enter coefficients for each constraint (including RHS):\n";
    for (int i = 0; i < Constraints; i++) {
        cout << "Constraint " << i + 1 << ":\n";
        for (int j = 0; j < Vars; j++) {
            cout << "x" << j + 1 << ": ";
            cin >> rawTable[i][j];
        }
        cout << "RHS: ";
        cin >> rawTable[i][Vars];
    }

    cout << "Enter coefficients of objective function (Minimize Z):\n";
    for (int j = 0; j < Vars; j++) {
        cout << "x" << j + 1 << ": ";
        cin >> rawTable[Constraints][j];
    }

    Matrix transposedRaw = transposeMatrix(rawTable);
    Matrix fullTable = reconstructFullTable(transposedRaw, Vars, Constraints);

    cout << "\n Initial Simplex Table:\n";
    printTable(fullTable);
    simplex(fullTable, Vars);

    return 0;
}
