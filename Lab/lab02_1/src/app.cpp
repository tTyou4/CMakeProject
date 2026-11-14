#include <app.hpp>

using namespace std;

void printMatrix(const vector<vector<double>>& A) {
    int n = A.size();
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << setw(12) << fixed << setprecision(4) << A[i][j] << " ";
        }
        cout << endl;
    }
}

bool gaussSolve(vector<vector<double>> A, vector<double> b, vector<double>& x, double& det) {
    int n = A.size();
    det = 1.0;
    int swaps = 0;

    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        if (abs(A[maxRow][i]) < 1e-12) {
            det = 0.0;
            return false;
        }

        if (maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(b[i], b[maxRow]);
            swaps++;
        }

        det *= A[i][i];

        for (int k = i + 1; k < n; ++k) {
            double factor = A[k][i] / A[i][i];
            for (int j = i; j < n; ++j) {
                A[k][j] -= factor * A[i][j];
            }
            b[k] -= factor * b[i];
        }
    }

    if (swaps % 2 != 0) {
        det = -det;
    }

    x.resize(n);
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        x[i] = (b[i] - sum) / A[i][i];
    }

    return true;
}

bool inverseMatrix(vector<vector<double>> A, vector<vector<double>>& inv) {
    int n = A.size();
    inv.assign(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        inv[i][i] = 1.0;
    }

    for (int i = 0; i < n; ++i) {
        int maxRow = i;
        for (int k = i + 1; k < n; ++k) {
            if (abs(A[k][i]) > abs(A[maxRow][i])) {
                maxRow = k;
            }
        }

        if (abs(A[maxRow][i]) < 1e-12) return false;

        if (maxRow != i) {
            swap(A[i], A[maxRow]);
            swap(inv[i], inv[maxRow]);
        }

        double pivot = A[i][i];
        for (int j = 0; j < n; ++j) {
            A[i][j] /= pivot;
            inv[i][j] /= pivot;
        }

        for (int k = 0; k < n; ++k) {
            if (k != i) {
                double factor = A[k][i];
                for (int j = 0; j < n; ++j) {
                    A[k][j] -= factor * A[i][j];
                    inv[k][j] -= factor * inv[i][j];
                }
            }
        }
    }
    return true;
}

int main() {
    vector<vector<double>> A = {
        { 36.4223, -20.6974,  14.4920,   32.8637 },
        { 37.1915,  41.2522,  -5.9842,  -14.5971 },
        { -7.2291, -49.9125,  27.3531,  -41.0175 },
        { 18.3893, -34.5912, -45.1528,  -48.1415 }
    };

    vector<double> b = { 14.7615, 41.4143, -6.8298, -25.7404 };
    int n = A.size();

    cout << "Matrix A:\n";
    printMatrix(A);
    cout << endl;

    for (int i = 0; i < n; ++i) {
        cout << "b[" << i + 1 << "] = " << setw(10) << fixed << setprecision(4) << b[i] << endl;
    }
    cout << endl;

    vector<double> x;
    double det = 0.0;
    gaussSolve(A, b, x, det);

    cout << "Solution (Gauss method):\n";
    for (int i = 0; i < n; ++i) {
        cout << "x[" << i + 1 << "] =  " << fixed << setprecision(8) << x[i] << endl;
    }
    cout << endl;

    vector<double> r(n);
    double max_norm_r = 0.0;
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += A[i][j] * x[j];
        }
        r[i] = sum - b[i];
        if (abs(r[i]) > max_norm_r) {
            max_norm_r = abs(r[i]);
        }
    }

    cout << "Residual r = Ax - b:\n";
    for (int i = 0; i < n; ++i) {
        double sum = 0.0;
        for (int j = 0; j < n; ++j) {
            sum += A[i][j] * x[j];
            cout << (A[i][j] >= 0 && j > 0 ? "+ " : "") << A[i][j] << "*x" << j + 1 << " ";
        }
        r[i] = sum - b[i];
        cout << "- (" << b[i] << ") = " << scientific << setprecision(4) << r[i] << endl;
        if (abs(r[i]) > max_norm_r) {
            max_norm_r = abs(r[i]);
        }
    }
    cout << "r (max norm) = " << scientific << setprecision(4) << max_norm_r << endl;
    cout << endl;

    cout << "Determinant det(A) = " << fixed << setprecision(6) << det << endl;
    cout << endl;

    vector<vector<double>> invA;
    inverseMatrix(A, invA);

    cout << "Inverse matrix A^-1:\n";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cout << setw(15) << fixed << setprecision(8) << invA[i][j] << " ";
        }
        cout << endl;
    }
    cout << endl;

    double max_err = 0.0;
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double elem = 0.0;
            for (int k = 0; k < n; ++k) {
                elem += A[i][k] * invA[k][j];
            }
            double target = (i == j) ? 1.0 : 0.0;
            if (abs(elem - target) > max_err) {
                max_err = abs(elem - target);
            }
        }
    }

    cout << "Check: max|A*A^-1 - E| = " << scientific << setprecision(4) << max_err << endl;

    return 0;
}