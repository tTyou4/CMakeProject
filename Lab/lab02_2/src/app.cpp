#include <app.hpp>
using namespace std;

void getGaussSolution(const vector<vector<double>>& A, const vector<double>& b, vector<double>& x) {
    int n = A.size();
    vector<vector<double>> A_copy = A;
    vector<double> b_copy = b;
    x.assign(n, 0.0);

    for (int i = 0; i < n; ++i) {
        int pivot = i;
        for (int j = i + 1; j < n; ++j) {
            if (abs(A_copy[j][i]) > abs(A_copy[pivot][i])) pivot = j;
        }
        swap(A_copy[i], A_copy[pivot]);
        swap(b_copy[i], b_copy[pivot]);

        for (int j = i + 1; j < n; ++j) {
            double factor = A_copy[j][i] / A_copy[i][i];
            for (int k = i; k < n; ++k) A_copy[j][k] -= factor * A_copy[i][k];
            b_copy[j] -= factor * b_copy[i];
        }
    }

    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += A_copy[i][j] * x[j];
        x[i] = (b_copy[i] - sum) / A_copy[i][i];
    }
}

int main() {
    int n = 4;
    vector<vector<double>> A = {
        { 36.4223, -20.6974,  14.4920,  32.8637 },
        { 37.1915,  41.2522,  -5.9842, -14.5971 },
        { -7.2291, -49.9125,  27.3531, -41.0175 },
        { 18.3893, -34.5912, -45.1528, -48.1415 }
    };
    vector<double> b = { 14.7615, 41.4143, -6.8298, -25.7404 };

    vector<double> x_Gauss;
    getGaussSolution(A, b, x_Gauss);

    vector<vector<double>> A_eq(n, vector<double>(n, 0.0));
    vector<double> b_eq(n, 0.0);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            for (int k = 0; k < n; ++k) {
                A_eq[i][j] += A[k][i] * A[k][j];
            }
        }
        for (int k = 0; k < n; ++k) {
            b_eq[i] += A[k][i] * b[k];
        }
    }

    cout << "Matrix A has NO diagonal dominance, so convergence is not guaranteed." << endl;
    cout << "Equivalent system A^T A x = A^T b is used (it is symmetric positive definite)." << endl << endl;

    char repeat = 'y';
    while (repeat == 'y' || repeat == 'Y') {
        double eps;
        cout << "Enter accuracy eps (for example 1e-3, 1e-6): ";
        if (!(cin >> eps)) break;

        double max_row_sum = 0.0;
        for (int i = 0; i < n; ++i) {
            double sum = 0.0;
            for (int j = 0; j < n; ++j) sum += abs(A_eq[i][j]);
            max_row_sum = max(max_row_sum, sum);
        }
        double tau = 1.0 / (max_row_sum * 1.5);

        vector<double> x_simple(n, 0.0);
        int iter_simple = 0;

        while (iter_simple < 10000) {
            vector<double> x_next(n);
            for (int i = 0; i < n; ++i) {
                double Ax = 0;
                for (int j = 0; j < n; ++j) Ax += A_eq[i][j] * x_simple[j];
                x_next[i] = x_simple[i] - tau * (Ax - b_eq[i]);
            }

            double diff = 0.0;
            for (int i = 0; i < n; ++i) diff = max(diff, abs(x_next[i] - x_simple[i]));

            if (isnan(diff) || isinf(diff) || diff > 1e10) {
                tau /= 2.0;
                iter_simple = 0;
                x_simple.assign(n, 0.0);
                continue;
            }

            x_simple = x_next;
            iter_simple++;
            if (diff < eps) break;
        }

        cout << "\nSimple iteration: " << iter_simple << " iterations" << endl;
        for (int i = 0; i < n; ++i) {
            cout << "x[" << i + 1 << "] = " << setw(12) << fixed << setprecision(8) << x_simple[i] << endl;
        }

        double norm_res_s = 0.0, norm_diff_s = 0.0;
        for (int i = 0; i < n; ++i) {
            double Ax = 0;
            for (int j = 0; j < n; ++j) Ax += A[i][j] * x_simple[j];
            norm_res_s = max(norm_res_s, abs(b[i] - Ax));
            norm_diff_s = max(norm_diff_s, abs(x_simple[i] - x_Gauss[i]));
        }
        cout << "||b - A x|| = " << scientific << setprecision(3) << norm_res_s
            << ",  ||x - x_Gauss|| = " << norm_diff_s << endl << endl;

        vector<double> x_seidel(n, 0.0);
        int iter_seidel = 0;

        while (iter_seidel < 10000) {
            vector<double> x_old = x_seidel;
            for (int i = 0; i < n; ++i) {
                double sum = 0.0;
                for (int j = 0; j < n; ++j) {
                    if (j != i) sum += A_eq[i][j] * x_seidel[j];
                }
                x_seidel[i] = (b_eq[i] - sum) / A_eq[i][i];
            }

            double diff = 0.0;
            for (int i = 0; i < n; ++i) diff = max(diff, abs(x_seidel[i] - x_old[i]));

            iter_seidel++;
            if (diff < eps) break;
        }

        cout << "Gauss-Seidel: " << iter_seidel << " iterations" << endl;
        for (int i = 0; i < n; ++i) {
            cout << "x[" << i + 1 << "] = " << setw(12) << fixed << setprecision(8) << x_seidel[i] << endl;
        }

        double norm_res_z = 0.0, norm_diff_z = 0.0;
        for (int i = 0; i < n; ++i) {
            double Ax = 0;
            for (int j = 0; j < n; ++j) Ax += A[i][j] * x_seidel[j];
            norm_res_z = max(norm_res_z, abs(b[i] - Ax));
            norm_diff_z = max(norm_diff_z, abs(x_seidel[i] - x_Gauss[i]));
        }
        cout << "||b - A x|| = " << scientific << setprecision(3) << norm_res_z
            << ",  ||x - x_Gauss|| = " << norm_diff_z << endl << endl;

        cout << "Repeat? (y/n): ";
        cin >> repeat;
        cout << endl;
    }

    return 0;
}