#include <app.hpp>
using namespace std;

int main() {
    int n = 5;

    vector<double> a = { 0.000,  45.171, -36.835, -48.532, -33.107 };
    vector<double> b = { -32.341, -47.216,  36.431,  20.790, -15.902 };
    vector<double> c = { -36.836, -44.407,  22.164, -28.244,   0.000 };
    vector<double> d = { 30.239,   2.654,  29.818, -35.517,  -9.759 };

    cout << "System (a, b, c | d):" << endl;
    for (int i = 0; i < n; ++i) {
        cout << setw(10) << fixed << setprecision(3) << a[i] << " "
            << setw(10) << fixed << setprecision(3) << b[i] << " "
            << setw(10) << fixed << setprecision(3) << c[i] << "  | "
            << setw(10) << fixed << setprecision(3) << d[i] << endl;
    }
    cout << endl;

    vector<double> alpha(n, 0.0);
    vector<double> beta(n, 0.0);

    alpha[0] = -c[0] / b[0];
    beta[0] = d[0] / b[0];

    for (int i = 1; i < n; ++i) {
        double denom = b[i] + a[i] * alpha[i - 1];
        if (i < n - 1) {
            alpha[i] = -c[i] / denom;
        }
        beta[i] = (d[i] - a[i] * beta[i - 1]) / denom;
    }

    vector<double> x(n, 0.0);
    x[n - 1] = beta[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        x[i] = alpha[i] * x[i + 1] + beta[i];
    }

    cout << "Solution (Thomas algorithm):" << endl;
    for (int i = 0; i < n; ++i) {
        cout << "x[" << i + 1 << "] = " << setw(13) << fixed << setprecision(8) << x[i] << endl;
    }
    cout << endl;

    double max_res = 0.0;
    for (int i = 0; i < n; ++i) {
        double Tx = b[i] * x[i];
        if (i > 0) Tx += a[i] * x[i - 1];
        if (i < n - 1) Tx += c[i] * x[i + 1];
        max_res = max(max_res, abs(d[i] - Tx));
    }
    cout << "Residual ||d - T x|| = " << scientific << setprecision(3) << max_res << endl;

    vector<vector<double>> T(n, vector<double>(n, 0.0));
    for (int i = 0; i < n; ++i) {
        T[i][i] = b[i];
        if (i > 0) T[i][i - 1] = a[i];
        if (i < n - 1) T[i][i + 1] = c[i];
    }

    vector<double> d_copy = d;
    vector<double> x_gauss(n, 0.0);
    for (int i = 0; i < n; ++i) {
        for (int j = i + 1; j < n; ++j) {
            double factor = T[j][i] / T[i][i];
            for (int k = i; k < n; ++k) T[j][k] -= factor * T[i][k];
            d_copy[j] -= factor * d_copy[i];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double sum = 0.0;
        for (int j = i + 1; j < n; ++j) sum += T[i][j] * x_gauss[j];
        x_gauss[i] = (d_copy[i] - sum) / T[i][i];
    }

    double diff_gauss = 0.0;
    for (int i = 0; i < n; ++i) diff_gauss = max(diff_gauss, abs(x[i] - x_gauss[i]));

    cout << "Difference from the Gauss method: " << scientific << setprecision(3) << diff_gauss << endl;

    return 0;
}