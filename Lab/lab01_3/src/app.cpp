#include <app.hpp>
using namespace std;

double fExact(double x) {
    return 0.25 * log((1 + x) / (1 - x)) + 0.5 * atan(x);
}

double seriesN(double x, int n) {
    double p = x, s = 0.0;
    double x4 = pow(x, 4);
    for (int i = 0; i <= n; i++) {
        s += p / (4 * i + 1);
        p *= x4;
    }
    return s;
}

double seriesEps(double x, double eps, int& cnt) {
    double p = x, s = 0.0, x4 = pow(x, 4);
    int i = 0;
    while (true) {
        double t = p / (4 * i + 1);
        s += t;
        if (fabs(t) < eps) { cnt = i + 1; return s; }
        p *= x4;
        i++;
    }
}

int main() {
    double a = 0.0, b = 0.9, eps = 1e-4;
    int k = 10, n = 8;
    cout << "Enter a, b, k, n, eps: ";
    cin >> a >> b >> k >> n >> eps;
    double h = (b - a) / (k - 1);

    cout << fixed << setprecision(6);
    cout << setw(6) << "x" << setw(14) << "S(n)" << setw(14) << "S(eps)"
        << setw(9) << "terms" << setw(14) << "exact" << endl;
    for (int j = 0; j < k; j++) {
        double x = a + j * h;
        int cnt;
        double s1 = seriesN(x, n);
        double s2 = seriesEps(x, eps, cnt);
        cout << setw(6) << setprecision(2) << x << setprecision(6)
            << setw(14) << s1 << setw(14) << s2 << setw(9) << cnt
            << setw(14) << fExact(x) << endl;
    }
    return 0;
}