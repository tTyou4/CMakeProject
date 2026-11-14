#include <app.hpp>
using namespace std;

int main() {
    const double x = 10.0;
    char again;
    do {
        double eps;
        cout << "Enter accuracy eps (0.1 ... 1e-6): ";
        cin >> eps;

        double a = 1.0, s = 1.0;
        int n = 0;
        while (true) {
            n++;
            a = a * x / n;    
            s += a;
            if (n > x && a < eps) break;
        }
        double ex = exp(x);
        cout << fixed << setprecision(8)
            << "S = " << s << ", terms: " << n + 1
            << ", e^10 = " << ex << ", |error| = "
            << scientific << setprecision(2) << fabs(s - ex) << endl;

        cout << "Repeat? (y/n): ";
        cin >> again;
    } while (again == 'y' || again == 'Y');
    return 0;
}