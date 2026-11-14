#include <app.hpp>
using namespace std;
int countDigits(string s) {
    string d;
    for (char c : s)
        if (isdigit((unsigned char)c)) d += c;
    size_t pos = d.find_first_not_of('0');
    return pos == string::npos ? 0 : (int)(d.size() - pos);
}

int topExp(double x) { return (int)floor(log10(fabs(x))); }

bool readNumber(const string& prompt, string& text, double& value) {
    cout << prompt;
    cin >> text;
    replace(text.begin(), text.end(), ',', '.');
    try {
        size_t used;
        value = stod(text, &used);
        return used == text.size();
    }
    catch (...) {
        return false;
    }
}

int main() {
    char again;
    do {
        int mode;
        cout << "\n1 - find n from the error, 2 - find the error from n: ";
        cin >> mode;

        if (mode != 1 && mode != 2) {
            cout << "Wrong choice, enter 1 or 2." << endl;
        }
        else {
            string xs;
            double x;
            if (!readNumber("Approximate number: ", xs, x) || x == 0) {
                cout << "Wrong number (must be a non-zero number)." << endl;
            }
            else {
                int e = topExp(x);
                int digits = countDigits(xs);

                if (mode == 1) {
                    string ds;
                    double delta;
                    if (!readNumber("Limiting absolute error: ", ds, delta) || delta <= 0) {
                        cout << "Wrong error (must be a positive number)." << endl;
                    }
                    else {
                        int kmin = (int)ceil(log10(2 * delta) - 1e-12);
                        int n = e - kmin + 1;
                        n = max(0, min(n, digits));
                        cout << "Number of correct digits: " << n << endl;
                    }
                }
                else {
                    int n;
                    cout << "Number of correct digits (1.." << digits << "): ";
                    cin >> n;
                    if (n < 1 || n > digits) {
                        cout << "n must be between 1 and " << digits
                            << " (the number has only " << digits << " digits)." << endl;
                    }
                    else {
                        double d = 0.5 * pow(10.0, e - n + 1);
                        cout << "Limiting absolute error: " << d << endl;
                        cout << "Limiting relative error: " << d / fabs(x) << endl;
                    }
                }
            }
        }

        if (!cin) {
            cin.clear();
            cin.ignore(10000, '\n');
        }
        cout << "Repeat? (y/n): ";
        cin >> again;
    } while (again == 'y' || again == 'Y');
    return 0;
}