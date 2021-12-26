#include "capd/capdlib.h"
#include <vector>

using namespace std;
using namespace capd;

template<class T>
void computeODECoeffs(vector <vector<T>> &w, vector <vector<T>> &coeffs, int order) { 
    //x^2, y^2, 2xy, -x-2xy, y^2-y, y^2-y-x^2
    for (int i = 0; i < order; i++) {
        // w[0] = x^2
        for (int j = 0; j < i + 1; j++) {
            w[i][0] += coeffs[j][0] * coeffs[i - j][0];
        }

        // w[1] = y^2
        for (int j = 0; j < i + 1; j++) {
            w[i][1] += coeffs[j][1] * coeffs[i - j][1];
        }

        // w[2] = 2xy
        for (int j = 0; j < i + 1; j++) {
            w[i][2] += coeffs[j][0] * coeffs[i - j][1];
        }
        w[i][2] = 2 * w[i][2];

        // w[3] = -x-2xy
        w[i][3] = (-1) * coeffs[i][0] - w[i][2];

        //w[4] = y^2-y
        w[i][4] = w[i][1] - coeffs[i][1];

        // w[5] = w[4]-w[0] = y^2-y-x^2
        w[i][5] = w[i][4] - w[i][0];

        T t = T(i+1);
        coeffs[i + 1][0] = coeffs[i][2] / t;
        coeffs[i + 1][1] = coeffs[i][3] / t;
        coeffs[i + 1][2] = w[i][3] / t;
        coeffs[i + 1][3] = w[i][5] / t;
    }
}

inline void computeEnclosureODECoeffs(vector <vector<interval>> &w, vector <vector<interval>> &coeffs, int order) { 
    computeODECoeffs(w,coeffs,order);
}

template<class T>
void computeODECoeffs(vector <vector<T>> &w, vector <vector<T>> &coeffs, vector <vector<T>> &M, vector <vector<T>> &dcoeffs, int order) {
    for (int i = 0; i < order; i++) {

        // w[0] = x^2
        for (int j = 0; j < i + 1; j++) {
            w[i][0] += coeffs[j][0] * coeffs[i - j][0];
            for (int k = 0; k < 4; k++) {
                M[i * 4 + k][0] += dcoeffs[4 * j][k] * coeffs[i - j][0] + coeffs[j][0] * dcoeffs[4 * (i - j)][k];
            }
        }

        // w[1] = y^2
        for (int j = 0; j < i + 1; j++) {
            w[i][1] += coeffs[j][1] * coeffs[i - j][1];
            for (int k = 0; k < 4; k++) {
                M[i * 4 + k][1] +=
                        dcoeffs[4 * j + 1][k] * coeffs[i - j][1] + coeffs[j][1] * dcoeffs[4 * (i - j) + 1][k];
            }
        }

        // w[2] = 2xy
        for (int j = 0; j < i + 1; j++) {
            w[i][2] += coeffs[j][0] * coeffs[i - j][1];
            for (int k = 0; k < 4; k++) {
                M[i * 4 + k][2] += dcoeffs[4 * j][k] * coeffs[i - j][1] + coeffs[j][0] * dcoeffs[4 * (i - j) + 1][k];
            }
        }
        w[i][2] = 2 * w[i][2];
        for (int k = 0; k < 4; k++) {
            M[i * 4 + k][2] = 2 * M[i * 4 + k][2];
        }

        // w[3] = -x-2xy
        w[i][3] = (-1) * coeffs[i][0] - w[i][2];
        for (int k = 0; k < 4; k++) {
            M[i * 4 + k][3] = (-1) * dcoeffs[4 * i][k] - M[4 * i + k][2];
        }

        //w[4] = y^2-y
        w[i][4] = w[i][1] - coeffs[i][1];
        for (int k = 0; k < 4; k++) {
            M[i * 4 + k][4] = M[4 * i + k][1] - dcoeffs[4 * i + 1][k];
        }

        // w[5] = w[4]-w[0] = y^2-y-x^2
        w[i][5] = w[i][4] - w[i][0];
        for (int k = 0; k < 4; k++) {
            M[4 * i + k][5] = M[4 * i + k][4] - M[4 * i + k][0];
        }

        T t = T(i+1);
        coeffs[i + 1][0] = coeffs[i][2] / t;
        coeffs[i + 1][1] = coeffs[i][3] / t;
        coeffs[i + 1][2] = w[i][3] / t;
        coeffs[i + 1][3] = w[i][5] / t;

        for (int k = 0; k < 4; k++) {
            dcoeffs[4 * (i + 1)][k] = dcoeffs[4 * i + 2][k] / t;
            dcoeffs[4 * (i + 1) + 1][k] = dcoeffs[4 * i + 3][k] / t;
            dcoeffs[4 * (i + 1) + 2][k] = M[4 * i + k][3] / t;
            dcoeffs[4 * (i + 1) + 3][k] = M[4 * i + k][5] / t;
        }
    }
}

inline void
computeEnclosureODECoeffs(vector <vector<interval>> &w, vector <vector<interval>> &coeffs, vector <vector<interval>> &M,
                          vector <vector<interval>> &dcoeffs, int order) {
    computeODECoeffs(w,coeffs,M,dcoeffs,order);
}


template<class T>
vector <T> taylor_poly(vector <vector<T>> &coeffs, interval h, int order) {

    vector <vector<T>> res(order + 1, vector<T>(4));

    for (int j = 0; j < 4; j++) {
        res[0][j] = coeffs[order][j];
        for (int i = 1; i < order + 1; i++) {
            res[i][j] = res[i - 1][j] * T(h) + coeffs[order - i][j];
        }
    }

    return res[order];
}

template<class T>
vector <vector<T>> taylor_poly_d(vector <vector<T>> &dcoeffs, interval h, int order) {

    vector <vector<T>> res(4 * (order + 1), vector<T>(4));

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res[i][j] = dcoeffs[4 * order + i][j];
        }
    }

    for (int i = 1; i < order + 1; i++) {
        for (int j = 0; j < 4; j++) {
            res[4 * i][j] = res[4 * (i - 1)][j] * T(h) + dcoeffs[4 * order - 4 * i][j];
            res[4 * i + 1][j] = res[4 * (i - 1) + 1][j] * T(h) + dcoeffs[4 * order - 4 * i + 1][j];
            res[4 * i + 2][j] = res[4 * (i - 1) + 2][j] * T(h) + dcoeffs[4 * order - 4 * i + 2][j];
            res[4 * i + 3][j] = res[4 * (i - 1) + 3][j] * T(h) + dcoeffs[4 * order - 4 * i + 3][j];
        }
    }

    vector <vector<T>> result(4, vector<T>(4));

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            result[i][j] = res[4 * order + i][j];
        }
    }
    return result;
}


template<typename VectorType>
VectorType
high_C0_enclose_1(VectorType u, double current_time, double time_step, int order, VectorType &previous_remainder) {

    vector <vector<interval>> x_tmp(order + 1, vector<interval>(6));
    vector <vector<interval>> x_coeffs(order + 2, vector<interval>(4));

    for (int i = 0; i < 4; i++) {
        x_coeffs[0][i] = u[i];
    }

    computeEnclosureODECoeffs(x_tmp, x_coeffs, order + 1); //[X]^1,...,[X]^(r+1)
    cout << "[X]^1,...,[X]^(r+1)" << endl;
    for (int i = 0; i < order + 2; i++) {
        for (int j = 0; j < 4; j++) {
            cout << x_coeffs[i][j];
        }
        cout << endl;
    }
    cout << endl;

    vector <interval> x_res(4);
    x_res = taylor_poly(x_coeffs, interval(0, time_step), order);    //[P]=[X]+[0,h][X]^1+...+[0,h]^r[X]^r

    cout << "[P]=[X]+[0,h][X]^1+...+[0,h]^r[X]^r" << endl;
    for (int i = 0; i < 4; i++) {
        cout << x_res[i];
    }
    cout << endl;
    cout << endl;


    double eps = 1e-30;
    interval factor(-2, 2);
    interval h = power(interval(0, time_step), order + 1);             

    cout << "[R]" << endl;
    vector <interval> R(4);
    if (current_time == 0) {
        for (int i = 0; i < 4; i++) {
            R[i] = factor * (x_coeffs[order + 1][i] * h + eps); //[R]=[-2,2]*[X]^(r+1)*[0,h]^(r+1)+[-2,2]*eps
            cout << R[i];
        }
    } else {
        for (int i = 0; i < 4; i++) {
            R[i] = factor * previous_remainder[i];  //[R]=[-2,2]*[0,h]^(r+1)[Y]^(r+1)
            cout << R[i];
        }
    }
    cout << endl;
    cout << endl;

    vector <interval> Y(4), Z(4);
    vector <vector<interval>> y_tmp(order + 1, vector<interval>(6));
    vector <vector<interval>> y_coeffs(order + 2, vector<interval>(4));

    int counter = 30;
    bool found = false;
    while (!found) {
        counter--;
        if (counter == 0) {
            throw "Nie mozna znalezc enclosure!";
        }

        cout << "[Y]=[P]+[R]" << endl;
        for (int i = 0; i < 4; i++) {
            Y[i] = x_res[i] + R[i]; //[Y] = [P]+[R]
            cout << Y[i];
        }
        cout << endl;
        cout << endl;


        for (int i = 0; i < 4; i++) {
            y_coeffs[0][i] = Y[i];
        }

        for (int i = 0; i < order + 1; i++) {
            for (int j = 0; j < 6; j++) {
                y_tmp[i][j] = interval(0);
            }

            for (int j = 0; j < 4; j++) {
                y_coeffs[i + 1][j] = interval(0);
            }
        }

        computeEnclosureODECoeffs(y_tmp, y_coeffs, order + 1); //[Y]^1,..,[Y]^(r+1)
        cout << "[Y], [Y]^1,...,[Y]^(r+1)" << endl;
        for (int i = 0; i < order + 2; i++) {
            for (int j = 0; j < 4; j++) {
                cout << y_coeffs[i][j];
            }
            cout << endl;
        }
        cout << endl;


        cout << "[Y]^(r+1)[0,h]^(r+1)" << endl;
        for (int i = 0; i < 4; i++) {
            Z[i] = y_coeffs[order + 1][i] * h; //[Y]^(r+1)[0,h]^(r+1)
            cout << Z[i];
        }
        cout << endl;

        found = true;
        /// tutaj jest inna strategia - przechodzimy po wszystkich 4-zmiennych i poprawiamy Y, jesli trzeba
        /// Dopiero po poprawieniu wszystkich wspolrzednych obliczamy ponownie wspolczynniki (jesli trzeba)
        for (int i = 0; i < 4; i++) {
            if (!Z[i].subset(R[i])) {
                cout << "Nie dziala ograniczenie dla wspolrzedniej: " << i << endl;
                cout << "[Y]^(r+1)[0,h]^(r+1): " << Z[i] << endl;
                cout << "[R]: " << R[i] << endl;
                found = false;
                R[i] = intervalHull(Z[i], R[i]) * factor;
            }
        }
    }

    // Update previous remainder
    for (int i = 0; i < 4; i++) {
        previous_remainder[i] = Z[i];
    }

    // Result
    vector <interval> xx_res(4);
    xx_res = taylor_poly(x_coeffs, interval(time_step, time_step), order);

    VectorType res(4);
    for (int i = 0; i < 4; i++) {
        res[i] = xx_res[i] + Z[i];
    }
    cout << res << endl;
    return res;
}

template<typename MatrixType, typename VectorType>
MatrixType high_C1_enclose_1(VectorType u, double current_time, double time_step, int order, VectorType enc,
                             MatrixType &previous_remainder) {

    // calculate V+[0,h]V^1+...+[0,h]V^r
    vector <vector<interval>> d_tmp(4 * (order + 1), vector<interval>(6));
    vector <vector<interval>> d_coeffs(4 * (order + 2), vector<interval>(4));

    vector <vector<interval>> x_tmp(order + 1, vector<interval>(6));
    vector <vector<interval>> x_coeffs(order + 2, vector<interval>(4));
    for (int i = 0; i < 4; i++) {
        d_coeffs[i][i] = interval(1);
        x_coeffs[0][i] = u[i];
    }

    computeEnclosureODECoeffs(x_tmp, x_coeffs, d_tmp, d_coeffs, order + 1);

    vector <vector<interval>> d_res(4, vector<interval>(4));

    d_res = taylor_poly_d(d_coeffs, interval(0, time_step), order); //[P]

    double eps = 1e-30;
    interval factor(-2, 2);
    interval h = power(interval(0, time_step), order + 1);

    vector <vector<interval>> R(4, vector<interval>(4));
    if (current_time == 0) {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                R[i][j] = factor * (d_coeffs[4 * (order + 1) + i][j] * h + eps);
            }
        }
    } else {
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                R[i][j] = factor * previous_remainder[i][j];
            }
        }
    }

    vector <vector<interval>> Y(4, vector<interval>(4)), Z(4, vector<interval>(4));
    vector <vector<interval>> y_tmp(4 * (order + 1), vector<interval>(6));
    vector <vector<interval>> y_coeffs(4 * (order + 2), vector<interval>(4));

    int counter = 30;
    bool found = false;
    while (!found) {
        counter--;
        if (counter == 0) {
            throw "Nie mozna znalezc enclosure!";
        }

        cout << "[Y]=[P]+[R]" << endl;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Y[i][j] = d_res[i][j] + R[i][j]; //[Y] = [P]+[R]
            }
        }
        cout << endl;
        cout << endl;


        for (int i = 0; i < 4; i++) {
            x_coeffs[0][i] = enc[i];
            for (int j = 0; j < 4; j++) {
                y_coeffs[i][j] = Y[i][j];
            }
        }

        for (int i = 0; i < 4 * (order + 1); i++) {
            for (int j = 0; j < 6; j++) {
                y_tmp[i][j] = interval(0);
            }
            for (int j = 0; j < 4; j++) {
                y_coeffs[i + 4][j] = interval(0);
            }
        }

        for (int i = 0; i < order + 1; i++) {
            for (int j = 0; j < 6; j++) {
                x_tmp[i][j] = interval(0);
            }

            for (int j = 0; j < 4; j++) {
                x_coeffs[i + 1][j] = interval(0);
            }
        }


        computeEnclosureODECoeffs(x_tmp, x_coeffs, y_tmp, y_coeffs, order + 1); //[Y]^1,..,[Y]^(r+1)
        cout << "[Y], [Y]^1,...,[Y]^(r+1)" << endl;
        for (int i = 0; i < 4 * (order + 2); i++) {
            for (int j = 0; j < 4; j++) {
                cout << y_coeffs[i][j];
            }
            cout << endl;
        }
        cout << endl;


        cout << "[Y]^(r+1)[0,h]^(r+1)" << endl;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                Z[i][j] = y_coeffs[4 * (order + 1) + i][j] * h; //[Y]^(r+1)[0,h]^(r+1)
                cout << Z[i];
            }
        }
        cout << endl;

        found = true;
        for (int i = 0; i < 4; i++) {
            for (int j = 0; j < 4; j++) {
                if (!Z[i][j].subset(R[i][j])) {
                    cout << "Nie dziala ograniczenie dla wspolrzedniej: " << i << j << endl;
                    cout << "[Y]^(r+1)[0,h]^(r+1): " << Z[i][j] << endl;
                    cout << "[R]: " << R[i][j] << endl;
                    found = false;
                    R[i][j] = intervalHull(Z[i][j], R[i][j]) * factor;
                }
            }
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            previous_remainder[i][j] = Z[i][j];
        }
    }
    // Result
    vector <vector<interval>> xx_res(4, vector<interval>(4));
    xx_res = taylor_poly_d(d_coeffs, interval(time_step, time_step), order);

    MatrixType res(4, 4);
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            res[i][j] = xx_res[i][j] + Z[i][j];
        }
    }

    return res;
}

