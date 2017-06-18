#include "math.h"
#include "stdlib.h"
#include "iostream"
#include "Helper.h"

using namespace std;

const double beta = 1e-6;
const double delta = 1;

int **matr_idx = NULL;
double **A1 = NULL, **A2 = NULL;
double **Axx = NULL, **Ayy = NULL, **Axy = NULL;
double *r = NULL, *r_prev = NULL, *fi = NULL, *w = NULL;
double **matr_tmp1 = NULL, **matr_tmp2 = NULL, **matr_tmp3 = NULL,
       **matr_tmp4 = NULL, **matr_tmp5 = NULL, **matr_tmp6 = NULL,
       **matr_tmp7 = NULL, **matr_tmp8 = NULL, **matr_tmp9 = NULL;
double *vec_tmp = NULL;
double **W = NULL, **R = NULL;

void alloc_arrays(int n, int max_idx,
    double** &Axx, double** &Ayy, double** &Axy, int** &matr_idx, double* &vec_tmp,
    double** &A1, double** &A2, double* &r, double* &r_prev, double* &fi, double* &w,
    double** &matr_tmp1, double** &matr_tmp2, double** &matr_tmp3,
    double** &matr_tmp4, double** &matr_tmp5, double** &matr_tmp6,
    double** &matr_tmp7, double** &matr_tmp8, double** &matr_tmp9, double** &W, double** &R) {

	matr_tmp1 = new double*[n - 3];
	matr_tmp2 = new double*[n - 3];
	matr_tmp3 = new double*[n - 3];
	matr_tmp4 = new double*[n - 3];
	matr_tmp5 = new double*[n - 3];
	matr_tmp6 = new double*[n - 3];
	matr_tmp7 = new double*[n - 3];
	matr_tmp8 = new double*[n - 3];
	matr_tmp9 = new double*[n - 3];
	for (int i = 0; i < n - 3; i++) {
        matr_tmp1[i] = new double[n - 3]();
        matr_tmp2[i] = new double[n - 3]();
        matr_tmp3[i] = new double[n - 3]();
        matr_tmp4[i] = new double[n - 3]();
        matr_tmp5[i] = new double[n - 3]();
        matr_tmp6[i] = new double[n - 3]();
        matr_tmp7[i] = new double[n - 3]();
        matr_tmp8[i] = new double[n - 3]();
        matr_tmp9[i] = new double[n - 3]();
	}

	matr_idx = new int*[n + 1];
	W = new double*[n + 1];
	R = new double*[n + 1];
	for (int i = 0; i < n + 1; i++) {
		matr_idx[i] = new int[n + 1]();
		W[i] = new double[n + 1]();
		R[i] = new double[n + 1]();
	}

	r = new double[max_idx]();
	r_prev = new double[max_idx];
	fi = new double[max_idx]();
	w = new double[max_idx]();
	vec_tmp = new double[max_idx]();

	A1 = new double*[max_idx];
	A2 = new double*[max_idx];
	Axx = new double*[max_idx];
	Ayy = new double*[max_idx];
	Axy = new double*[max_idx];
	for (int i = 0; i < max_idx; i++) {
		A1[i] = new double[max_idx]();
		A2[i] = new double[max_idx]();
		Axx[i] = new double[max_idx]();
        Ayy[i] = new double[max_idx]();
        Axy[i] = new double[max_idx]();
	}
}

void release_arrays(int n, int max_idx,
    double** &Axx, double** &Ayy, double** &Axy, int** &matr_idx, double* &vec_tmp,
    double** &A1, double** &A2, double* &r, double* &r_prev, double* &fi, double* &w,
    double** &matr_tmp1, double** &matr_tmp2, double** &matr_tmp3,
    double** &matr_tmp4, double** &matr_tmp5, double** &matr_tmp6,
    double** &matr_tmp7, double** &matr_tmp8, double** &matr_tmp9, double** &W, double** &R) {

    delete[] r;
    delete[] r_prev;
    delete[] fi;
    delete[] w;

    for (int i = 0; i <  n - 3; i++) {
		delete[] matr_tmp1[i];
		delete[] matr_tmp2[i];
		delete[] matr_tmp3[i];
		delete[] matr_tmp4[i];
		delete[] matr_tmp5[i];
		delete[] matr_tmp6[i];
		delete[] matr_tmp7[i];
		delete[] matr_tmp8[i];
		delete[] matr_tmp9[i];
	}
    delete[] matr_tmp1;
    delete[] matr_tmp2;
    delete[] matr_tmp3;
    delete[] matr_tmp4;
    delete[] matr_tmp5;
    delete[] matr_tmp6;
    delete[] matr_tmp7;
    delete[] matr_tmp8;
    delete[] matr_tmp9;

    for (int i = 0; i < n + 1; i++) {
        delete[] matr_idx[i];
        delete[] W[i];
        delete[] R[i];
    }
    delete[] matr_idx;
    delete[] W;
    delete[] R;

    for (int i = 0; i < max_idx; i++) {
        delete[] A1[i];
        delete[] A2[i];
        delete[] Axx[i];
		delete[] Ayy[i];
		delete[] Axy[i];
    }
    delete[] A1;
    delete[] A2;
    delete[] Axx;
	delete[] Ayy;
	delete[] Axy;
}

void vector_to_matrix(int max_idx, int n, double* &v, double** &M) {
    int k = 0;
    for (int i = 0; i < n - 3; i++) {
        for (int j = 0; j < n - 3; j++) {
            M[i][j] = v[k++];
        }
    }
}

void matrix_mul_matrix(int n, double** &M1, double** &M2, double** &matr_tmp, double c = 1) {
    double tmp;
    for (int i = 0; i < n - 3; i++) {
        for (int j = 0; j < n - 3; j++) {
            tmp = 0;
            for (int k = 0; k < n - 3; k++) {
                tmp += M1[i][k] * M2[k][j];
            }
            matr_tmp[i][j] = c * tmp;
        }
    }
}

void matrix_mul_vector(int max_idx, double** &M, double* &v, double* &vec_tmp, double c = 1, double q0 = 0) {
    double tmp;
    for (int i = 0; i < max_idx; i++) {
		tmp = 0;
		for (int j = 0; j < max_idx; j++) {
			if (q0 == 0) {
                tmp += M[i][j] * v[j];
			} else {
                tmp += M[i][j] * q0;
			}
		}
		vec_tmp[i] = c * tmp;
	}
}

void calc_vec_fi(int max_idx, int n, double E, double h, double** &A2,
    double** &matr_tmp1, double** &matr_tmp2, double** &matr_tmp3, double** &matr_tmp4, double** &matr_tmp5,
    double** &Axx, double** &Ayy, double** &Axy, double* &w, double* &fi, double* &vec_tmp) {

    matrix_mul_vector(max_idx, Axy, w, vec_tmp);
    vector_to_matrix(max_idx, n, vec_tmp, matr_tmp1);
    matrix_mul_matrix(n, matr_tmp1, matr_tmp1, matr_tmp4);

    matrix_mul_vector(max_idx, Axx, w, vec_tmp);
    vector_to_matrix(max_idx, n, vec_tmp, matr_tmp2);
    matrix_mul_vector(max_idx, Ayy, w, vec_tmp);
    vector_to_matrix(max_idx, n, vec_tmp, matr_tmp3);
    matrix_mul_matrix(n, matr_tmp2, matr_tmp3, matr_tmp5);

    int k = 0;
    for (int i = 0; i < n - 3; i++) {
        for (int j = 0; j < n - 3; j++) {
            vec_tmp[k++] = matr_tmp4[i][j] - matr_tmp5[i][j];
        }
    }

    matrix_mul_vector(max_idx, A2, vec_tmp, fi, E * h);
}

void calc_vec_w(int max_idx, int n, double D, double q0, double** &A2,
    double** &matr_tmp1, double** &matr_tmp2, double** &matr_tmp3, double** &matr_tmp4,
    double** &matr_tmp5, double** &matr_tmp6, double** &matr_tmp7, double** &matr_tmp8, double** &matr_tmp9,
    double** &Axx, double** &Ayy, double** &Axy, double* &w, double* &fi, double* &vec_tmp) {

    matrix_mul_vector(max_idx, Ayy, fi, vec_tmp);
    vector_to_matrix(max_idx, n, vec_tmp, matr_tmp4);
    matrix_mul_vector(max_idx, Axx, fi, vec_tmp);
    vector_to_matrix(max_idx, n, vec_tmp, matr_tmp5);
    matrix_mul_vector(max_idx, Axy, fi, vec_tmp);
    vector_to_matrix(max_idx, n, vec_tmp, matr_tmp6);

    matrix_mul_matrix(n, matr_tmp4, matr_tmp2, matr_tmp7);
    matrix_mul_matrix(n, matr_tmp5, matr_tmp3, matr_tmp8);
    matrix_mul_matrix(n, matr_tmp6, matr_tmp1, matr_tmp9);

    int k = 0;
    for (int i = 0; i < n - 3; i++) {
        for (int j = 0; j < n - 3; j++) {
            vec_tmp[k++] = (q0 - r[k]) + matr_tmp7[i][j] + matr_tmp8[i][j] - 2 * matr_tmp9[i][j];
        }
    }

    matrix_mul_vector(max_idx, A2, vec_tmp, w, 1 / D);
}

void generalized_reaction_method(int num_iters, int max_idx, int n,
    double q0, double D, double E, double h, double** &A2, double* &w,
    double* &r, double* &r_prev, double* &fi, double* &vec_tmp,
    double** &matr_tmp1, double** &matr_tmp2, double** &matr_tmp3,
    double** &matr_tmp4, double** &matr_tmp5, double** &matr_tmp6,
    double** &matr_tmp7, double** &matr_tmp8, double** &matr_tmp9,
    double** &Axx, double** &Ayy, double** &Axy) {

    matrix_mul_vector(max_idx, A2, vec_tmp, w, 1 / D, q0);
    cout << "w0:" << endl;
	for (int i = 0; i < max_idx; i++) {
        cout << w[i] << " ";
	}
	cout << endl << endl;
	cout << "r0:" << endl;
	for (int i = 0; i < max_idx; i++) {
        cout << r[i] << " ";
	}
	cout << endl << endl;

	for (int k = 1; k <= num_iters; k++) {
        for (int i = 0; i < max_idx; i++) {
            r[i] = r_prev[i] + beta * (w[i] - delta);
            if (r[i] < 0) { r[i] = 0; }
            r_prev[i] = r[i];
        }

        calc_vec_fi(max_idx, n, E, h, A2, matr_tmp1, matr_tmp2, matr_tmp3,
            matr_tmp4, matr_tmp5, Axx, Ayy, Axy, w, fi, vec_tmp);
        calc_vec_w(max_idx, n, D, q0, A2, matr_tmp1, matr_tmp2, matr_tmp3, matr_tmp4,
            matr_tmp5, matr_tmp6, matr_tmp7, matr_tmp8, matr_tmp9,Axx, Ayy, Axy, w, fi, vec_tmp);

        cout << "r" << k << endl;
        for (int i = 0; i < max_idx; i++) {
            cout << r[i] << " ";
        }
        cout << endl << endl;
        cout << "fi" << k << endl;
        for (int i = 0; i < max_idx; i++) {
            cout << fi[i] << " ";
        }
        cout << endl << endl;
        cout << "w" << k << endl;
        for (int i = 0; i < max_idx; i++) {
            cout << w[i] << " ";
        }
        cout << endl << endl;
	}
}

int main(int argc, char* argv[]) {
	int n = 10;
	double a = 100, b = 100;
	double h = 1;
	double q0 = 10;
	double E = 2 * pow(10, 6);
    double nu = 0.3;
    int num_iters = 5;

    Helper helper;
	helper.set_args(argc, argv, n, a, b, h, q0, E, nu);
	if (argc > 8) { num_iters = atoi(argv[8]); }

	double D = E * pow(h, 3) / (12 * (1 - pow(nu,2)));
	int max_idx = (n - 3) * (n - 3);
	double hx = a / n, hy = b / n;

	alloc_arrays(n, max_idx, Axx, Ayy, Axy,
                 matr_idx, vec_tmp, A1, A2, r, r_prev, fi, w, matr_tmp1, matr_tmp2, matr_tmp3,
                 matr_tmp4, matr_tmp5, matr_tmp6, matr_tmp7, matr_tmp8, matr_tmp9, W, R);

	helper.init_matr_idx(n, matr_idx);
    helper.init_A1(n, hx, hy, matr_idx, A1);
	helper.init_A2(max_idx, A2);
	helper.init_Axx_Ayy_Axy(n, hx, hy, matr_idx, Axx, Ayy, Axy);
	helper.front_process(max_idx, A1, A2);
	helper.back_process(max_idx, A1, A2);

	generalized_reaction_method(num_iters, max_idx, n, q0, D, E, h, A2, w, r, r_prev, fi, vec_tmp,
       matr_tmp1, matr_tmp2, matr_tmp3, matr_tmp4, matr_tmp5, matr_tmp6, matr_tmp7,
       matr_tmp8, matr_tmp9, Axx, Ayy, Axy);

	helper.init_w_r(n, w, r, W, R);
	helper.save_matrix_to_json(n, W, "w");
	helper.save_matrix_to_json(n, R, "r");

	release_arrays(n, max_idx, Axx, Ayy, Axy,
                 matr_idx, vec_tmp, A1, A2, r, r_prev, fi, w, matr_tmp1, matr_tmp2, matr_tmp3,
                 matr_tmp4, matr_tmp5, matr_tmp6, matr_tmp7, matr_tmp8, matr_tmp9, W, R);
	return 0;
}
