#include "math.h"
#include "stdlib.h"
#include "Helper.h"

const double beta = 1, delta = 1;

int **matr_idx = NULL;
double **A1 = NULL, **A2 = NULL;
double *vec_tmp = NULL, *_w = NULL, *_r = NULL, *_r1 = NULL;
double **w = NULL, **r = NULL;

void alloc_arrays(int n, int max_idx, int** &matr_idx, double** &A1, double** &A2,
    double* &vec_tmp, double* &_w, double* &_r, double* &_r1, double** &w, double** &r) {
	matr_idx = new int*[n + 1];
	w = new double*[n + 1];
	r = new double*[n + 1];
	for (int i = 0; i < n + 1; i++) {
		matr_idx[i] = new int[n + 1]();
		w[i] = new double[n + 1]();
		r[i] = new double[n + 1]();
	}
	A1 = new double*[max_idx];
	A2 = new double*[max_idx];
	for (int i = 0; i < max_idx; i++) {
		A1[i] = new double[max_idx]();
		A2[i] = new double[max_idx]();
	}
	vec_tmp = new double[max_idx];
	_w = new double[max_idx];
	_r = new double[max_idx]();
	_r1 = new double[max_idx]();
}

void release_arrays(int n, int max_idx, int** &matr_idx, double** &A1, double** &A2,
    double* &vec_tmp, double* &_w, double* &_r, double* &_r1, double** &w, double** &r) {
	for (int i = 0; i < n + 1; i++) {
		delete[] matr_idx[i];
		delete[] w[i];
		delete[] r[i];
	}
	delete[] matr_idx;
	delete[] w;
	delete[] r;
	for (int i = 0; i < max_idx; i++) {
		delete[] A1[i];
		delete[] A2[i];
	}
	delete[] A1;
	delete[] A2;
	delete[] vec_tmp;
	delete[] _w;
	delete[] _r;
	delete[] _r1;
}

void calc_vec_w(int i, int max_idx, double D, double q0,
    double** &A2, double* &_w, double* &_r) {
	double tmp = 0;
    for (int j = 0; j < max_idx; j++) {
        tmp += A2[i][j] * (q0 - _r[j]);
    }
    tmp *= 1 / D;
    _w[i] = tmp;
}

void calc_w_r(double* &_r, double* &_r1, double* &_w, int i, int max_idx,
        double D, double q0, double** &A2) {
    _r[i] = _r1[i] + beta * (_w[i] - delta);
    if (_r[i] < 0) { _r[i] = 0; }
    calc_vec_w(i, max_idx, D, q0, A2, _w, _r);
    _r1[i] = _r[i];
}

void generalized_reaction_method(int num_iters, int max_idx, double D, double q0,
    double** &A2, double* &vec_tmp, double* &_r, double* &_r1, double* &_w) {
	for (int i = 0; i < max_idx; i++) {
        calc_vec_w(i, max_idx, D, q0, A2, _w, _r);
        vec_tmp[i] = _w[i];
	}
	if (num_iters == 0) {
        int kol = 0;
        double norma = 1, sum1 = 0, sum2 = 0, eps = 1e-3;
        while (norma >= eps) {
            sum1 = 0; sum2 = 0;
            for (int i = 0; i < max_idx; i++) {
                calc_w_r(_r, _r1, _w, i, max_idx, D, q0, A2);
                sum1 += pow(_w[i] - vec_tmp[i], 2);
                sum2 += pow(_w[i], 2);
                norma = sqrt(sum1) / sqrt(sum2);
                vec_tmp[i] = _w[i];
            }
            kol++;
        }
        cout << "kol = " << kol << endl;
	} else {
        for (int k = 1; k <= num_iters; k++) {
            for (int i = 0; i < max_idx; i++) {
                calc_w_r(_r, _r1, _w, i, max_idx, D, q0, A2);
            }
        }
	}
}

int main(int argc, char* argv[]) {
	int n, num_iters;
	double a, b, h, q0, E, nu;

	Helper helper;
	helper.set_args(argc, argv, n, a, b, h, q0, E, nu);
	if (argc > 8) { num_iters = atoi(argv[8]); }

	double D = E * pow(h, 3) / (12 * (1 - pow(nu,2)));
	int max_idx = (n - 3) * (n - 3);
	double hx = a / n, hy = b / n;

	alloc_arrays(n, max_idx, matr_idx, A1, A2, vec_tmp, _w, _r, _r1, w, r);

	helper.init_matr_idx(n, matr_idx);
    helper.init_A1(n, hx, hy, matr_idx, A1);
	helper.init_A2(max_idx, A2);
	helper.front_process(max_idx, A1, A2);
	helper.back_process(max_idx, A1, A2);

	generalized_reaction_method(num_iters, max_idx, D, q0, A2, vec_tmp, _r, _r1, _w);

	helper.init_w_r(n, _w, _r, w, r);
	helper.save_matrix_to_json(n, w, "w");
	helper.save_matrix_to_json(n, r, "r");

	release_arrays(n, max_idx, matr_idx, A1, A2, vec_tmp, _w, _r, _r1, w, r);
	return 0;
}
