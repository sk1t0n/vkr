#include "math.h"
#include "Helper.h"

int **matr_idx = NULL;
double **A1 = NULL, **A2 = NULL;
double *_w = NULL;
double **w = NULL;

void alloc_arrays(int n, int max_idx, int** &matr_idx,
    double** &A1, double** &A2, double* &_w, double** &w) {
	matr_idx = new int*[n + 1];
	w = new double*[n + 1];
	for (int i = 0; i < n + 1; i++) {
		matr_idx[i] = new int[n + 1]();
		w[i] = new double[n + 1]();
	}
	A1 = new double*[max_idx];
	A2 = new double*[max_idx];
	for (int i = 0; i < max_idx; i++) {
		A1[i] = new double[max_idx]();
		A2[i] = new double[max_idx]();
	}
	_w = new double[max_idx];
}

void release_arrays(int n, int max_idx, int** &matr_idx,
    double** &A1, double** &A2, double* &_w, double** &w) {
	for (int i = 0; i < n + 1; i++) {
		delete[] matr_idx[i];
		delete[] w[i];
	}
	delete[] matr_idx;
	delete[] w;
	for (int i = 0; i < max_idx; i++) {
		delete[] A1[i];
		delete[] A2[i];
	}
	delete[] A1;
	delete[] A2;
	delete[] _w;
}

void init_w(int n, int max_idx, double q0, double D,
    double** &A2, double* &_w, double** &w) {
	double s;
	int k = 0;
	for (int i = 0; i < max_idx; i++) {
		s = 0;
		for (int j = 0; j < max_idx; j++) {
            s += A2[i][j] * (q0 / D);
		}
		_w[i] = s;
	}
	for (int i = 2; i <= n - 2; i++) {
        for (int j = 2; j <= n - 2; j++) {
            w[i][j] = _w[k++];
        }
	}
}

int main(int argc, char* argv[]) {
	int n;
	double a, b, h, q0, E, nu;

	Helper helper;
	helper.set_args(argc, argv, n, a, b, h, q0, E, nu);

    double D = E * pow(h, 3) / (12 * (1 - pow(nu,2)));
	int max_idx = (n - 3) * (n - 3);
	double hx = a / n, hy = b / n;

	alloc_arrays(n, max_idx, matr_idx, A1, A2, _w, w);

	helper.init_matr_idx(n, matr_idx);
    helper.init_A1(n, hx, hy, matr_idx, A1);
	helper.init_A2(max_idx, A2);
	helper.front_process(max_idx, A1, A2);
	helper.back_process(max_idx, A1, A2);

	init_w(n, max_idx, q0, D, A2, _w, w);
	helper.save_matrix_to_json(n, w, "w");

	release_arrays(n, max_idx, matr_idx, A1, A2, _w, w);
	return 0;
}
