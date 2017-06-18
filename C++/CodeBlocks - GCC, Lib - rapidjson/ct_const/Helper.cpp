#include "Helper.h"
#include "math.h"
#include "rapidjson/writer.h"
#include "rapidjson/stringbuffer.h"
#include "fstream"

using namespace rapidjson;

Helper::Helper(void) {}

Helper::~Helper(void) {}

void Helper::set_args(int argc, char* argv[], int &n, double &a, double &b,
    double &h, double &q0, double &E, double &nu) {
	if (argc > 1) { n = atoi(argv[1]); }
	if (argc > 2) { a = atof(argv[2]); }
	if (argc > 3) { b = atof(argv[3]); }
	if (argc > 4) { h = atof(argv[4]); }
	if (argc > 5) { q0 = atof(argv[5]); }
	if (argc > 6) { E = atof(argv[6]); }
	if (argc > 7) { nu = atof(argv[7]); }
}

void Helper::init_matr_idx(int n, int** &matr_idx) {
	int k = 1;
	for (int i = 2; i <= n - 2; i++) {
        for (int j = 2; j <= n - 2; j++) {
            matr_idx[i][j] = k++;
        }
	}
}

void Helper::init_A1(int n, double hx, double hy, int** &matr_idx, double** &A1) {
	double hx4 = pow(hx, 4), hy4 = pow(hy, 4), hx2hy2 = pow(hx, 2) * pow(hy, 2);
	int k = 0;
	for (int i = 2; i <= n - 2; i++) {
		for (int j = 2; j <= n - 2; j++) {
			A1[k][matr_idx[i][j] - 1] = 6 / hx4 + 6 / hy4 + 8 / hx2hy2;
			if (matr_idx[i - 2][j] != 0) { A1[k][matr_idx[i - 2][j] - 1] = 1 / hx4; }
			if (matr_idx[i + 2][j] != 0) { A1[k][matr_idx[i + 2][j] - 1] = 1 / hx4; }
			if (matr_idx[i][j - 2] != 0) { A1[k][matr_idx[i][j - 2] - 1] = 1 / hy4; }
			if (matr_idx[i][j + 2] != 0) { A1[k][matr_idx[i][j + 2] - 1] = 1 / hy4; }
			if (matr_idx[i - 1][j] != 0) { A1[k][matr_idx[i - 1][j] - 1] = -4 / hx4 -4 / hx2hy2; }
			if (matr_idx[i + 1][j] != 0) { A1[k][matr_idx[i + 1][j] - 1] = -4 / hx4 -4 / hx2hy2; }
			if (matr_idx[i][j - 1] != 0) { A1[k][matr_idx[i][j - 1] - 1] = -4 / hy4 -4 / hx2hy2; }
			if (matr_idx[i][j + 1] != 0) { A1[k][matr_idx[i][j + 1] - 1] = -4 / hy4 -4 / hx2hy2; }
			if (matr_idx[i - 1][j - 1] != 0) { A1[k][matr_idx[i - 1][j - 1] - 1] = 2 / hx2hy2; }
			if (matr_idx[i - 1][j + 1] != 0) { A1[k][matr_idx[i - 1][j + 1] - 1] = 2 / hx2hy2; }
			if (matr_idx[i + 1][j - 1] != 0) { A1[k][matr_idx[i + 1][j - 1] - 1] = 2 / hx2hy2; }
			if (matr_idx[i + 1][j + 1] != 0) { A1[k][matr_idx[i + 1][j + 1] - 1] = 2 / hx2hy2; }
			k++;
		}
	}
}

void Helper::init_A2(int max_idx, double** &A2) {
    for (int i = 0; i < max_idx; i++) {
        A2[i][i] = 1.0;
    }
}

void Helper::init_Axx_Ayy_Axy(int n, double hx, double hy, int** &matr_idx,
                              double** &Axx, double** &Ayy, double** &Axy) {
	double hx2 = pow(hx, 2), hy2 = pow(hy, 2);
	int k = 0;
	for (int i = 2; i <= n - 2; i++) {
		for (int j = 2; j <= n - 2; j++) {
			Axx[k][matr_idx[i][j] - 1] = -2 / hx2;
			Ayy[k][matr_idx[i][j] - 1] = -2 / hy2;
			Axy[k][matr_idx[i][j] - 1] = 1 / (hx * hy);
			if (matr_idx[i - 1][j] != 0) {
                Axx[k][matr_idx[i - 1][j] - 1] = 1 / hx2;
                Axy[k][matr_idx[i - 1][j] - 1] = -1 / (2 * hx * hy);
            }
			if (matr_idx[i + 1][j] != 0) {
                Axx[k][matr_idx[i + 1][j] - 1] = 1 / hx2;
                Axy[k][matr_idx[i + 1][j] - 1] = -1 / (2 * hx * hy);
            }
			if (matr_idx[i][j - 1] != 0) {
                Ayy[k][matr_idx[i][j - 1] - 1] = 1 / hy2;
                Axy[k][matr_idx[i][j - 1] - 1] = -1 / (2 * hx * hy);
            }
			if (matr_idx[i][j + 1] != 0) {
                Ayy[k][matr_idx[i][j + 1] - 1] = 1 / hy2;
                Axy[k][matr_idx[i][j + 1] - 1] = -1 / (2 * hx * hy);
            }
			if (matr_idx[i - 1][j - 1] != 0) { Axy[k][matr_idx[i - 1][j - 1] - 1] = 1 / (2 * hx * hy); }
			if (matr_idx[i + 1][j + 1] != 0) { Axy[k][matr_idx[i + 1][j + 1] - 1] = 1 / (2 * hx * hy); }
			k++;
		}
	}
}

void Helper::front_process(int max_idx, double** &A1, double** &A2) {
	double tmp;
	for (int k = 0; k < max_idx; k++) {
		if (A1[k][k] == 0.0) {
            for (int j = 0; j < max_idx; j++) {
				tmp = A1[k][j]; A1[k][j] = A1[k + 1][j]; A1[k + 1][j] = tmp;
				tmp = A2[k][j]; A2[k][j] = A2[k + 1][j]; A2[k + 1][j] = tmp;
			}
		}
		if (A1[k][k] != 0.0 && A1[k][k] != 1.0) {
			tmp = A1[k][k];
			for (int j = k; j < max_idx; j++) {
                A1[k][j] /= tmp;
			}
			for (int j = 0; j < max_idx; j++) {
                A2[k][j] /= tmp;
			}
			if (k == max_idx - 1) { break; }
		}
		for (int i = k + 1; i < max_idx; i++) {
			tmp = A1[i][k];
			for (int j = k; j < max_idx; j++) {
                A1[i][j] -= A1[k][j] * tmp;
			}
			for (int j = 0; j < max_idx; j++) {
                A2[i][j] -= A2[k][j] * tmp;
			}
		}
	}
}

void Helper::back_process(int max_idx, double** &A1, double** &A2) {
	double tmp;
	for (int k = max_idx - 1; k > 0; k--) {
        for (int i = k - 1; i > -1; i--) {
			tmp = A1[i][k];
			for (int j = k; j > i; j--) {
                A1[i][j] -= A1[k][j] * tmp;
			}
			for (int j = 0; j < max_idx; j++) {
                A2[i][j] -= A2[k][j] * tmp;
			}
		}
	}
}

void Helper::init_w_r(int n, double* &w, double* &r, double** &W, double** &R) {
	int k = 0;
	for (int i = 2; i <= n - 2; i++) {
		for (int j = 2; j <= n - 2; j++) {
			W[i][j] = w[k];
			R[i][j] = r[k];
			k++;
		}
	}
}

void Helper::save_matrix_to_json(int n, double** matrix, string matrix_name) {
	ofstream fout((matrix_name + ".json").c_str());
	StringBuffer sb;
	Writer<StringBuffer> writer(sb);
	writer.StartObject();
	writer.Key(matrix_name.c_str());
	writer.StartArray();
	for (int i = 0; i < n + 1; i++) {
		writer.StartArray();
		for (int j = 0; j < n + 1; j++) {
			writer.Double(matrix[i][j]);
		}
		writer.EndArray();
	}
	writer.EndArray();
	writer.EndObject();

	fout << sb.GetString() << endl;
	fout.close();
}
