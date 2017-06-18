#pragma once
#include "iostream"

using namespace std;

class Helper {
public:
	Helper(void);
	~Helper(void);
	void set_args(int argc, char* argv[], int &n, double &a, double &b, double &h, double &q0, double &E, double &nu);
	void init_matr_idx(int n, int** &matr_idx);
	void init_A1(int n, double hx, double hy, int** &matr_idx, double** &A1);
	void init_A2(int max_idx, double** &A2);
	void init_Axx_Ayy_Axy(int n, double hx, double hy, int** &matr_idx, double** &Axx, double** &Ayy, double** &Axy);
	void front_process(int max_idx, double** &A1, double** &A2);
	void back_process(int max_idx, double** &A1, double** &A2);
	void init_w_r(int n, double* &w, double* &r, double** &W, double** &R);
	void save_matrix_to_json(int n, double** matrix, string matrix_name);
};
