#include "SolveTrans.h"
#include <complex>

extern Complex CMath::exp_C(const Complex &rhs);
extern Complex CMath::sin_C(const Complex &rhs);
extern Complex CMath::cos_C(const Complex &rhs);

const double M_Pi = 3.14159265358979323846;

void SolveTrans::QramSolver(int &order, SpMat & matrix, vector<double>& N, const double & time)
{
	vector<Complex> alpha(order/2);
	vector<Complex> theta(order/2);

	Complex theta_1;
	for (int i = 1; i <= order / 2; ++i)
	{
		theta_1 = M_Pi * (double(2 * i - 1) / order);
		theta[i - 1] = order * (0.1309 - 0.1194 * (theta_1 * theta_1) + Complex(0,0.2500) * theta_1);
		alpha[i - 1] = Complex(0, 1) * CMath::exp_C(theta[i - 1]) * (-.1194 * 2 * theta_1 + Complex(0, 0.2500));
	}

	int dim = N.size();
	vector<Complex> N_temp(dim);
	//vector<Complex> x(dim);

	vector < Complex > N_eol;
	N_eol.resize(dim);

	for (int j = 0; j < order / 2; ++j)
	{
		for (int i = 0; i < dim; ++i)
		{
			N_temp[i] = N[i];
		}
		matrix.LUEliminationForCram(theta[j], time, N_temp);
		for (int i = 0; i < dim; ++i)
		{
			N_eol[i] += 2.0 * alpha[j] * N_temp[i];
		}
	}
	for (unsigned int i = 0; i < dim; ++i) {
		N[i] = N_eol[i].real_;
	}
}

void SolveTrans::QramSolver(int &order, SpMat & matrix, vector<double>& N, vector<vector <double> > &F, const double & time)
{
	vector<Complex> alpha(order / 2);
	vector<Complex> theta(order / 2);

	Complex theta_1;
	for (int i = 1; i <= order / 2; ++i)
	{
		theta_1 = M_Pi * (double(2 * i - 1) / order);
		theta[i - 1] = order * (0.1309 - 0.1194 * (theta_1 * theta_1) + Complex(0, 0.2500) * theta_1);
		alpha[i - 1] = Complex(0, 1) * CMath::exp_C(theta[i - 1]) * (-.1194 * 2 * theta_1 + Complex(0, 0.2500));
	}

	int dim = N.size();
	vector<Complex> N_temp(dim);
	//vector<Complex> x(dim);

	vector < Complex > N_eol;
	N_eol.resize(dim);

	for (int j = 0; j < order / 2; ++j)
	{
		//Complex temp = time / theta[j];
		vector<Complex> temp;
		for (int jj = 0; jj < F.size(); ++jj)
		{
			temp.push_back(CMath::Factorial(jj) * pow(time, jj + 1) / CMath::pow_C(theta[j], jj + 1));
		}
		for (int i = 0; i < dim; ++i)
		{
			//N_temp[i] = N[i] + F[i] * temp;
			N_temp[i] = N[i];
			for (int jj = 0; jj < F.size(); ++jj)
			{
				N_temp[i] += F[jj][i] * temp[jj];
			}
		}

		matrix.LUEliminationForCram(theta[j], time, N_temp);
		for (int i = 0; i < dim; ++i)
		{
			N_eol[i] += 2.0 * alpha[j] * N_temp[i];
		}
	}
	for (unsigned int i = 0; i < dim; ++i) {
		N[i] = N_eol[i].real_;
	}
}

void SolveTrans::QramSolver(const int &order, SpMat &matrix, const SpMat &TransMatrixReprocess, SpMat &TransMatrixStockage, vector < double > &N, const double &time) {

	vector<Complex> alpha(order / 2);
	vector<Complex> theta(order / 2);

	Complex theta_1;
	for (int i = 1; i <= order / 2; ++i)
	{
		theta_1 = M_Pi * (double(2 * i - 1) / order);
		theta[i - 1] = order * (0.1309 - 0.1194 * (theta_1 * theta_1) + Complex(0, 0.2500) * theta_1);
		alpha[i - 1] = Complex(0, 1) * CMath::exp_C(theta[i - 1]) * (-.1194 * 2 * theta_1 + Complex(0, 0.2500));
	}

	int dim = N.size();

	// 判断dim是否为奇数，若为奇数，意味着采用了利用增广矩阵方法求解常添料率的方法
	// 若为偶数，则采用数值积分方法求解常添料率燃耗方程
	if (dim % 2 == 0) {
		vector<Complex > N_core(dim / 2);
		//vector <Complex > N_temp(dim / 2);
		vector <Complex > N_stockage(dim / 2);
		vector <Complex > N_eol;
		N_eol.resize(dim);

		//vector <double> N_eol_real;
		//N_eol_real.resize(dim);

		for (int j = 0; j < order / 2; ++j) {
			for (unsigned int i = 0; i < dim / 2; ++i) {
				N_core[i] = N[i];
				N_stockage[i] = N[i + dim / 2];
				//N_temp[i] = N[i] * alpha[j];
			}
			matrix.LUEliminationForCram(theta[j], time, N_core);

			for (unsigned int i = 0; i < dim / 2; ++i) {
				N_eol[i] += 2.0 * alpha[j] * N_core[i];
			}

			for (int ii = 0; ii < dim / 2; ++ii) {
				N_stockage[ii] = N_stockage[ii] - TransMatrixReprocess.diagonal_val_[ii] * time * N_core[ii];
			}
			TransMatrixStockage.LUEliminationForCram(theta[j], time, N_stockage);

			for (unsigned int i = 0; i < dim / 2; ++i) {
				N_eol[i + dim / 2] += 2.0  * alpha[j] * N_stockage[i];
			}

		}
		for (unsigned int i = 0; i < dim; ++i) {
			//N_eol[i] += N[i] * alpha0;
			N[i] = N_eol[i].real_;
		}

		//return N_eol_real;
	}
	else {
		vector<Complex > N_core((dim + 1) / 2); // 堆芯核素浓度向量多一个表示添料率的伪核素，其数值等于1
												//vector <Complex > N_temp(dim / 2);
		vector <Complex > N_stockage((dim - 1) / 2);
		vector <Complex > N_eol;
		N_eol.resize(dim);

		//vector <double> N_eol_real;
		//N_eol_real.resize(dim);


		for (int j = 0; j < order / 2; ++j) {
			for (int i = 0; i < (dim + 1) / 2; ++i) {
				N_core[i] = N[i] * alpha[j];
				if (i < (dim - 1) / 2) {
					N_stockage[i] = N[i + (dim + 1) / 2] * alpha[j];
				}
			}
			matrix.LUEliminationForCram(theta[j], time, N_core);

			for (int i = 0; i < (dim + 1) / 2; ++i) {
				N_eol[i] += 2.0 * N_core[i];
			}

			for (int ii = 0; ii < (dim - 1) / 2; ++ii) {
				N_stockage[ii] = N_stockage[ii] - TransMatrixReprocess.diagonal_val_[ii] * time * N_core[ii];
			}
			TransMatrixStockage.LUEliminationForCram(theta[j], time, N_stockage);

			for (unsigned int i = 0; i < (dim - 1) / 2; ++i) {
				N_eol[i + (dim + 1) / 2] += 2.0 * N_stockage[i];
			}

		}
		for (unsigned int i = 0; i < dim; ++i) {
			//N_eol[i] += N[i] * alpha0;
			N[i] = N_eol[i].real_;
		}

		//return N_eol_real;
	}

}

// 采用拉普拉斯变换方法，求解追踪堆外核素演化的非齐次燃耗方程
void SolveTrans::QramSolver(const int &order, SpMat &matrix, const SpMat &TransMatrixReprocess, SpMat &TransMatrixStockage, vector < double > &N, vector <double> &F, const double &time) {
	vector<Complex> alpha(order / 2);
	vector<Complex> theta(order / 2);

	Complex theta_1;
	for (int i = 1; i <= order / 2; ++i)
	{
		theta_1 = M_Pi * (double(2 * i - 1) / order);
		theta[i - 1] = order * (0.1309 - 0.1194 * (theta_1 * theta_1) + Complex(0, 0.2500) * theta_1);
		alpha[i - 1] = Complex(0, 1) * CMath::exp_C(theta[i - 1]) * (-.1194 * 2 * theta_1 + Complex(0, 0.2500));
	}

	int dim = N.size();

	// 判断dim是否为奇数，若为奇数，意味着采用了利用增广矩阵方法求解常添料率的方法
	// 若为偶数，则采用数值积分方法求解常添料率燃耗方程
	vector<Complex > N_core(dim / 2);
	//vector <Complex > N_temp(dim / 2);
	vector <Complex > N_stockage(dim / 2);
	vector <Complex > N_eol;
	N_eol.resize(dim);

	//vector <double> N_eol_real;
	//N_eol_real.resize(dim);

	for (int j = 0; j < order / 2; ++j) {
		Complex temp = time / theta[j];

		for (unsigned int i = 0; i < dim / 2; ++i) {
			N_core[i] = (N[i] + F[i] * temp) * alpha[j];
			N_stockage[i] = (N[i + dim / 2] + F[i + dim / 2] * temp) * alpha[j];
			//N_temp[i] = N[i] * alpha[j];
		}
		matrix.LUEliminationForCram(theta[j], time, N_core);

		for (unsigned int i = 0; i < dim / 2; ++i) {
			N_eol[i] += 2.0 * N_core[i];
		}

		for (int ii = 0; ii < dim / 2; ++ii) {
			N_stockage[ii] = N_stockage[ii] - TransMatrixReprocess.diagonal_val_[ii] * time * N_core[ii];
		}
		TransMatrixStockage.LUEliminationForCram(theta[j], time, N_stockage);

		for (unsigned int i = 0; i < dim / 2; ++i) {
			N_eol[i + dim / 2] += 2.0 * N_stockage[i];
		}

	}
	for (unsigned int i = 0; i < dim; ++i) {
		//N_eol[i] += N[i] * alpha0;
		N[i] = N_eol[i].real_;
	}

}

