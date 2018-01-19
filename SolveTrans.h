#ifndef _Solve_Trans_H
#define _Solve_Trans_H

#include "SparseMatrix.h"
#include "NuclList.h"
#include <algorithm>
#include <stack>
//#include "gmp.h"
//#include <boost/multiprecision/cpp_dec_float.hpp>

using namespace std;
//using namespace boost::multiprecision;

class SolveTrans /*燃耗方程的求解方法。各种方法的实现在该类中完成*/
{
public:
	/* CRAM方法求解燃耗方程 */
	/* PFD —— partial fraction decomposition 是CRAM方法的最流行的数值实现算法*/
	int pfd_cram_order = 16; // PFD的CRAM方法的阶数，默认16阶
	vector<double>PfdCramSolver(SpMat &matrix, const vector <double> &N, const double &time);
	vector<double>PfdCramSolver(SpMat &matrix, const SpMat &TransMatrixReprocess, SpMat &TransMatrixStockage, const vector < double > &N, const double &time);
	
	/* IPF —— incomplete partial fractions 是CRAM方法数值实现的一类算法*/
	int ipf_cram_order = 16; // IPF的CRAM方法的阶数，默认16阶
	vector<double>IpfCramSolver(const int &order, SpMat &matrix, const vector <double> &N, const double &time);
	vector<double>IpfCramSolver(const int &order, SpMat &matrix, const SpMat &TransMatrixReprocess, SpMat &TransMatrixStockage, const vector < double > &N, const double &time);

	//--------------------------------------------------------------------------------------//


	/* -----------------------  TTA —— 线性子链法的最基本方法，即Bateman方法 ---------------------------- */
	//mpz_t num;
	long double cutoff_ = 1.0e-20; //截断标志，缺省值为1.0e-20
	long double initial_n_; // 核素初始浓度
	vector< double >initial_n_vector_; //核素初始浓度
	vector< long double > end_n_; // 核素结束浓度

	// ----------------- 标准TTA方法（已遗弃）------------------- //
	//vector<int> node_matirx_index_list_; // 存储节点核素在矩阵中的列号
	//vector<int> node_chain_index_list_;  // 存储节点在链中的序号
	//vector<int> node_daughter_number_list_; // 存储节点的子核个数
	//vector<long double> node_beff_list_; // beff(0) = 1 ; beff(k+1) = beff(k)*lamdaji
	//vector<vector<long double> > node_alpha_list_; // alpha(0) = exp(-lamda_0*t) ; alpha(k+1,i)=alpha(k,i)*1/(lamda_k+1 - lamda_i), alpha(k+1,k+1)
	
	
	//vector<long double> chain_beff_list_;
	//vector<vector<long double> > chain_alpha_list_;
	
	//----------- 采用下述文章中的新的TTA方法的递推公式 -----------//
	//-- “Improvements to the Transmutation Trajectory Analysis of depletion evaluation” ---//
	stack<int> node_a_;
	stack<vector<int> > node_b_;
	stack<vector<vector<long double> > > node_gamma_;

	int chain_a_;
	vector<int> chain_b_;
	vector<vector<long double> > chain_gamma_;
	vector<int> chain_nuclide_id;
	vector<long double> chain_lamda_list_;
	//--------------------------------------------------------//

	vector<int> node_visited_list_; // 燃耗链上的核素访问次数，如果访问一次，则对应位置的元素加一。

	//---------- 连续添料相关参数与函数 ------------------//
	long double tot_feeding_rate_; // 总添料率，用于确定假核素的初始核素浓度
	long double epsilon_ = 1.0e-10; // 用于添料率计算浓度极小量
	vector< long double > feed_rate_; // 各个核素添料率
	vector<int> feed_nuclide_id_;

	//SparseMatrixMCS matrix_;
	vector<vector<int> > matrix_col_index_;
	vector<vector<double> > matrix_col_val_;
	vector<double> matrix_diagonal_val_;

	//void SearchOneChainForDecay(int n0_index, long double time, int& chain_point_count);

	void SearchOneChainForDepletion(const int & n0_index, const long double & time);
	
	void CalOneTree(const int & n0_index, const long double & time) ;
	long double CalCutoffFlag(const int & next_daughter_index, const long double & beff_ji, const long double & time);

	vector<double>TtaSolver(const SparseMatrixMCS &matrix, const double &time);

	vector<double>TtaSolverForFeeding(const SparseMatrixMCS &matrix, const double &time);
	
	void TtaInitialize(const vector <double> &N)
	{
		node_visited_list_.resize(N.size());
		end_n_.resize(N.size());
		initial_n_vector_ = N;

		// ----------- 计算截断判断准则 ----------- //
		long double intial_tot = 0.0;
		for (int i = 0; i < N.size(); ++i)
		{
			if (N[i] != 0.0)
			{
				intial_tot += N[i];
			}
		}
		cutoff_ = cutoff_ * intial_tot;
		// ------------------------------------ //
	}
};

#endif
