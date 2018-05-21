/**
*  @file     SolveTrans.h
*  @brief    燃耗方程核心求解模块
*  @author   夏少鹏
*  @email    xiashaopeng@sinap.ac.cn
*  @version  1.0.0
*  @date     01/30/2018
*/

#ifndef _Solve_Trans_H
#define _Solve_Trans_H

#include "SparseMatrix.h"
#include "NuclList.h"
#include <algorithm>
#include <stack>

using namespace std;

/**
 * @brief MODEC的燃耗方程核心求解类
 * 包含了CRAM算法和TTA算法两大类燃耗算法
*/
class SolveTrans {
  public:
    /* CRAM方法求解燃耗方程 */
    /* PFD —— partial fraction decomposition 是CRAM方法的最流行的数值实现算法*/
    int pfd_cram_order = 16;	///< PFD格式的CRAM方法的阶数，默认16阶

    /**
     * @brief PFD格式的CRAM求解器
     * @param[in] matrix The matrix storing the whole coefficients of transmutation
     * @param[in] N The nuclide concentration vector at beginning of the sub-step
     * @param[in] time The sub-step time of burnup calculations
     * @return The nuclide concentration vector at the end of the sub-step
    */
    void PfdCramSolver(SpMat &matrix, vector <double> &N, const double &time);

    /**
     * @brief PFD格式的CRAM求解器重载，用于追踪堆外核素演化
     * @param[in] matrix The matrix storing the whole coefficients of transmutation
    * @param[in] TransMatrixReprocess The matrix storing the coefficients of decay in reprocess unit
           * @param[in] N The nuclide concentration vector at beginning of the sub-step
           * @param[in] time The sub-step time of burnup calculations
           * @return The nuclide concentration vector at the end of the sub-step
          */
    void PfdCramSolver(SpMat &matrix, const SpMat &TransMatrixReprocess, SpMat &TransMatrixStockage, vector < double > &N, const double &time);

    /* IPF —— incomplete partial fractions 是CRAM方法数值实现的一类算法*/
    int ipf_cram_order = 16;	///< IPF格式的CRAM方法的阶数，默认16阶

    /**
    * @brief IPF格式的CRAM求解器
    */
    void IpfCramSolver32(SpMat &matrix, vector <double> &N, const double &time);
	
	/**
    * @brief IPF格式的CRAM求解器
    */
    void IpfCramSolver48(SpMat &matrix, vector <double> &N, const double &time);

    /**
    * @brief IPF格式的求解器重载，用于追踪堆外核素演化
    */
    void IpfCramSolver32(SpMat &matrix, const SpMat &TransMatrixReprocess, SpMat &TransMatrixStockage, vector < double > &N, const double &time);

    /**
    * @name TTA公式主要计算参数
    * @{
    */
    long double cutoff_std_ = 1.0e-20;			///< TTA标准截断浓度
    long double cutoff_ ;					///< TTA方法实际截断浓度: cutoff_ = cutoff_std * (总的初始核素浓度)
    long double initial_n_; 				///< 总的核素初始浓度
    vector< double >initial_n_vector_;			///< 各核素初始浓度
    vector< long double > end_n_;				///< 各核素燃耗步末浓度

    stack<int> node_a_;					///< 递推公式中系数a的存储向量(只在分叉点处存储)
    stack<vector<int> > node_b_;				///< 递推公式中系数b的存储向量(只在分叉点处存储)
    stack<vector<vector<long double> > > node_gamma_;	///< 递推公式中系数gamma的存储向量(只在分叉点处存储)

    int chain_a_;						///< 当前待求核素的系数a
    vector<int> chain_b_;					///< 当前待求核素的系数b
    vector<vector<long double> > chain_gamma_;		///< 当前待求核素的系数gamma
    vector<int> chain_nuclide_id;				///< 当前求解的燃耗链上各核素id
    vector<long double> chain_lamda_list_;			///< 当前求解的燃耗链上各核素的等效衰变系数

    vector<int> node_visited_list_; 			///< 燃耗链上的核素访问次数，如果访问一次，则对应位置的元素加一。
    /** @} */

    /**
    * @name TTA在线添料模块相关计算参数
    * @{
    */
    long double tot_feeding_rate_; 				///< 总添料率，用于确定假核素的初始核素浓度
    long double epsilon_ = 1.0e-10; 			///< 用于维持添料率恒定而人为设定的一个接近于0的浮点数
    vector< long double > feed_rate_; 			///< 各个添料率核素的添料率
    vector<int> feed_nuclide_id_;				///< 各个添料率核素的id

    vector<vector<int> > matrix_col_index_;			///< 邻接燃耗矩阵行索引
    vector<vector<double> > matrix_col_val_;		///< 对应行索引的各个系数
    vector<double> matrix_diagonal_val_;			///< 邻接燃耗矩阵对角元素
    /** @} */

    /**
    * @brief 递归搜索燃耗链
    */
    void SearchOneChainForDepletion(const int & n0_index, const long double & time);

    /**
    * @brief 计算由一个顶端核素延伸的完整燃耗链
    */
    void CalOneTree(const int & n0_index, const long double & time) ;

    /**
    * @brief 截断条件判断
    */
    long double CalCutoffFlag(const int & next_daughter_index, const long double & beff_ji, const long double & time);

    /**
     * @brief 齐次燃耗方程的TTA求解器
     * @param[in] matrix	The adjacency matrix for TTA solver
     * @param[in] time	Burnup time for TTA solving
     * @return	Return the nuclide concentration vector at the end of the burnup time
    */
    vector<double>TtaSolver(const SparseMatrixMCS &matrix, const double &time);

    /**
     * @brief 齐次燃耗方程的TTA求解器重载
     * @param[in] matrix    The adjacency matrix for TTA solver
     * @param[in] initial_n_vector	The nuclide concentration vector at the beginning of the burnup time
     * @param[in] time      Burnup time for TTA solving
     */
    void TtaSolver(const SparseMatrixMCS &matrix, vector<double> &initial_n_vector, const double &time);

    /**
     * @brief 非齐次燃耗方程的TTA求解器
     * @param[in] matrix    The adjacency matrix for TTA solver
     * @param[in] time      Burnup time for TTA solving
     * @return      Return the nuclide concentration vector at the end of the burnup time
     */
    vector<double>TtaSolverForFeeding(const SparseMatrixMCS &matrix, const double &time);

    /**
     * @brief 非齐次燃耗方程的TTA求解器重载
     * @param[in] matrix    The adjacency matrix for TTA solver
     * @param[in] initial_n_vector The nuclide concentration vector at the beginning of the burnup time
     * @param[in] time      Burnup time for TTA solving
     */
    void TtaSolverForFeeding(const SparseMatrixMCS &matrix, vector<double> &initial_n_vector, const double &time);

    /**
     * @brief TTA求解器初始化
     * @param[in] N	The nuclide concentration vector at the beginning of the burnup time
     */
    void TtaInitialize(const vector <double> &N) {
        node_visited_list_.resize(N.size());
        end_n_.resize(N.size());
        initial_n_vector_ = N;

        // ----------- 计算截断判断准则 ----------- //
        long double intial_tot = 0.0;
        for (int i = 0; i < N.size(); ++i) {
            if (N[i] != 0.0) {
                intial_tot += N[i];
            }
        }
        cutoff_ = cutoff_std_ * intial_tot;
        // ------------------------------------ //
    }

    /**
     * @brief TTA求解器初始化重载
     * @param[in] size The nuclide concentration vector at the beginning of the burnup time
     */
    void TtaInitialize(int size) {
        node_visited_list_.resize(size);
        end_n_.resize(size);
    }

    /** @} */
};

#endif
