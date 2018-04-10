#include "SolveTrans.h"


// TTA用于衰变计算，没有闭合环路，基本不存在完全相等的衰变常数，故采用递推关系式求解系数最为高效 //
//void SolveTrans::SearchOneChainForDecay(int n0_index, long double time, int& chain_point_count)
//{
//	long double beff = chain_beff_list_.back();
//	vector<long double> alpha = chain_alpha_list_.back();
//
//	vector<int> next_index_list = matrix_col_index_[n0_index];
//
//	if (next_index_list[0] == -1)
//	{
//		return;
//	}
//	if (next_index_list.size() == 1)
//	{
//		int next_daughter_index = next_index_list[0];
//		long double lamda_ii = long double(chain_lamda_list_.back());
//		long double lamda_jj = -matrix_diagonal_val_[next_daughter_index].real();
//		long double beff_ji = matrix_col_val_[n0_index][0].real();
//
//		beff = beff * beff_ji / lamda_ii; // 按照标准的公式来计算TTA的相关系数
//
//		long double alpha_new = 1;
//		long double alpha_sum = 0;
//		for (int ii = 0; ii < alpha.size(); ++ii)
//		{
//			if (lamda_jj == chain_lamda_list_[ii]) // 这里的lamda系数的处理可能存在问题
//			{
//				lamda_jj = lamda_jj*(1 + 1e-5);
//			}
//
//			alpha[ii] = alpha[ii] * long double(lamda_ii / (lamda_jj - chain_lamda_list_[ii]));
//
//			alpha_sum += alpha[ii];
//
//			alpha_new = alpha_new * long double(chain_lamda_list_[ii] / (chain_lamda_list_[ii] - lamda_jj));
//
//		}
//		alpha_sum += alpha_new * exp(-lamda_jj*time);
//		alpha.push_back(alpha_new * exp(-lamda_jj*time));
//
//		end_n_[next_daughter_index] += initial_n_ * alpha_sum * beff;
//
//		chain_beff_list_.push_back(beff);
//		chain_alpha_list_.push_back(alpha);
//		chain_lamda_list_.push_back(lamda_jj);
//
//		chain_point_count++;
//
//		SearchOneChainForDecay(next_daughter_index, time, chain_point_count);
//	}
//	if (next_index_list.size() > 1)
//	{
//		node_matirx_index_list_.push_back(n0_index);
//		node_chain_index_list_.push_back(chain_point_count);
//		node_alpha_list_.push_back(alpha); // 第一个核素作为首个节点，其alpha系数需要存储
//		node_beff_list_.push_back(beff); // 第一个核素作为首个节点，其beff系数需要存储
//		node_daughter_number_list_.push_back(next_index_list.size());
//
//		int node_location = node_chain_index_list_.back(); //
//
//		for (int i = 0; i < next_index_list.size(); ++i)
//		{
//			if (i > 0)
//			{
//				for (int j = node_location; j < chain_point_count; ++j)
//				{
//					chain_lamda_list_.pop_back();
//					chain_beff_list_.pop_back();
//					chain_alpha_list_.pop_back();
//				}
//				beff = chain_beff_list_.back();
//				alpha = chain_alpha_list_.back();
//
//				chain_point_count = node_location;
//			}
//			int next_daughter_index = next_index_list[i];
//
//			long double lamda_ii = long double(chain_lamda_list_.back());
//			long double lamda_jj = -matrix_diagonal_val_[next_daughter_index].real();
//			long double beff_ji = matrix_col_val_[n0_index][i].real();
//
//			beff = beff * beff_ji / lamda_ii;
//
//			long double alpha_new = 1;
//			long double alpha_sum = 0;
//
//			for (int ii = 0; ii < alpha.size(); ++ii)
//			{
//				if (lamda_jj == chain_lamda_list_[ii]) // 这里的lamda系数的处理可能存在问题
//				{
//					lamda_jj = lamda_jj*(1 + 1e-5);
//				}
//
//				alpha[ii] = alpha[ii] * long double(lamda_ii / (lamda_jj - chain_lamda_list_[ii]));
//
//				alpha_sum += alpha[ii];
//
//				alpha_new = long double(alpha_new * chain_lamda_list_[ii] / (chain_lamda_list_[ii] - lamda_jj));
//			}
//			alpha_sum += alpha_new * exp(-lamda_jj*time);
//			alpha.push_back(alpha_new * exp(-lamda_jj*time));
//
//			end_n_[next_daughter_index] += initial_n_ * alpha_sum * beff;
//
//			chain_beff_list_.push_back(beff);
//			chain_alpha_list_.push_back(alpha);
//			chain_lamda_list_.push_back(lamda_jj);
//
//			chain_point_count++;
//
//			SearchOneChainForDecay(next_daughter_index, time, chain_point_count);
//		}
//		return;
//	}
//
//}

// 对于燃耗问题，不考虑保留每步的alpha和beff（内存爆炸），而是只保留每个核素的衰变和截面信息，再每步时计算alpha和beff
//void SolveTrans::SearchOneChainForDepletion(SparseMatrixMCS &matrix_, int n0_index, long double time, int& chain_point_count)
//{
//	vector<int> next_index_list = matrix_.RowListIndex(n0_index);
//	long double lamda_ii = -matrix_.diagonal_val_[n0_index].real();
//	if (next_index_list[0] == -1)
//	{
//		return;
//	}
//	if (next_index_list.size() == 1)
//	{
//		int next_daughter_index = next_index_list[0];
//		node_visited_list_[next_daughter_index] += 1;
//		long double n0_i1 = 0;
//		if (node_visited_list_[next_daughter_index] == 1)
//		{
//			n0_i1 = initial_n_vector_[next_daughter_index].real();
//		}
//		else
//		{
//			n0_i1 = 0;
//		}
//
//		long double lamda_jj = -matrix_.diagonal_val_[next_daughter_index].real();
//		long double beff_ji = matrix_.ElementCCS(next_daughter_index, n0_index).real();
//
//		//beff = beff * beff_ji / lamda_ii; // 按照标准的公式来计算TTA的相关系数
//
//		int visited_number = 0;
//		for (int ii = 0; ii < chain_nuclide_id.size(); ++ii)
//		{
//			if (next_daughter_index == chain_nuclide_id[ii]) // 这里的lamda系数的处理可能存在问题
//			{
//				visited_number ++;
//			}
//		}
//
//		if (visited_number > 0)
//		{
//			long double judge_cutoff = CalCutoffFlag(next_daughter_index, lamda_ii, time);
//			n0_i1 = 0;
//			if (abs(judge_cutoff) <= cutoff_)
//			{
//				return;
//			}
//		}
//
//		long double temp_n = 0; // 子核的核素浓度
//		long double temp_gamma_j0 = 0;
//		long double temp_gamma_i1 = 0;
//		int equal_j = -1;
//
//		/*int chain_a = chain_a_;
//		vector<int> chain_b = chain_b_;*/
//		vector<vector<long double> > chain_gamma = chain_gamma_;
//
//		for (int j = 0; j < chain_a_; ++j)
//		{
//			if (lamda_jj == chain_lamda_list_[j])
//			{
//				equal_j = j;
//			}
//			else
//			{
//				for (int k = chain_b_[j]; k >= 0; --k)
//				{
//					if (k == chain_b_[j])
//					{
//						chain_gamma_[j][k] = chain_gamma_[j][k] * beff_ji / (lamda_jj - chain_lamda_list_[j]);
//					}
//					else
//					{
//						chain_gamma_[j][k] = (chain_gamma_[j][k] * beff_ji - (k + 1)*chain_gamma_[j][k + 1])
//							/ (lamda_jj - chain_lamda_list_[j]);
//					}
//					if (k == 0)
//					{
//						temp_gamma_j0 += chain_gamma_[j][k];
//					}
//					temp_n += chain_gamma_[j][k] * pow(time, k)*exp(-chain_lamda_list_[j] * time);
//				}
//			}
//		}
//
//		temp_gamma_i1 = n0_i1 - temp_gamma_j0;
//
//		if (equal_j == -1)
//		{
//			chain_a_ += 1;
//			chain_b_.push_back(0);
//			chain_gamma_.push_back({ n0_i1 - temp_gamma_j0 });
//			temp_n += chain_gamma_.back()[0] * exp(-lamda_jj * time);
//		}
//		else
//		{
//			chain_b_[equal_j] += 1;
//			for (int k = 0; k <= chain_b_[equal_j]; ++k)
//			{
//				if (k == chain_b_[equal_j])
//				{
//					chain_gamma_[equal_j].push_back(beff_ji*chain_gamma[equal_j][k - 1] / k);
//					temp_n += chain_gamma_[equal_j].back() * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
//				}
//				else if (k > 0 && k < chain_b_[equal_j])
//				{
//					chain_gamma_[equal_j][k] = beff_ji*chain_gamma[equal_j][k - 1] / k;
//					temp_n += chain_gamma_[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
//				}
//				else
//				{
//					chain_gamma_[equal_j][k] = n0_i1 - temp_gamma_j0;
//					temp_n += chain_gamma_[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
//				}
//
//			}
//		}
//
//		end_n_[next_daughter_index] += temp_n;
//
//		if (equal_j == -1)
//		{
//			chain_lamda_list_.push_back(lamda_jj);
//			chain_nuclide_id.push_back(next_daughter_index);
//		}
//
//		chain_point_count++;
//
//		SearchOneChainForDepletion(matrix_, next_daughter_index, time, chain_point_count);
//	}
//	if (next_index_list.size() > 1)
//	{
//		node_a_.push_back(chain_a_);
//		node_b_.push_back(chain_b_);
//		//node_chain_index_list_.push_back(chain_point_count);
//		//node_alpha_list_.push_back(alpha); // 第一个核素作为首个节点，其alpha系数需要存储
//		//node_beff_list_.push_back(beff); // 第一个核素作为首个节点，其beff系数需要存储
//		node_gamma_.push_back(chain_gamma_);
//
//		//int node_location = node_chain_index_list_.back(); //
//
//		for (int i = 0; i < next_index_list.size(); ++i)
//		{
//			if (i > 0)
//			{
//				chain_a_ = node_a_.back();
//				chain_b_ = node_b_.back();
//				chain_gamma_ = node_gamma_.back();
//
//				int size = chain_lamda_list_.size();
//				for (int j = chain_a_; j < size; ++j)
//				{
//					chain_lamda_list_.pop_back();
//					chain_nuclide_id.pop_back();
//				}
//
//				//chain_point_count = node_location;
//			}
//
//			int next_daughter_index = next_index_list[i];
//			node_visited_list_[next_daughter_index] += 1;
//			long double n0_i1 = 0;
//			if (node_visited_list_[next_daughter_index] == 1)
//			{
//				n0_i1 = initial_n_vector_[next_daughter_index].real();
//			}
//			else
//			{
//				n0_i1 = 0;
//			}
//
//			long double lamda_jj = -matrix_.diagonal_val_[next_daughter_index].real();
//			long double beff_ji = matrix_.ElementCCS(next_daughter_index, n0_index).real();
//
//			//beff = beff * beff_ji / lamda_ii; // 按照标准的公式来计算TTA的相关系数
//
//			int visited_number = 0;
//			for (int ii = 0; ii < chain_nuclide_id.size(); ++ii)
//			{
//				if (next_daughter_index == chain_nuclide_id[ii]) // 这里的lamda系数的处理可能存在问题
//				{
//					visited_number++;
//				}
//			}
//
//			if (visited_number > 0)
//			{
//				long double judge_cutoff = CalCutoffFlag(next_daughter_index, lamda_ii, time);
//				n0_i1 = 0;
//				if (abs(judge_cutoff) <= cutoff_)
//				{
//					node_a_.pop_back();
//					node_b_.pop_back();
//					node_gamma_.pop_back();
//
//					return;
//				}
//			}
//
//			long double temp_n = 0; // 子核的核素浓度
//			long double temp_gamma_j0 = 0;
//			int equal_j = -1;
//
//			vector<vector<long double> > chain_gamma = chain_gamma_;
//
//			for (int j = 0; j < chain_a_; ++j)
//			{
//				if (lamda_jj == chain_lamda_list_[j])
//				{
//					equal_j = j;
//				}
//				else
//				{
//					for (int k = chain_b_[j]; k >= 0; --k)
//					{
//						if (k == chain_b_[j])
//						{
//							chain_gamma_[j][k] = chain_gamma[j][k] * beff_ji / (lamda_jj - chain_lamda_list_[j]);
//						}
//						else
//						{
//							chain_gamma_[j][k] = (chain_gamma[j][k] * beff_ji - (k + 1)*chain_gamma_[j][k + 1])
//								/ (lamda_jj - chain_lamda_list_[j]);
//						}
//						if (k == 0)
//						{
//							temp_gamma_j0 += chain_gamma_[j][k];
//						}
//						temp_n += chain_gamma_[j][k] * pow(time, k)*exp(-chain_lamda_list_[j] * time);
//					}
//				}
//			}
//
//			if (equal_j == -1)
//			{
//				chain_a_ += 1;
//				chain_b_.push_back(0);
//				chain_gamma_.push_back({ n0_i1 - temp_gamma_j0 });
//				chain_lamda_list_.push_back(lamda_jj);
//
//				chain_nuclide_id.push_back(next_daughter_index);
//				temp_n += chain_gamma_[chain_a_ - 1][0] * exp(-lamda_jj * time);
//			}
//			else
//			{
//				chain_b_[equal_j] += 1;
//				for (int k = 0; k <= chain_b_[equal_j]; ++k)
//				{
//					if (k == chain_b_[equal_j])
//					{
//						chain_gamma_[equal_j].push_back(beff_ji*chain_gamma[equal_j][k - 1] / k);
//						temp_n += chain_gamma_[equal_j].back() * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
//					}
//					else if (k > 0 && k < chain_b_[equal_j])
//					{
//						chain_gamma_[equal_j][k] = beff_ji*chain_gamma[equal_j][k - 1] / k;
//						temp_n += chain_gamma_[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
//					}
//					else
//					{
//						chain_gamma_[equal_j][k] = n0_i1 - temp_gamma_j0;
//						temp_n += chain_gamma_[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
//					}
//				}
//			}
//
//			end_n_[next_daughter_index] += temp_n;
//
//			chain_point_count++;
//
//			SearchOneChainForDepletion(matrix_, next_daughter_index, time, chain_point_count);
//		}
//
//		node_a_.pop_back();
//		node_b_.pop_back();
//		node_gamma_.pop_back();
//
//		return;
//	}
//
//}

// 采用人为构造堆栈的方式，使用循环替代递归，目的是提高计算效率（是否能提高还需检验）
// 将堆栈定义为局部变量，这样不需要人为释放堆栈空间
// 将存储TTA系数的向量仍然设为类成员，原因是方便其他函数调用，主要是截断判断函数，因此其需要手动清除
void SolveTrans::SearchOneChainForDepletion(const int & n0_index, const long double & time) {
    
    // 重新定义局部变量
    // int chain_a_;
    // vector<int> chain_b_;
    // vector<vector<long double> > chain_gamma_;
    // vector<long double> chain_lamda_list_;
    
    stack<int> node_a_;					///< 递推公式中系数a的存储向量(只在分叉点处存储)
    stack<vector<int> > node_b_;				///< 递推公式中系数b的存储向量(只在分叉点处存储)
    stack<vector<vector<long double> > > node_gamma_;	///< 递推公式中系数gamma的存储向量(只在分叉点处存储)
    stack<vector<long double> > node_lamda_list_;

    // 定义一个元素为vector的栈，表示整个燃耗链
    // 其中，vector第一个元素对应母核素，之后的元素对应各个子核素
    stack <vector<int> > nuclide_chain;
    // 对应nuclide_chain栈，定义一个用于存储母核到各个子核素的堆栈
    stack <vector<double> > nuclide_to_daughter;

    // 定义一个哨兵数组向量，用于描述该核素是否在栈中
    // 如果出现重复的核素，该位置的元素需设为true
    // vector<bool> check_repeat_nuclide;
    // check_repeat_nuclide.resize(matrix_diagonal_val_.size());
    
    // 定义一个int型数组，用来存储母核素及其子核素id，并且母核素位于数组的第一个位置。
    vector<int> nuclide_group (matrix_col_index_[n0_index]);
    nuclide_group.insert(nuclide_group.begin(), n0_index);

    nuclide_chain.push(nuclide_group);

    // 初始化堆栈
    nuclide_to_daughter.push(matrix_col_val_[n0_index]);

    // 初始化各个系数
    chain_a_ = 1;
    chain_b_.push_back(0);
    chain_gamma_.push_back({ initial_n_ });
    chain_lamda_list_.push_back(-matrix_diagonal_val_[n0_index]);

    node_a_.push(chain_a_);
    node_b_.push(chain_b_);
    node_gamma_.push(chain_gamma_);
    node_lamda_list_.push(chain_lamda_list_);

    // 根据根结点系数计算根节点核素浓度，并将哨兵数组该核素位置加一，表示访问过一次了
    end_n_[n0_index] += initial_n_ * exp(-matrix_diagonal_val_[n0_index]*time);
    node_visited_list_[n0_index] += 1;
    //check_repeat_nuclide[n0_index] = true;

    // 如果该核素之前访问过，则意味着该核素的所有子核素都已经计算过一遍，直接返回
    if (node_visited_list_[n0_index] > 1) {
        return;
    }

    while (!nuclide_chain.empty()) {
        int mother_nuclide (nuclide_chain.top()[0]);

        // 将栈顶数组中的最后一个元素作为下一个循环的母核素，同时在栈顶数组中删除该核素
        // 如果数组中只剩下一个元素，意味着母核素的所有子核素都已经找到，
        // 这时候需要从堆栈中弹出，并且直接进入下一个循环
        if (nuclide_chain.top().size() == 1) {
            nuclide_chain.pop();
            nuclide_to_daughter.pop();
            node_a_.pop();
            node_b_.pop();
            node_gamma_.pop();
            node_lamda_list_.pop();

            chain_a_ = node_a_.top();
            chain_b_ = node_b_.top();
            chain_gamma_ = node_gamma_.top();
            chain_lamda_list_ = node_lamda_list_.top();

            continue;
        }
        int next_daughter_nuclide (nuclide_chain.top().back());
        nuclide_chain.top().pop_back();
        
        node_visited_list_[next_daughter_nuclide] += 1;
        // ------ TTA实现程序主体 ------- //
        long double lamda_jj(-matrix_diagonal_val_[next_daughter_nuclide]);

        long double beff_ji(nuclide_to_daughter.top().back());
        nuclide_to_daughter.top().pop_back();

        long double n0_i1;
        long double temp_n(0.0); // 子核的核素浓度
        long double temp_gamma_j0(0.0);
        int equal_j(-1);

        // 判断该核素是否曾经计算过，如果计算过，则该核素的初始浓度需要设为0
        if (node_visited_list_[next_daughter_nuclide] == 1) {
            n0_i1 = initial_n_vector_[next_daughter_nuclide];
        } else {
            n0_i1 = 0;
        }
        // 截断判断
        if (nuclide_chain.size() >= 5) {
            long double judge_cutoff(CalCutoffFlag(next_daughter_nuclide, -matrix_diagonal_val_[mother_nuclide], time));
            if (abs(judge_cutoff) <= cutoff_) {
                continue;
            }
        }

        for (int j = 0; j < chain_a_; ++j) {
            if (lamda_jj == chain_lamda_list_[j]) {
                equal_j = j;
            } else {
                for (int k = chain_b_[j]; k >= 0; --k) {
                    if (k == chain_b_[j]) {
                        chain_gamma_[j][k] = chain_gamma_[j][k] * beff_ji / (lamda_jj - chain_lamda_list_[j]);
                    } else {
                        chain_gamma_[j][k] = (chain_gamma_[j][k] * beff_ji - (k + 1)*chain_gamma_[j][k + 1])
                                             / (lamda_jj - chain_lamda_list_[j]);
                    }
                    if (k == 0) {
                        temp_gamma_j0 += chain_gamma_[j][k];
                    }
                    temp_n += chain_gamma_[j][k] * pow(time, k)*exp(-chain_lamda_list_[j] * time);
                }
            }
        }

        if (equal_j == -1) {
            chain_a_ += 1;
            chain_b_.push_back(0);
            chain_gamma_.push_back({  n0_i1 - temp_gamma_j0 });
            chain_lamda_list_.push_back(lamda_jj);

            //chain_nuclide_id.push_back(next_daughter_nuclide);
            temp_n += chain_gamma_[chain_a_ - 1][0] * exp(-lamda_jj * time);
        } else {
            vector<long double> chain_gamma(chain_gamma_[equal_j]);

            chain_b_[equal_j] += 1;
            for (int k = 0; k <= chain_b_[equal_j]; ++k) {
                if (k == chain_b_[equal_j]) {
                    chain_gamma_[equal_j].push_back(beff_ji*chain_gamma[k - 1] / k);
                    temp_n += chain_gamma_[equal_j].back() * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
                } else if (k > 0 && k < chain_b_[equal_j]) {
                    chain_gamma_[equal_j][k] = beff_ji*chain_gamma[k - 1] / k;
                    temp_n += chain_gamma_[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
                } else {
                    chain_gamma_[equal_j][k] = n0_i1 - temp_gamma_j0;
                    temp_n += chain_gamma_[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
                }
            }
        }
        if (temp_n >= 0) {
            end_n_[next_daughter_nuclide] += temp_n;
        }
        // 将next_mother_nuclide核素对应的nuclide_group压入堆栈
        vector<int> nuclide_group (matrix_col_index_[next_daughter_nuclide]);
        nuclide_group.insert(nuclide_group.begin(), next_daughter_nuclide);
        nuclide_chain.push(nuclide_group);

        nuclide_to_daughter.push(matrix_col_val_[next_daughter_nuclide]);

        node_a_.push(chain_a_);
        node_b_.push(chain_b_);
        node_gamma_.push(chain_gamma_);
        node_lamda_list_.push(chain_lamda_list_);
    }    
}

long double SolveTrans::CalCutoffFlag(const int & next_daughter_index, const long double & beff_ji, const long double & time) {
    long double lamda_jj = 0.0;// -matrix_diagonal_val_[next_daughter_index];
    //long double alpha_new = 1;
    //long double alpha_sum = 0;
    long double n0_i1;
    if (node_visited_list_[next_daughter_index] == 1) {
        n0_i1 = initial_n_vector_[next_daughter_index];
    } else {
        n0_i1 = 0;
    }
    int chain_a = chain_a_;
    vector<int> chain_b = chain_b_;
    vector<vector<long double> > chain_gamma = chain_gamma_;

    long double temp_n = 0; // 子核的核素浓度
    long double temp_gamma_j0 = 0;
    int equal_j = -1;
    for (int j = 0; j < chain_a; ++j) {
        if (lamda_jj == chain_lamda_list_[j]) {
            equal_j = j;
        } else {
            for (int k = chain_b[j]; k >= 0; --k) {
                if (k == chain_b[j]) {
                    chain_gamma[j][k] = chain_gamma[j][k] * beff_ji / (lamda_jj - chain_lamda_list_[j]);
                } else {
                    chain_gamma[j][k] = (chain_gamma[j][k] * beff_ji - (k + 1)*chain_gamma[j][k + 1])
                                        / (lamda_jj - chain_lamda_list_[j]);
                }
                if (k == 0) {
                    temp_gamma_j0 += chain_gamma[j][k];
                }
                temp_n += chain_gamma[j][k] * pow(time, k)*exp(-chain_lamda_list_[j] * time);
            }
        }
    }

    if (equal_j == -1) {
        temp_n += (n0_i1 - temp_gamma_j0) * exp(-lamda_jj * time);
    } else {
        chain_b[equal_j] += 1;
        for (int k = 0; k <= chain_b[equal_j]; ++k) {
            if (k == chain_b[equal_j]) {
                chain_gamma[equal_j].push_back(beff_ji*chain_gamma_[equal_j][k - 1] / k);
                temp_n += chain_gamma[equal_j].back() * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
            } else if (k > 0 && k < chain_b[equal_j]) {
                chain_gamma[equal_j][k] = beff_ji*chain_gamma_[equal_j][k - 1] / k;
                temp_n += chain_gamma[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
            } else {
                chain_gamma[equal_j][k] = n0_i1 - temp_gamma_j0;
                temp_n += chain_gamma[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
            }
        }
    }
    if(temp_n > 0.0) {
        return temp_n;
    } else {
        return 0.0;
    }
}

void SolveTrans::CalOneTree(const int &n0_index, const long double &time) {

    // ----------------- search one chain depletion using recursion method  ------------------ //
    SearchOneChainForDepletion(n0_index, time); //

    //------------------------------- 完成一条燃耗树的计算，清理相关向量 ------------------------------//
    chain_lamda_list_.clear();
    chain_b_.clear();
    chain_gamma_.clear();
}

vector<double> SolveTrans::TtaSolver(const SparseMatrixMCS &matrix,const double &time ) {
    // for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix)
    // end_n_[ix] = 0.0;

    matrix_col_index_ = matrix.col_index_;
    matrix_col_val_ = matrix.col_val_;
    matrix_diagonal_val_ = matrix.diagonal_val_;

    for (int i = 0; i < initial_n_vector_.size(); ++i) {
        if (initial_n_vector_[i] != 0.0) {
            initial_n_ = initial_n_vector_[i];
            CalOneTree(i, time);
        }
    }
    for (vector<int>::size_type ix = 0; ix != node_visited_list_.size(); ++ix)
        node_visited_list_[ix] = 0;

    vector<double> end_out;
    end_out.resize(end_n_.size());
    for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix) {
        end_out[ix] = end_n_[ix];
        end_n_[ix] = 0.0;
    }

    return end_out;
}

void SolveTrans::TtaSolver(const SparseMatrixMCS &matrix, vector<double> &initial_n_vector, const double &time) {
    long double initial_tot = 0.0;
    for (int i = 0; i < initial_n_vector.size() ; ++i) {
        initial_tot += initial_n_vector[i];
    }
    cutoff_ = cutoff_std_ * initial_tot;

    initial_n_vector_.resize(initial_n_vector.size());
    initial_n_vector_ = initial_n_vector;

    // initialize the sort vector
    vector<int> sort_vector;
    sort_vector.resize(initial_n_vector.size());
    for (vector<int>::size_type ix = 0; ix != sort_vector.size(); ++ix)
        sort_vector[ix] = ix;


    matrix_col_index_ = matrix.col_index_;
    matrix_col_val_ = matrix.col_val_;
    matrix_diagonal_val_ = matrix.diagonal_val_;

    // employ sort algrothim to sort the nuclide concentrations and store the sorting in sort_vector
    double key;
    int i, key_s;
    for(int j = 1; j < initial_n_vector.size(); ++j) {
        key = initial_n_vector[j];
        key_s = sort_vector[j];
        i = j - 1;
        while( i > -1 && initial_n_vector[i] < key) {
            initial_n_vector[i + 1] = initial_n_vector[i];
            sort_vector[i + 1] = sort_vector[i];
            i -= 1;
        }
        initial_n_vector[i + 1] = key;
        sort_vector[i + 1] = key_s;
    }

    int index;
    for (int i = 0; i < initial_n_vector_.size(); ++i) {
        index = sort_vector[i];
        if ( initial_n_vector_[index] != 0.0 ) {
            initial_n_ = initial_n_vector_[index];
            CalOneTree(index, time);
        }
    }

    for (vector<int>::size_type ix = 0; ix != node_visited_list_.size(); ++ix)
        node_visited_list_[ix] = 0;

    for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix) {
        initial_n_vector[ix] = end_n_[ix];
        end_n_[ix] = 0.0;
    }
}

vector<double> SolveTrans::TtaSolverForFeeding(const SparseMatrixMCS &matrix, const double &time) {
    // for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix)
    // end_n_[ix] = 0.0;

    matrix_col_index_ = matrix.col_index_;
    matrix_col_val_ = matrix.col_val_;
    matrix_diagonal_val_ = matrix.diagonal_val_;

    matrix_diagonal_val_.back() = - epsilon_ / time;

    long double dummy_initial_nuclide(tot_feeding_rate_ * time / epsilon_);
    initial_n_vector_[initial_n_vector_.size() - 1] = dummy_initial_nuclide; // 伪核素的初始值


    for (unsigned int i = 0; i < feed_rate_.size(); ++i) {
        //int nuclide_id = feed_nuclide_id_[i];
        matrix_col_index_.back().push_back(feed_nuclide_id_[i]);
        matrix_col_val_.back().push_back(feed_rate_[i] / dummy_initial_nuclide);
    }

    int size = initial_n_vector_.size();
    for (int i = 0; i < size ; ++i) {
        if (initial_n_vector_[i] != 0.0) {
            initial_n_ = initial_n_vector_[i];
            CalOneTree(i, time);
        }
    }
    for (vector<int>::size_type ix = 0; ix != node_visited_list_.size(); ++ix)
        node_visited_list_[ix] = 0;

    vector<double> end_out;
    end_out.resize(end_n_.size());
    for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix) {
        end_out[ix] = end_n_[ix];
        end_n_[ix] = 0.0;
    }

    return end_out;
}

void SolveTrans::TtaSolverForFeeding(const SparseMatrixMCS &matrix, vector<double> &initial_n_vector, const double &time) {
    // for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix)
    // end_n_[ix] = 0.0;

    //node_visited_list_.resize(initial_n_vector.size());
    //end_n_.resize(initial_n_vector.size());

    long double initial_tot = 0.0;
    for (int i = 0; i < initial_n_vector.size() - 1 ; ++i) {
        initial_tot += initial_n_vector[i];
    }
    cutoff_ = cutoff_std_ * initial_tot;

    initial_n_vector_ = initial_n_vector;

    matrix_col_index_ = matrix.col_index_;
    matrix_col_val_ = matrix.col_val_;
    matrix_diagonal_val_ = matrix.diagonal_val_;

    matrix_diagonal_val_.back() = - epsilon_ / time;

    long double dummy_initial_nuclide(tot_feeding_rate_ * time / epsilon_);
    initial_n_vector_.back() = dummy_initial_nuclide; // 伪核素的初始值


    // initialize the sort vector
    vector<int> sort_vector;
    sort_vector.resize(initial_n_vector.size());
    for (vector<int>::size_type ix = 0; ix != sort_vector.size(); ++ix)
        sort_vector[ix] = ix;

    // employ the sort algrothim to sort nuclide concentrations and store the sorting ing the sort_vector
    double key;
    int i, key_s;
    for(int j = 1; j < initial_n_vector.size() - 1; ++j) {
        key = initial_n_vector[j];
        key_s = sort_vector[j];
        i = j - 1;
        while ( i > -1 && initial_n_vector[i] < key ) {
            initial_n_vector[i + 1] = initial_n_vector[i];
            sort_vector[i + 1] = sort_vector [i];
            i -= 1;
        }
        initial_n_vector[i + 1] = key;
        sort_vector[i + 1] = key_s;
    }

    for (unsigned int i = 0; i < feed_rate_.size(); ++i) {
        //int nuclide_id = feed_nuclide_id_[i];
        matrix_col_index_.back().push_back(feed_nuclide_id_[i]);
        matrix_col_val_.back().push_back(feed_rate_[i] / dummy_initial_nuclide);
    }

    int index;
    for (int i = 0; i < initial_n_vector_.size() ; ++i) {
        index = sort_vector[i];
        if (initial_n_vector_[index] != 0.0) {
            initial_n_ = initial_n_vector_[index];
            CalOneTree(index, time);
        }
    }
    for (vector<int>::size_type ix = 0; ix != node_visited_list_.size(); ++ix)
        node_visited_list_[ix] = 0;


    for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix) {
        initial_n_vector[ix] = end_n_[ix];
        end_n_[ix] = 0.0;
    }
}
