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

void SolveTrans::SearchOneChainForDepletion(const int & n0_index, const long double & time)
{
	vector<int> next_index_list(matrix_col_index_[n0_index]);
	long double lamda_ii(-matrix_diagonal_val_[n0_index]);
	if (next_index_list.size() == 0)
	{
		return;
	}
	if (next_index_list.size() > 1)
	{
		node_a_.push(chain_a_);
		node_b_.push(chain_b_);
		node_gamma_.push(chain_gamma_);
	}

	for (int i = 0; i < next_index_list.size(); ++i)
	{
		if (i > 0)
		{
			chain_a_ = node_a_.top();
			chain_b_ = node_b_.top();
			chain_gamma_ = node_gamma_.top();

			int size = chain_lamda_list_.size();
			for (int j = chain_a_; j < size; ++j)
			{
				chain_lamda_list_.pop_back();
				chain_nuclide_id.pop_back();
			}
		}

		int next_daughter_index(next_index_list[i]);
		
		node_visited_list_[next_daughter_index] += 1;
		long double n0_i1(0.0);
		if (node_visited_list_[next_daughter_index] == 1)
		{
			n0_i1 = initial_n_vector_[next_daughter_index];
		}
		else
		{
			n0_i1 = 0;
		}

		long double lamda_jj(-matrix_diagonal_val_[next_daughter_index]);
		long double beff_ji(matrix_col_val_[n0_index][i]);

		int visited_number(0);
		for (int ii = 0; ii < chain_nuclide_id.size(); ++ii)
		{
			if (next_daughter_index == chain_nuclide_id[ii]) // 这里的lamda系数的处理可能存在问题
			{
				visited_number++;
			}
		}

		if (visited_number > 0 || chain_nuclide_id.size()>=5)
		{
			long double judge_cutoff(CalCutoffFlag(next_daughter_index, lamda_ii, time));
			n0_i1 = 0;
			if (abs(judge_cutoff) <= cutoff_ && next_index_list.size() > 1)
			{
				node_a_.pop();
				node_b_.pop();
				node_gamma_.pop();
				
				return;
			}
			if (abs(judge_cutoff) <= cutoff_ && next_index_list.size() == 1)
			{
				return;
			}
		}

		long double temp_n(0.0); // 子核的核素浓度
		long double temp_gamma_j0(0.0);
		int equal_j(-1);

		

		for (int j = 0; j < chain_a_; ++j)
		{
			if (lamda_jj == chain_lamda_list_[j])
			{
				equal_j = j;
			}
			else
			{
				for (int k = chain_b_[j]; k >= 0; --k)
				{
					if (k == chain_b_[j])
					{
						chain_gamma_[j][k] = chain_gamma_[j][k] * beff_ji / (lamda_jj - chain_lamda_list_[j]);
					}
					else
					{
						chain_gamma_[j][k] = (chain_gamma_[j][k] * beff_ji - (k + 1)*chain_gamma_[j][k + 1])
							/ (lamda_jj - chain_lamda_list_[j]);
					}
					if (k == 0)
					{
						temp_gamma_j0 += chain_gamma_[j][k];
					}
					temp_n += chain_gamma_[j][k] * pow(time, k)*exp(-chain_lamda_list_[j] * time);
				}
			}
		}

		if (equal_j == -1)
		{
			chain_a_ += 1;
			chain_b_.push_back(0);
			chain_gamma_.push_back({  n0_i1 - temp_gamma_j0 });
			chain_lamda_list_.push_back(lamda_jj);

			chain_nuclide_id.push_back(next_daughter_index);
			temp_n += chain_gamma_[chain_a_ - 1][0] * exp(-lamda_jj * time);
		}
		else
		{
			vector<long double> chain_gamma(chain_gamma_[equal_j]);

			chain_b_[equal_j] += 1;
			for (int k = 0; k <= chain_b_[equal_j]; ++k)
			{
				if (k == chain_b_[equal_j])
				{
					chain_gamma_[equal_j].push_back(beff_ji*chain_gamma[k - 1] / k);
					temp_n += chain_gamma_[equal_j].back() * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
				}
				else if (k > 0 && k < chain_b_[equal_j])
				{
					chain_gamma_[equal_j][k] = beff_ji*chain_gamma[k - 1] / k;
					temp_n += chain_gamma_[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
				}
				else
				{
					chain_gamma_[equal_j][k] = n0_i1 - temp_gamma_j0;
					temp_n += chain_gamma_[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
				}
			}
		}
		if (temp_n >= 0)
		{
			end_n_[next_daughter_index] += temp_n;
		}
		SearchOneChainForDepletion(next_daughter_index, time);
	}

	if (next_index_list.size() > 1)
	{
		node_a_.pop();
		node_b_.pop();
		node_gamma_.pop();
	}
		
	return;

}

long double SolveTrans::CalCutoffFlag(const int & next_daughter_index, const long double & beff_ji, const long double & time)
{
	long double lamda_jj = 0.0;// -matrix_diagonal_val_[next_daughter_index];
	//long double alpha_new = 1;
	//long double alpha_sum = 0;
	long double n0_i1;
	if (node_visited_list_[next_daughter_index] == 1)
	{
		n0_i1 = initial_n_vector_[next_daughter_index];
	}
	else
	{
		n0_i1 = 0;
	}
	int chain_a = chain_a_;
	vector<int> chain_b = chain_b_;
	vector<vector<long double> > chain_gamma = chain_gamma_;

	long double temp_n = 0; // 子核的核素浓度
	long double temp_gamma_j0 = 0;
	int equal_j = -1;
	for (int j = 0; j < chain_a; ++j)
	{
		if (lamda_jj == chain_lamda_list_[j])
		{
			equal_j = j;
		}
		else
		{
			for (int k = chain_b[j]; k >= 0; --k)
			{
				if (k == chain_b[j])
				{
					chain_gamma[j][k] = chain_gamma[j][k] * beff_ji / (lamda_jj - chain_lamda_list_[j]);
				}
				else
				{
					chain_gamma[j][k] = (chain_gamma[j][k] * beff_ji - (k + 1)*chain_gamma[j][k + 1])
						/ (lamda_jj - chain_lamda_list_[j]);
				}
				if (k == 0)
				{
					temp_gamma_j0 += chain_gamma[j][k];
				}
				temp_n += chain_gamma[j][k] * pow(time, k)*exp(-chain_lamda_list_[j] * time);
			}
		}
	}

	if (equal_j == -1)
	{
		temp_n += (n0_i1 - temp_gamma_j0) * exp(-lamda_jj * time);
	}
	else
	{
		chain_b[equal_j] += 1;
		for (int k = 0; k <= chain_b[equal_j]; ++k)
		{
			if (k == chain_b[equal_j])
			{
				chain_gamma[equal_j].push_back(beff_ji*chain_gamma_[equal_j][k - 1] / k);
				temp_n += chain_gamma[equal_j].back() * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
			}
			else if (k > 0 && k < chain_b[equal_j])
			{
				chain_gamma[equal_j][k] = beff_ji*chain_gamma_[equal_j][k - 1] / k;
				temp_n += chain_gamma[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
			}
			else
			{
				chain_gamma[equal_j][k] = n0_i1 - temp_gamma_j0;
				temp_n += chain_gamma[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
			}
		}
	}
	if(temp_n > 0.0)
	{
		return temp_n;
	}
	else
	{
		return 0.0;
	}
}

void SolveTrans::CalOneTree(const int &n0_index, const long double &time)
{
	//--------------------- 计算树根核素的alpha系数和beff系数 ---------------------//
	long double lamda = -matrix_diagonal_val_[n0_index];
	//-------------------------------------------------------------------------------//
	chain_a_ = 1;
	chain_b_.push_back(0);
	chain_gamma_.push_back({ initial_n_ });
	//-------------------------------------------------------------------------------//


	int chain_point_count = 0; // 计数，如果找到一个子核，加1

	node_visited_list_[n0_index] += 1;
	if (node_visited_list_[n0_index] > 1)
	{
		chain_b_.clear();
		chain_gamma_.clear();
		chain_lamda_list_.clear();
		chain_nuclide_id.clear();

		node_a_.empty();
		node_b_.empty();
		node_gamma_.empty();

		return;
	}

	end_n_[n0_index] += initial_n_ * exp(-lamda*time); // 得到树根核素的浓度

	chain_lamda_list_.push_back(lamda);
	
	chain_nuclide_id.push_back(n0_index);
	//----------------------------------------------------------------------------------//

	//--------------------循环方法实现树的搜索 ----------------------- //

	//int stack_depth(1); // 栈的深度初始化为1
	//int stack_upper;
	//vector<int> stack_next_to_upper;
	//vector<vector<int> > stack_next_to_upper_list;
	//vector<vector<long double> > stack_chain_gamma;
	//vector<vector<vector<long double> > > chain_gamma_list;
	//chain_gamma_list.push_back(chain_gamma_);
	//vector<int> stack_chain_b;
	//vector< vector<int> > chain_b_list;
	//chain_b_list.push_back(chain_b_);
	//int stack_chain_a;
	//vector< int > chain_a_list;
	//chain_a_list.push_back(chain_a_);
	//do 
	//{
	//	stack_upper = chain_nuclide_id.back();
	//	stack_next_to_upper = matrix_col_index_[stack_upper]; // 找到与栈顶元素相邻的元素集
	//	stack_chain_gamma = chain_gamma_list.back();
	//	stack_chain_b = chain_b_list.back();
	//	stack_chain_a = chain_a_list.back();
	//	if (stack_next_to_upper[0] == -1) // 不存在相邻元素
	//	{
	//		// stack_upper出栈
	//		chain_a_ = node_a_.top();
	//		chain_b_ = node_b_.top();
	//		chain_gamma_ = node_gamma_.top();
	//		int size = chain_lamda_list_.size();
	//		for (int j = chain_a_; j < size; ++j)
	//		{
	//			chain_lamda_list_.pop_back();
	//			chain_nuclide_id.pop_back();
	//		}
	//	}
	//	else
	//	{
	//		//node_visited_list_[stack_upper] += 1;
	//		for (int i = 0; i < stack_next_to_upper.size(); ++i)
	//		{
	//			int next_daughter_index(stack_next_to_upper[i]);
	//			if (node_visited_list_[next_daughter_index] == 1) //已经访问过该子节点
	//			{
	//				continue;
	//			}
	//			else // 第一次访问该子节点
	//			{
	//				node_visited_list_[next_daughter_index] += 1;
	//				long double lamda_ii(-matrix_diagonal_val_[stack_upper]); // 母核的lamda
	//				long double lamda_jj(-matrix_diagonal_val_[next_daughter_index]); // 待求子核素的lamda
	//				long double lamda_ij(matrix_col_val_[stack_upper][i]); // 母核到子核的产生率系数
	//				long double n0_i1(initial_n_vector_[next_daughter_index]); // n0_i1为待求子核素的初始浓度，如果访问过则为0
	//				if (node_visited_list_[next_daughter_index] == 1)
	//				{
	//					n0_i1 = initial_n_vector_[next_daughter_index];
	//				}
	//				else
	//				{
	//					n0_i1 = 0;
	//				}
	//				int visited_number(0);
	//				for (int ii = 0; ii < chain_nuclide_id.size(); ++ii)
	//				{
	//					if (next_daughter_index == chain_nuclide_id[ii]) // 这里的lamda系数的处理可能存在问题
	//					{
	//						visited_number++;
	//					}
	//				}
	//				if (visited_number > 0 || chain_nuclide_id.size() >= 5)
	//				{
	//					long double judge_cutoff(CalCutoffFlag(next_daughter_index, lamda_ii, time));
	//					//n0_i1 = 0;
	//					if (abs(judge_cutoff) <= cutoff_ )
	//					{
	//						continue;
	//					}
	//				}
	//				long double temp_n(0.0); // 待求子核的核素浓度
	//				long double temp_gamma_j0(0.0);
	//				int equal_j(-1);
	//				for (int j = 0; j < stack_chain_a; ++j)
	//				{
	//					if (lamda_jj == chain_lamda_list_[j])
	//					{
	//						equal_j = j;
	//					}
	//					else
	//					{
	//						for (int k = stack_chain_b[j]; k >= 0; --k)
	//						{
	//							if (k == stack_chain_b[j])
	//							{
	//								stack_chain_gamma[j][k] = stack_chain_gamma[j][k] * lamda_ij / (lamda_jj - chain_lamda_list_[j]);
	//							}
	//							else
	//							{
	//								stack_chain_gamma[j][k] = (stack_chain_gamma[j][k] * lamda_ij - (k + 1)*stack_chain_gamma[j][k + 1])
	//									/ (lamda_jj - chain_lamda_list_[j]);
	//							}
	//							if (k == 0)
	//							{
	//								temp_gamma_j0 += stack_chain_gamma[j][k];
	//							}
	//							temp_n += stack_chain_gamma[j][k] * pow(time, k)*exp(-chain_lamda_list_[j] * time);
	//						}
	//					}
	//				}
	//				if (equal_j == -1)
	//				{
	//					stack_chain_a += 1;
	//					stack_chain_b.push_back(0);
	//					stack_chain_gamma.push_back({ n0_i1 - temp_gamma_j0 });
	//					chain_lamda_list_.push_back(lamda_jj);
	//					chain_nuclide_id.push_back(next_daughter_index);
	//					temp_n += stack_chain_gamma[stack_chain_a - 1][0] * exp(-lamda_jj * time);
	//					//temp_nn_vec.push_back((n0_i1 - temp_test2) * exp(-lamda_jj * time));
	//				}
	//				else
	//				{
	//					vector<long double> chain_gamma(stack_chain_gamma[equal_j]);
	//					stack_chain_b[equal_j] += 1;
	//					for (int k = 0; k <= stack_chain_b[equal_j]; ++k)
	//					{
	//						if (k == stack_chain_b[equal_j])
	//						{
	//							stack_chain_gamma[equal_j].push_back(lamda_ij*chain_gamma[k - 1] / k);
	//							temp_n += stack_chain_gamma[equal_j].back() * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
	//							//temp_nn_vec.push_back(chain_gamma_[equal_j][k] * pow(time, long double(k))*exp(-chain_lamda_list_[equal_j] * time));
	//						}
	//						else if (k > 0 && k < stack_chain_b[equal_j])
	//						{
	//							stack_chain_gamma[equal_j][k] = lamda_ij*chain_gamma[k - 1] / k;
	//							temp_n += stack_chain_gamma[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
	//							//temp_nn_vec.push_back(chain_gamma_[equal_j][k] * pow(time, long double(k))*exp(-chain_lamda_list_[equal_j] * time));
	//						}
	//						else
	//						{
	//							stack_chain_gamma[equal_j][k] = n0_i1 - temp_gamma_j0;
	//							temp_n += stack_chain_gamma[equal_j][k] * pow(time, k)*exp(-chain_lamda_list_[equal_j] * time);
	//							//temp_nn_vec.push_back(chain_gamma_[equal_j][k] * pow(time, long double(k))*exp(-chain_lamda_list_[equal_j] * time));
	//						}
	//					}
	//				}
	//				if (temp_n > 0)
	//				{
	//					end_n_[next_daughter_index] += temp_n;
	//				}
	//				chain_a_list.push_back(stack_chain_a);
	//				chain_b_list.push_back(stack_chain_b);
	//				chain_gamma_list.push_back(stack_chain_gamma);
	//				chain_nuclide_id.push_back(next_daughter_index); // 将该子节点的核素ID压入堆栈中
	//				chain_lamda_list_.push_back(lamda_ii);
	//				break; // 跳出循环
	//			}
	//		}
	//	}
	//} while (stack_depth == 0); // 栈的长度小于0时，意味着完成计算
	//stack<int> stack_node_nuclide; // 分叉点的核素存储栈
	//stack<int> stack_point_nuclide; // 燃耗子链上的核素存储栈
	//--------------------------------------------------------------//

	// ----------------- 递归方法实现树的搜索 ------------------ //
	SearchOneChainForDepletion(n0_index, time); //

	//------------------------------- 完成一条燃耗树的计算，清理相关向量 ------------------------------//

	//node_matirx_index_list_.clear(); // 存储节点核素在矩阵中的列号
	//node_chain_index_list_.clear();  // 存储节点在链中的序号
	//node_daughter_number_list_.clear(); // 存储节点的子核个数
	//node_beff_list_.clear(); // beff(0) = 1 ; beff(k+1) = beff(k)*lamdaji
	//node_alpha_list_.clear(); // alpha(0) = exp(-lamda_0*t) ; alpha(k+1,i)=alpha(k,i)*1/(lamda_k+1 - lamda_i), alpha(k+1,k+1)
	chain_lamda_list_.clear();
	//chain_beff_list_.clear();
	//chain_alpha_list_.clear();
	chain_nuclide_id.clear();
	//--------------------------------------------------------------------------//
	chain_b_.clear();
	chain_gamma_.clear();
	node_a_.empty();
	node_b_.empty();
	node_gamma_.empty();
}

vector<double> SolveTrans::TtaSolver(const SparseMatrixMCS &matrix,const double &time )
{
	// for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix)
		// end_n_[ix] = 0.0;

	matrix_col_index_ = matrix.col_index_;
	matrix_col_val_ = matrix.col_val_;
	matrix_diagonal_val_ = matrix.diagonal_val_;

	for (int i = 0; i < initial_n_vector_.size(); ++i)
	{
		if (initial_n_vector_[i] != 0.0)
		{
			initial_n_ = initial_n_vector_[i];
			CalOneTree(i, time);
		}
	}
	for (vector<int>::size_type ix = 0; ix != node_visited_list_.size(); ++ix)
		node_visited_list_[ix] = 0;

	vector<double> end_out;
	end_out.resize(end_n_.size());
	for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix)
	{
		end_out[ix] = end_n_[ix];
		end_n_[ix] = 0.0;
	}
		
	return end_out;
}

vector<double> SolveTrans::TtaSolverForFeeding(const SparseMatrixMCS &matrix, const double &time)
{
	// for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix)
	// end_n_[ix] = 0.0;

	matrix_col_index_ = matrix.col_index_;
	matrix_col_val_ = matrix.col_val_;
	matrix_diagonal_val_ = matrix.diagonal_val_;

	matrix_diagonal_val_.back() = - epsilon_ / time;

	long double dummy_initial_nuclide(tot_feeding_rate_ * time / epsilon_);
	initial_n_vector_[initial_n_vector_.size() - 1] = dummy_initial_nuclide; // 伪核素的初始值


	for (unsigned int i = 0; i < feed_rate_.size(); ++i)
	{
		//int nuclide_id = feed_nuclide_id_[i];
		matrix_col_index_.back().push_back(feed_nuclide_id_[i]);
		matrix_col_val_.back().push_back(feed_rate_[i] / dummy_initial_nuclide);
	}

	int size = initial_n_vector_.size();
	for (int i = 0; i < size ; ++i)
	{
		if (initial_n_vector_[i] != 0.0)
		{
			initial_n_ = initial_n_vector_[i];
			CalOneTree(i, time);
		}
	}
	for (vector<int>::size_type ix = 0; ix != node_visited_list_.size(); ++ix)
		node_visited_list_[ix] = 0;

	vector<double> end_out;
	end_out.resize(end_n_.size());
	for (vector<int>::size_type ix = 0; ix != end_n_.size(); ++ix)
	{
		end_out[ix] = end_n_[ix];
		end_n_[ix] = 0.0;
	}

	return end_out;
}