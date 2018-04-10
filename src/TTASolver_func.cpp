#include "SolveTrans.h"

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
    // 定义一个long double型数组，用来存储母核素相关系数信息
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
