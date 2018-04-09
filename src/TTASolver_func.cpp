#include "SolveTrans.h"

void SolveTrans::SearchOneChainForDepletion(const int & n0_index, const long double & time) {
    // 定义一个元素为vector的栈，表示整个燃耗链
    // 其中，vector第一个元素对应母核素，之后的元素对应各个子核素
    stack < vector<int> > nuclide_chain;
    
    // 定义一个int型数组，用来存储母核素及其子核素id，并且母核素位于数组的第一个位置。
    vector<int> nuclide_group (matrix_col_index_[n0_index]);
    nuclide_group.insert(nuclide_group.begin(), n0_index);

    nuclide_chain.push(nuclide_group);
    while (!nuclide_chain.empty()) {
        int mother_nuclide (nuclide_chain.top()[0]);

        // 将栈顶数组中的最后一个元素作为下一个循环的母核素，同时在栈顶数组中删除该核素
        // 如果数组中只剩下一个元素，意味着母核素的所有子核素都已经找到，
        // 这时候需要从堆栈中弹出，并且直接进入下一个循环
        if (nuclide_chain.top().size() == 1) {
            nuclide_chain.pop();
            continue;
        }
        int next_mother_nuclide (nuclide_chain.top().back());
        nuclide_chain.top().pop_back();
        
        // 将next_mother_nuclide核素对应的nuclide_group压入堆栈
        vector<int> nuclide_group (matrix_col_index_[next_mother_nuclide]);
        nuclide_group.insert(nuclide_group.begin(), next_mother_nuclide);
        nuclide_chain.push(nuclide_group);
    }    
}