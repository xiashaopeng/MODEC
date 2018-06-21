#include "ModecClass.h"

using namespace std;

//////////////////////////////////   在线后处理的系数加入   /////////////////////////////////////////////////
void ModecClass::AddOnlineReprocessingCoeffi() {
    if (solver_selection_ != 0) {
        int group_num = remove_rate_vector_.size();
        vector<int> eleList;
        for (int i = 0; i < group_num; ++i) {
            int ele_num = remove_element_vector_[i].size();
            for (int j = 0; j < ele_num; ++j) {
                int ID = remove_element_vector_[i][j];
                eleList = ModecNuclideLibrary.GetEleIndex(ID);

                for (unsigned int ii = 0; ii < eleList.size(); ++ii) {
                    TransMatrixDecay.AddElement(eleList[ii], eleList[ii], -remove_rate_vector_[i]);
                    if (if_keeping_eutectic_stable_ == true && ID < 90) {
                        int Li7Index = ModecNuclideLibrary.GetNuclIndex(30070);
                        int Li6Index = ModecNuclideLibrary.GetNuclIndex(30060);
                        TransMatrixDecay.AddElement(Li7Index, eleList[ii], remove_rate_vector_[i] * 0.99995);
                        TransMatrixDecay.AddElement(Li6Index, eleList[ii], remove_rate_vector_[i] * 0.00005);
                    }
                    if (if_variable_feeding_ == true && ID >= 90) {
                        int feed_size = variable_feeding_nuclide_id_vector_.size(); // 应该等于2
                        for (int a = 0; a < feed_size; ++a) {
                            double feedratio;
                            if (a == 0) {
                                feedratio = variable_feeding_ratio_;
                            }
                            if (a == 1) {
                                feedratio = 1 - variable_feeding_ratio_;
                            }
                            for (unsigned int aa = 0; aa < variable_feeding_nuclide_id_vector_[a].size(); ++aa) {
                                int Feed_ID = variable_feeding_nuclide_id_vector_[a][aa];
                                double _val = variable_feeding_nuclide_ratio_vector_[a][aa] * feedratio * remove_rate_vector_[i];
                                TransMatrixDecay.AddElement(ModecNuclideLibrary.GetNuclIndex(Feed_ID), eleList[ii], _val);
                            }
                        }
                    }

                    if (if_tracking_stockage == true) {
                        TransMatrixReprocess.AddElement(eleList[ii], eleList[ii], remove_rate_vector_[i]);
                    }
                }
            }
        }

    } else {
        int group_num = remove_rate_vector_.size();
        vector<int> eleList;
        for (int i = 0; i < group_num; ++i) {
            int ele_num = remove_element_vector_[i].size();
            for (int j = 0; j < ele_num; ++j) {
                int ID = remove_element_vector_[i][j];
                eleList = ModecNuclideLibrary.GetEleIndex(ID);

                for (unsigned int ii = 0; ii < eleList.size(); ++ii) {
                    TtaMatrixDecay.AddElementCCS(eleList[ii], eleList[ii], -remove_rate_vector_[i]);
                    if (if_keeping_eutectic_stable_ == true && ID < 90) {
                        int Li7Index = ModecNuclideLibrary.GetNuclIndex(30070);
                        int Li6Index = ModecNuclideLibrary.GetNuclIndex(30060);
                        TtaMatrixDecay.AddElementCCS(Li7Index, eleList[ii], remove_rate_vector_[i] * 0.99995);
                        TtaMatrixDecay.AddElementCCS(Li6Index, eleList[ii], remove_rate_vector_[i] * 0.00005);
                    }
                    if (if_variable_feeding_ == true && ID >= 90) {
                        int feed_size = variable_feeding_nuclide_id_vector_.size(); // 应该等于2
                        for (int a = 0; a < feed_size; ++a) {
                            double feedratio;
                            if (a == 0) {
                                feedratio = variable_feeding_ratio_;
                            }
                            if (a == 1) {
                                feedratio = 1 - variable_feeding_ratio_;
                            }
                            for (unsigned int aa = 0; aa < variable_feeding_nuclide_id_vector_[a].size(); ++aa) {
                                int Feed_ID = variable_feeding_nuclide_id_vector_[a][aa];
                                double _val = variable_feeding_nuclide_ratio_vector_[a][aa] * feedratio * remove_rate_vector_[i];
                                TtaMatrixDecay.AddElementCCS(ModecNuclideLibrary.GetNuclIndex(Feed_ID), eleList[ii], _val);
                            }
                        }
                    }

                    if (if_tracking_stockage == true) {
                        TransMatrixReprocess.AddElement(eleList[ii], eleList[ii], remove_rate_vector_[i]);
                    }
                }
            }
        }

    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////   在线添料系数的加入考虑   /////////////////////////////////////////////////

// 1. 将添料率常数转换成其次项（结合液态熔盐堆重金属浓度守恒以及临界运行的特点）
//    添料核素的总添料率等于系统的总的裂变率

void ModecClass::ContinuouslyFeeding() {

    int Li7Index, Li6Index;
    double Li7con, Li6con, Li7ratio;
    if (if_keeping_eutectic_stable_ == true) {
        Li7Index = ModecNuclideLibrary.GetNuclIndex(30070);
        Li6Index = ModecNuclideLibrary.GetNuclIndex(30060);
        Li7con = ModecNuclideLibrary.nuclide_library_vector_[0][Li7Index];
        Li6con = ModecNuclideLibrary.nuclide_library_vector_[0][Li6Index];

        Li7ratio = Li7con / (Li7con + Li6con);
    }

    int size = ModecNuclideLibrary.nuclide_library_vector_[7].size(); // 裂变截面向量
    for (int i = 0; i < size; ++i) {
        if (ModecNuclideLibrary.nuclide_library_vector_[7][i] > 0.0) {
            int feed_size = variable_feeding_nuclide_id_vector_.size(); // 应该等于2
            for (int a = 0; a < feed_size; ++a) {
                double feedratio;
                if (a == 0) {
                    feedratio = variable_feeding_ratio_;
                }
                if (a == 1) {
                    feedratio = 1 - variable_feeding_ratio_;
                }
                for (unsigned int aa = 0; aa < variable_feeding_nuclide_id_vector_[a].size(); ++aa) {
                    int Feed_ID = variable_feeding_nuclide_id_vector_[a][aa];
                    int FeedIndex = ModecNuclideLibrary.GetNuclIndex(Feed_ID);
                    double _val = variable_feeding_nuclide_ratio_vector_[a][aa] * feedratio * ModecNuclideLibrary.nuclide_library_vector_[7][i];
                    TransMatrixCrossSection.AddElement(FeedIndex, i, _val);
                }
            }

            if (if_keeping_eutectic_stable_ == true) {
                TransMatrixCrossSection.AddElement(Li7Index, i, -ModecNuclideLibrary.nuclide_library_vector_[7][i] * 2 * Li7ratio); // 认为裂变产额为2
                TransMatrixCrossSection.AddElement(Li6Index, i, -ModecNuclideLibrary.nuclide_library_vector_[7][i] * 2 * (1 - Li7ratio));
            }
        }
    }
}
