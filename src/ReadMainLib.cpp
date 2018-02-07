#include "ModecClass.h"

using namespace std;
void ModecClass::DecayToSpMat() {
    /*fb  10000
    fb1 10000 +1
    fp  -10000
    fp1 -10000 +1
    fa  -20040
    ft  -1
    fsf
    fbn 10000 -10
    fbb 20000
    fn  -10
    fba -10040
    */
    int nuclid, units;
    int _row, _col;
    double halfl, coeff, lamda, q, ampc, wmpc;
    double fb, fb1, fp, fp1, fa, ft, fsf, fbn, fbb, fn, fba;
    int alphaID = ModecNuclideLibrary.GetNuclIndex(20040);

    //int mat_dim = ModecNuclideLibrary.nuclide_library_vector_[0].size();

    ifstream modec_lib;
    char title[100];
    modec_lib.open(decay_library_name_);
    if (!modec_lib) {
        InfoMessage::ErrorMessage("Position: void ModecClass::DecayToSpMat; \n Error: cannot open decay library file of MODEC! go and check the existence of the lib file.", 1);
    }
    modec_lib.getline(title, 100, '\n');
omit:
    while (!modec_lib.eof()) {
        char line[100], line2[100], line3[100];
        char cnucl[8], cunit[6];
        char chalfl[13], cdata[13];
        modec_lib.getline(line, 100, '\n');
        for (int i = 0; i < 8; ++i) {
            cnucl[i] = line[i];
        }
        nuclid = atoi(cnucl);

        for (int i = 0; i < 6; ++i) {
            cunit[i] = line[8 + i];
        }
        units = atoi(cunit);
        if (units == 6) { // 稳定核素，没有衰变反应
            modec_lib.getline(line2, 100, '\n');
            modec_lib.getline(line3, 100, '\n');
            goto omit;
        }

        _row = ModecNuclideLibrary.GetNuclIndex(nuclid);

        if (_row == -1) {
            int i = 1;
        }

        if (units == 1) {
            coeff = log(2);
        } else if (units == 2) {
            coeff = log(2) / 60;
        } else if (units == 3) {
            coeff = log(2) / 3600;
        } else if (units == 4) {
            coeff = log(2) / 86400;
        } else if (units == 5) {
            coeff = log(2) / (365.25 * 86400);
        } else if (units == 7) {
            coeff = log(2) / (365.25 * 86400 * 1.0e3);
        } else if (units == 8) {
            coeff = log(2) / (365.25 * 86400 * 1.0e6);
        } else if (units == 9) {
            coeff = log(2) / (365.25 * 86400 * 1.0e9);
        }


        for (int i = 0; i < 13; ++i) {
            chalfl[i] = line[14 + i];
        }
        halfl = atof(chalfl);
        lamda = coeff / halfl; // 衰变常数lamda
        TransMatrixDecay.AddElement(_row, _row, -lamda);
        ModecNuclideLibrary.nuclide_library_vector_[1][_row] = lamda; // 将核素的衰变常数存储到NuciLib中

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line[27 + i];
        }
        fb1 = atof(cdata);
        if (fb1 != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid + 10000) / 10 * 10 + 1);
            double _val = lamda*fb1;
            TransMatrixDecay.AddElement(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line[40 + i];
        }
        fp = atof(cdata);
        if (fp != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid - 10000) / 10 * 10);
            double _val = lamda*fp;
            TransMatrixDecay.AddElement(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line[53 + i];
        }
        fp1 = atof(cdata);
        if (fp1 != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid - 10000) / 10 * 10 + 1);
            double _val = lamda*fp1;
            TransMatrixDecay.AddElement(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line[66 + i];
        }
        fa = atof(cdata);
        if (fa != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid - 20040) / 10 * 10); // 阿尔法衰变的子核为基态
            double _val = lamda*fa;
            TransMatrixDecay.AddElement(_col, _row, _val);
            TransMatrixDecay.AddElement(alphaID, _row, _val); // 阿尔法衰变产生的He-4元素
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line[79 + i];
        }
        ft = atof(cdata);
        if (ft != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex(nuclid - 1);
            double _val = lamda*ft;
            TransMatrixDecay.AddElement(_col, _row, _val);
        }
        /* 第2排衰变数据 */
        modec_lib.getline(line2, 100, '\n');
        for (int i = 0; i < 13; ++i) {
            cdata[i] = line2[14 + i];
        }
        fsf = atof(cdata);
        ModecNuclideLibrary.nuclide_library_vector_[2][_row] = fsf; // 将核素的自发裂变份额存储到NuciLib中

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line2[27 + i];
        }
        fbn = atof(cdata);
        if (fbn != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid + 10000 - 10) / 10 * 10);
            double _val = lamda*fbn;
            TransMatrixDecay.AddElement(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line2[40 + i];
        }
        q = atof(cdata);
        ModecNuclideLibrary.nuclide_library_vector_[3][_row] = q; // 将核素的衰变热系数存储到NuciLib中

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line2[66 + i];
        }
        ampc = atof(cdata);
        ModecNuclideLibrary.nuclide_library_vector_[4][_row] = ampc; // 将核素的衰变热系数存储到NuciLib中

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line2[79 + i];
        }
        wmpc = atof(cdata);
        ModecNuclideLibrary.nuclide_library_vector_[5][_row] = wmpc; // 将核素的衰变热系数存储到NuciLib中
        /* 第3排衰变数据 */
        modec_lib.getline(line3, 100, '\n');
        for (int i = 0; i < 13; ++i) {
            cdata[i] = line3[27 + i];
        }
        fb = atof(cdata);
        if (fb != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid + 10000) / 10 * 10);
            double _val = lamda*fb;
            TransMatrixDecay.AddElement(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line3[40 + i];
        }
        fbb = atof(cdata);
        if (fbb != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid + 20000) / 10 * 10);
            double _val = lamda*fbb;
            TransMatrixDecay.AddElement(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line3[53 + i];
        }
        fn = atof(cdata);
        if (fn != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid - 10) / 10 * 10);
            double _val = lamda*fn;
            TransMatrixDecay.AddElement(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line3[66 + i];
        }
        fba = atof(cdata);
        if (fba != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid - 10040) / 10 * 10);
            double _val = lamda*fba;
            TransMatrixDecay.AddElement(_col, _row, _val);
            TransMatrixDecay.AddElement(alphaID, _row, _val); // 阿尔法衰变产生的He-4元素
        }

    }

    modec_lib.close();

    if (if_tracking_stockage == true) {
        TransMatrixStockage = TransMatrixDecay;
    }
}

void ModecClass::DecayToSpMatForTta() {
    /*fb  10000
    fb1 10000 +1
    fp  -10000
    fp1 -10000 +1
    fa  -20040
    ft  -1
    fsf
    fbn 10000 -10
    fbb 20000
    fn  -10
    fba -10040
    */
    int nuclid, units;
    int _row, _col;
    double halfl, coeff, lamda, q, ampc, wmpc;
    double fb, fb1, fp, fp1, fa, ft, fsf, fbn, fbb, fn, fba;
    int alphaID = ModecNuclideLibrary.GetNuclIndex(20040);

    //int mat_dim = ModecNuclideLibrary.nuclide_library_vector_[0].size();

    ifstream modec_lib;
    char title[100];
    modec_lib.open(decay_library_name_);
    if (!modec_lib) {
        InfoMessage::ErrorMessage("Position: void ModecClass::DecayToSpMatForTta; \n Error: cannot open decay library file of MODEC! go and check the existence of the lib file.", 1);
    }
    modec_lib.getline(title, 100, '\n');
omit:
    while (!modec_lib.eof()) {
        char line[100], line2[100], line3[100];
        char cnucl[8], cunit[6];
        char chalfl[13], cdata[13];
        modec_lib.getline(line, 100, '\n');
        for (int i = 0; i < 8; ++i) {
            cnucl[i] = line[i];
        }
        nuclid = atoi(cnucl);

        for (int i = 0; i < 6; ++i) {
            cunit[i] = line[8 + i];
        }
        units = atoi(cunit);
        if (units == 6) { // 稳定核素，没有衰变反应
            modec_lib.getline(line2, 100, '\n');
            modec_lib.getline(line3, 100, '\n');
            goto omit;
        }

        _row = ModecNuclideLibrary.GetNuclIndex(nuclid);

        if (_row == -1) {
            int i = 1;
        }

        if (units == 1) {
            coeff = log(2);
        } else if (units == 2) {
            coeff = log(2) / 60;
        } else if (units == 3) {
            coeff = log(2) / 3600;
        } else if (units == 4) {
            coeff = log(2) / 86400;
        } else if (units == 5) {
            coeff = log(2) / (365.25 * 86400);
        } else if (units == 7) {
            coeff = log(2) / (365.25 * 86400 * 1.0e3);
        } else if (units == 8) {
            coeff = log(2) / (365.25 * 86400 * 1.0e6);
        } else if (units == 9) {
            coeff = log(2) / (365.25 * 86400 * 1.0e9);
        }


        for (int i = 0; i < 13; ++i) {
            chalfl[i] = line[14 + i];
        }
        halfl = atof(chalfl);
        lamda = coeff / halfl; // 衰变常数lamda
        TtaMatrixDecay.AddElementCCS(_row, _row, -lamda);
        ModecNuclideLibrary.nuclide_library_vector_[1][_row] = lamda; // 将核素的衰变常数存储到NuciLib中

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line[27 + i];
        }
        fb1 = atof(cdata);
        if (fb1 != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid + 10000) / 10 * 10 + 1);
            double _val = lamda*fb1;
            TtaMatrixDecay.AddElementCCS(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line[40 + i];
        }
        fp = atof(cdata);
        if (fp != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid - 10000) / 10 * 10);
            double _val = lamda*fp;
            TtaMatrixDecay.AddElementCCS(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line[53 + i];
        }
        fp1 = atof(cdata);
        if (fp1 != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid - 10000) / 10 * 10 + 1);
            double _val = lamda*fp1;
            TtaMatrixDecay.AddElementCCS(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line[66 + i];
        }
        fa = atof(cdata);
        if (fa != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid - 20040) / 10 * 10); // 阿尔法衰变的子核为基态
            double _val = lamda*fa;
            TtaMatrixDecay.AddElementCCS(_col, _row, _val);
            TtaMatrixDecay.AddElementCCS(alphaID, _row, _val); // 阿尔法衰变产生的He-4元素
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line[79 + i];
        }
        ft = atof(cdata);
        if (ft != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex(nuclid - 1);
            double _val = lamda*ft;
            TtaMatrixDecay.AddElementCCS(_col, _row, _val);
        }
        /* 第2排衰变数据 */
        modec_lib.getline(line2, 100, '\n');
        for (int i = 0; i < 13; ++i) {
            cdata[i] = line2[14 + i];
        }
        fsf = atof(cdata);
        ModecNuclideLibrary.nuclide_library_vector_[2][_row] = fsf; // 将核素的自发裂变份额存储到NuciLib中

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line2[27 + i];
        }
        fbn = atof(cdata);
        if (fbn != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid + 10000 - 10) / 10 * 10);
            double _val = lamda*fbn;
            TtaMatrixDecay.AddElementCCS(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line2[40 + i];
        }
        q = atof(cdata);
        ModecNuclideLibrary.nuclide_library_vector_[3][_row] = q; // 将核素的衰变热系数存储到NuciLib中

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line2[66 + i];
        }
        ampc = atof(cdata);
        ModecNuclideLibrary.nuclide_library_vector_[4][_row] = ampc; // 将核素的衰变热系数存储到NuciLib中

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line2[79 + i];
        }
        wmpc = atof(cdata);
        ModecNuclideLibrary.nuclide_library_vector_[5][_row] = wmpc; // 将核素的衰变热系数存储到NuciLib中
        /* 第3排衰变数据 */
        modec_lib.getline(line3, 100, '\n');
        for (int i = 0; i < 13; ++i) {
            cdata[i] = line3[27 + i];
        }
        fb = atof(cdata);
        if (fb != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid + 10000) / 10 * 10);
            double _val = lamda*fb;
            TtaMatrixDecay.AddElementCCS(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line3[40 + i];
        }
        fbb = atof(cdata);
        if (fbb != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid + 20000) / 10 * 10);
            double _val = lamda*fbb;
            TtaMatrixDecay.AddElementCCS(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line3[53 + i];
        }
        fn = atof(cdata);
        if (fn != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid - 10) / 10 * 10);
            double _val = lamda*fn;
            TtaMatrixDecay.AddElementCCS(_col, _row, _val);
        }

        for (int i = 0; i < 13; ++i) {
            cdata[i] = line3[66 + i];
        }
        fba = atof(cdata);
        if (fba != 0.0) {
            _col = ModecNuclideLibrary.GetNuclIndex((nuclid - 10040) / 10 * 10);
            double _val = lamda*fba;
            TtaMatrixDecay.AddElementCCS(_col, _row, _val);
            TtaMatrixDecay.AddElementCCS(alphaID, _row, _val); // 阿尔法衰变产生的He-4元素
        }

    }

    modec_lib.close();
}

void ModecClass::XSfromTriton() {
    int _row, _col;
    double _val;
    double ave_E; // 平均裂变能量
    ifstream triton;
    triton.open("TRITON.out");
    if (!triton) {
        InfoMessage::ErrorMessage("Position: void ModecClass::XSfromTriton; \n Error: cannot open TRITON.out! go and check the existence of the file.", 1);
    }
    while ((!triton.fail()) && (!triton.eof())) {
        char buffer[256];
        triton.getline(buffer, 256, '\n');
        /********************* 读入具有裂变份额的核素的平均裂变能量，用于有效裂变份额的计算 *********************/
        if ((buffer[10] == '*') && (buffer[12] == 'A')) {
            for (int i = 0; i < 30; ++i) {
                char data[50];
                triton.getline(data, 256, '\n');
                char c_ave_E[12];
                for (int ii = 0; ii < 12; ++ii) {
                    c_ave_E[ii] = data[26 + ii];
                }
                ave_E = atof(c_ave_E);
                ModecNuclideLibrary.effective_fission_yields_vector_[i].push_back(ave_E);
            }
        }
        /*************************************************************************************************/

        if ((buffer[0] == '1') && (buffer[1] == 'c')) {
            triton.getline(buffer, 256, '\n');
in:
            while (1) {
                char line[256];
                triton.getline(line, 256, '\n');

                if (line[1] == 't')goto out;

                char cnuclid1[7], cnuclid2[7];
                for (int i = 0; i <= 6; ++i) {
                    cnuclid1[i] = line[20 + i];
                }
                int nuclid1 = atoi(cnuclid1);
                if (nuclid1 == 952441) {
                    int i = 1;
                }
                if (nuclid1 > 992550) {
                    goto in;
                }
                _col = ModecNuclideLibrary.GetNuclIndex(nuclid1);
                if (_col == -1) {
                    goto in;
                }
                char ccs[11];
                for (int i = 0; i <= 10; ++i) {
                    ccs[i] = line[39 + i];
                }
                _val = atof(ccs);//*1e-24;
                if (_val == 0) {
                    goto in;
                }
                if (line[30] == 't') {
                    TransMatrixCrossSection.AddElement(_col, _col, -_val);
                    ModecNuclideLibrary.nuclide_library_vector_[8][_col] = _val;
                } else if (line[30] == 's') {
                    ModecNuclideLibrary.nuclide_library_vector_[7][_col] = _val;
                } else if (line[30] == ' ') {
                    for (int i = 0; i <= 6; ++i) {
                        cnuclid2[i] = line[30 + i];
                    }
                    int nuclid2 = atoi(cnuclid2);
                    if (nuclid2 > 992550) {
                        goto in;
                    }
                    _row = ModecNuclideLibrary.GetNuclIndex(nuclid2);
                    if (_row == -1) {
                        goto in;
                    }
                    TransMatrixCrossSection.AddElement(_row, _col, _val);
                }

            }
        }
    }
out:
    triton.close();
    return;
}

void ModecClass::CalculateEffectiveFissionYields() {
    vector<int> nucl_id;
    vector<int> fp_id;
    vector<vector<double> > fissionyield;
    int num_nucl, num_fp;
    char cint[8], cdouble[13];
    string title;
    char line[100];
    ifstream FYLib("BurnFissionYield.lib", ios::in);
    if (!FYLib) {
        InfoMessage::ErrorMessage("Position: void ModecClass::CalculateEffectiveFissionYields; \n Error: cannot open FP-field library file of MODEC! go and check the existence of the lib file.", 1);
    }
    FYLib.getline(line, 100, '\n');

    //      read 'num_nucl' and 'num_fp'     //
    FYLib.getline(line, 100, '\n');
    for (int i = 0; i < 8; ++i) {
        cint[i] = line[i];
    }
    num_nucl = atoi(cint);
    for (int i = 0; i < 8; ++i) {
        cint[i] = line[8 + i];
    }
    num_fp = atoi(cint);
    ////////////////////////////////////////////

    nucl_id.resize(num_nucl);
    fp_id.resize(num_fp);


    ////////// 为有效裂变份额vector定义存储空间 //////////
    //ModecNuclideLibrary.effective_fission_yields_vector_.resize(num_nucl);
    for (unsigned int i = 0; i < ModecNuclideLibrary.effective_fission_yields_vector_.size(); ++i) {
        ModecNuclideLibrary.effective_fission_yields_vector_[i].resize(num_fp + 1); // 留出第一个位置用来存储裂变能量
    }
    ////////////////////////////////////////////////

    int count = 0;
    for (int i = 1; i <= 3; ++i) {
        FYLib.getline(line, 100, '\n');
        for (int j = 0; j < 8; ++j) {
            for (int ii = 0; ii < 8; ++ii) {
                cint[ii] = line[8 * j + ii];
            }
            nucl_id[count] = ModecNuclideLibrary.GetNuclIndex(atoi(cint));
            count++;
        }
    }


    FYLib.getline(line, 100, '\n');
    //count = 0;
    for (int j = 0; j < 6; ++j) {
        for (int ii = 0; ii < 8; ++ii) {
            cint[ii] = line[8 * j + ii];
        }
        nucl_id[count] = ModecNuclideLibrary.GetNuclIndex(atoi(cint));
        count++;
    }
    FYLib.getline(line, 100, '\n');

    count = 0;
    while (1) {
        FYLib.getline(line, 100, '\n');
        if (FYLib.eof()) {
            break;
        }
        char cnum[6];
        for (int ii = 0; ii < 6; ++ii) {
            cnum[ii] = line[8 + ii];
        }
        int num = atoi(cnum);
        for (int i = 0; i < 13; ++i) {
            cdouble[i] = line[13 + i];
        }
        double energy = atof(cdouble);

        fissionyield.resize(num);
        for (int i = 0; i < num; ++i) {
            fissionyield[i].resize(num_fp + 1); // fissionyield[i][0] 用来存储裂变能量
        }

        fissionyield[0][0] = energy;
        int ncount = 1;
        while (1) {
            FYLib.getline(line, 100, '\n');

            if (count == 0) {
                for (int i = 0; i < 8; ++i) {
                    cint[i] = line[i];
                }
                fp_id[ncount - 1] = ModecNuclideLibrary.GetNuclIndex(atoi(cint));
            }

            for (int i = 0; i < 13; ++i) {
                cdouble[i] = line[8 + i];
            }
            fissionyield[0][ncount] = atof(cdouble);
            ncount++;
            if (ncount > num_fp) {
                FYLib.getline(line, 100, '\n');
                break;
            }

            if (count == 0) {
                for (int i = 0; i < 8; ++i) {
                    cint[i] = line[21 + i];
                }
                fp_id[ncount - 1] = ModecNuclideLibrary.GetNuclIndex(atoi(cint));
            }

            for (int i = 0; i < 13; ++i) {
                cdouble[i] = line[29 + i];
            }
            fissionyield[0][ncount] = atof(cdouble);
            ncount++;
            if (ncount > num_fp) {
                FYLib.getline(line, 100, '\n');
                break;
            }

            if (count == 0) {
                for (int i = 0; i < 8; ++i) {
                    cint[i] = line[42 + i];
                }
                fp_id[ncount - 1] = ModecNuclideLibrary.GetNuclIndex(atoi(cint));
            }
            for (int i = 0; i < 13; ++i) {
                cdouble[i] = line[50 + i];
            }
            fissionyield[0][ncount] = atof(cdouble);
            ncount++;
            if (ncount > num_fp) {
                FYLib.getline(line, 100, '\n');
                break;
            }
        }

        //ncount = 1;
        for (int i = 1; i < num; ++i) {
            FYLib.getline(line, 100, '\n');
            for (int ii = 0; ii < 13; ++ii) {
                cdouble[ii] = line[13 + ii];
            }
            double energy = atof(cdouble);
            fissionyield[i][0] = energy;

            ncount = 1;
            while (1) {
                FYLib.getline(line, 100, '\n');
                //if (line[0] == '-') { break; }
                for (int i = 0; i < 13; ++i) {
                    cdouble[i] = line[8 + i];
                }
                fissionyield[i][ncount] = atof(cdouble);
                ncount++;
                if (ncount > num_fp) {
                    FYLib.getline(line, 100, '\n');
                    break;
                }

                for (int i = 0; i < 13; ++i) {
                    cdouble[i] = line[29 + i];
                }
                fissionyield[i][ncount] = atof(cdouble);
                ncount++;
                if (ncount > num_fp) {
                    FYLib.getline(line, 100, '\n');
                    break;
                }

                for (int i = 0; i < 13; ++i) {
                    cdouble[i] = line[50 + i];
                }
                fissionyield[i][ncount] = atof(cdouble);
                ncount++;
                if (ncount > num_fp) {
                    FYLib.getline(line, 100, '\n');
                    break;
                }
            }

        }

        ///////////////////////////////////////////////////////////////////
        //    calculate effective fission yields  //
        double ave_FE = ModecNuclideLibrary.effective_fission_yields_vector_[count][0];
        if (ave_FE <= fissionyield[0][0]) {
            for (int i = 1; i <= num_fp; ++i) {
                ModecNuclideLibrary.effective_fission_yields_vector_[count][i] = fissionyield[0][i];
            }
        } else if (ave_FE >= fissionyield[num - 1][0]) {
            for (int i = 1; i <= num_fp; ++i) {
                ModecNuclideLibrary.effective_fission_yields_vector_[count][i] = fissionyield[num - 1][i];
            }
        } else {
            if (num == 2) {
                double Eb = fissionyield[0][0];
                double Ee = fissionyield[num - 1][0];
                double E = ave_FE;
                for (int i = 1; i <= num_fp; ++i) {
                    double coeff = (fissionyield[num - 1][i] - fissionyield[0][i])
                                   / (fissionyield[num - 1][0] - fissionyield[0][0]);
                    double b = fissionyield[0][i] - coeff * fissionyield[0][0];

                    ModecNuclideLibrary.effective_fission_yields_vector_[count][i] =
                        coeff*ModecNuclideLibrary.effective_fission_yields_vector_[count][0] + b;
                }
            } else if (num == 3) {
                if (ave_FE > fissionyield[0][0] && ave_FE <= fissionyield[1][0]) {
                    double Eb = fissionyield[0][0];
                    double Ee = fissionyield[1][0];
                    double E = ave_FE;
                    for (int i = 1; i <= num_fp; ++i) {
                        double coeff = (fissionyield[1][i] - fissionyield[0][i])
                                       / (fissionyield[1][0] - fissionyield[0][0]);
                        double b = fissionyield[0][i] - coeff * fissionyield[0][0];

                        ModecNuclideLibrary.effective_fission_yields_vector_[count][i] =
                            coeff*ModecNuclideLibrary.effective_fission_yields_vector_[count][0] - b;
                    }
                } else {
                    double Eb = fissionyield[1][0];
                    double Ee = fissionyield[num - 1][0];
                    double E = ave_FE;
                    for (int i = 1; i <= num_fp; ++i) {
                        double coeff = (fissionyield[num - 1][i] - fissionyield[1][i])
                                       / (fissionyield[num - 1][0] - fissionyield[1][0]);
                        double b = fissionyield[1][i] - coeff * fissionyield[1][0];

                        ModecNuclideLibrary.effective_fission_yields_vector_[count][i] =
                            coeff*ModecNuclideLibrary.effective_fission_yields_vector_[count][0] - b;
                    }
                }
            }
        }
        count++;
        if (count == 30) {
            int i = 1;
        }
        //if (count >= 30) {
        //	break; }
    }
    FYLib.close();

    /////////////////////////////////////////////////////////////////////////////////
    ////////////  将裂变产物份额对应的裂变产物产生项填入裂变矩阵中		/////////////////
    for (unsigned int i = 0; i < ModecNuclideLibrary.effective_fission_yields_vector_.size(); ++i) {
        int _col = nucl_id[i];
        double fissionXs = ModecNuclideLibrary.nuclide_library_vector_[7][i]; // 母核裂变截面
        for (unsigned int j = 0; j < fp_id.size(); ++j) {
            int _row = fp_id[j];
            double _val = fissionXs*ModecNuclideLibrary.effective_fission_yields_vector_[i][j + 1];
            TransMatrixCrossSection.AddElement(_row, _col, _val);
        }
    }
    return;
}
