#include "ModecClass.h"

using namespace std;

////////////////////////////////////////////读取Depth库///////////////////////////////////////////////

/// 既读取衰变矩阵，也读取裂变矩阵
void ModecClass::ReadFromDepthLib() {
    char line[200];
    int NuclID, ProdID, itemp;
    double dtemp;
    int _col, _row;
    double lamda;
    int alphaID = ModecNuclideLibrary.GetNuclIndex(20040);
    int protonID = ModecNuclideLibrary.GetNuclIndex(10010);
    int FisNucNum, PrdctNum; // 有裂变场的裂变核素个数及裂变产物个数
    int FPid; // 裂变产物ID
    double FPyield; // 裂变产物份额
    //vector<int> FisNucIndex;

    double halfl, fb1, fp, fp1, fa, ft;       //// first line
    double fsf, fbn, fb, fbb, fneut, fba;    //// second line
    double q, abund, ampc, wmpc;			//// third line (only for ORIGEN-s decay library)
    double SNG, SN2N, SN3N_SNA, SNF_SNP, SNGX, SN2NX;       //// fourth line cross-section
    double fisrate;

    //int mat_dim = ModecNuclideLibrary.nuclide_library_vector_[0].size(); // 矩阵维度

    ifstream ReadDepth;
    ReadDepth.open(depth_library_name_, ios::in);
    if (!ReadDepth) {
        InfoMessage::ErrorMessage( "Position: void ModecClass::ReadFromDepthLib; \n Error: cannot open DEPTH library file! go and check the existence of the lib file.",1);
    }
    while (!ReadDepth.eof()) {
        if (ReadDepth.peek() == '*' || ReadDepth.peek() == '!' || ReadDepth.peek() == '/') {
            ReadDepth.getline(line, 200);
        } else {
            break;
        }
    }

    // read decay and cross-section data from library
    while (!ReadDepth.eof()) {
        ReadDepth >> NuclID;
        if (NuclID == -1) {
            ReadDepth.getline(line, 200);
            break;
        }
        _col = ModecNuclideLibrary.GetNuclIndex(NuclID);
        if (_col == -1) {
            stringstream ss;
            string id;
            ss << NuclID;
            ss >> id;
            info_message_ = "Position: void ModecClass::ReadFromDepthLib; \n Warning: unknown nuclide id " + id + ".";
            InfoMessage::ErrorMessage(info_message_,0);
        }
        if (!(ReadDepth >> line)) {
            info_message_ = "Position: void ModecClass::ReadFromDepthLib; \n Warning: unknown flag " + string(line)+ ".";
            InfoMessage::ErrorMessage(info_message_, 0);
        }
        /// first line ///
        ReadDepth >> lamda >> fb1 >> fp >> fp1 >> fa >> ft;

        /// second line ///
        ReadDepth >> fsf >> fbn >> fb >> fbb >> fneut >> fba;

        /// third line ///
        ReadDepth >> q >> abund >> ampc >> wmpc >> dtemp >> dtemp;

        /// fourth line ///
        ReadDepth >> SNG >> SN2N >> SN3N_SNA >> SNF_SNP >> SNGX >> SN2NX;
        double totXS = SNG + SN2N + SN3N_SNA + SNF_SNP + SNGX + SN2NX;

        /// tag line ///
        ReadDepth >> itemp;

        if (itemp != -1) {
            InfoMessage::ErrorMessage("Position: void ModecClass::ReadFromDepthLib; \n Error: error flag in DEPTH lib!",1);
        }

        // 衰变矩阵TransMatrixDecay建立
        if (lamda != 0.0 && _col != -1) {
            //lamda = log(2) / halfl;

            if (fb1 != 0.0) {
                ProdID = (NuclID + 10000) / 10 * 10 + 1;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fb1;
                TransMatrixDecay.AddElement(_row, _col, _val);
            }

            if (fp != 0.0) {
                ProdID = (NuclID - 10000) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fp;
                TransMatrixDecay.AddElement(_row, _col, _val);
            }

            if (fp1 != 0.0) {
                ProdID = (NuclID - 10000) / 10 * 10 + 1;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fp1;
                TransMatrixDecay.AddElement(_row, _col, _val);
            }

            if (fa != 0.0) {
                ProdID = (NuclID - 20040) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fa;
                TransMatrixDecay.AddElement(_row, _col, _val);

                TransMatrixDecay.AddElement(alphaID, _col, _val); // 阿尔法衰变产生的He-4元素
            }

            if (ft != 0.0) {
                ProdID = NuclID - 1;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*ft;
                TransMatrixDecay.AddElement(_row, _col, _val);
            }

            if (fbn != 0.0) {
                ProdID = (NuclID + 10000 - 10) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fbn;
                TransMatrixDecay.AddElement(_row, _col, _val);
            }

            if (fb != 0.0) {
                ProdID = (NuclID + 10000) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fb;
                TransMatrixDecay.AddElement(_row, _col, _val);
            }

            if (fbb != 0.0) {
                ProdID = (NuclID + 20000) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fbb;
                TransMatrixDecay.AddElement(_row, _col, _val);
            }

            if (fneut != 0.0) {
                ProdID = (NuclID - 10) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fneut;
                TransMatrixDecay.AddElement(_row, _col, _val);
            }

            if (fba != 0.0) {
                ProdID = (NuclID - 10040) / 10 * 10;
                if (ProdID > 992550)
                    _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fba;
                TransMatrixDecay.AddElement(_row, _col, _val);
                TransMatrixDecay.AddElement(alphaID, _col, _val); // 阿尔法衰变产生的He-4元素
            }

            if (_col != -1) {
                TransMatrixDecay.AddElement(_col, _col, -lamda);
                ModecNuclideLibrary.nuclide_library_vector_[1][_col] = lamda;
                ModecNuclideLibrary.nuclide_library_vector_[2][_col] = fsf;
                ModecNuclideLibrary.nuclide_library_vector_[3][_col] = q;
                ModecNuclideLibrary.nuclide_library_vector_[4][_col] = ampc;
                ModecNuclideLibrary.nuclide_library_vector_[5][_col] = wmpc;
            }

        }

        if (if_tracking_stockage == true) {
            TransMatrixStockage = TransMatrixDecay;
        }
        ////////////////////////

        // 截面矩阵matrix_XS建立

        //1: (n,g) cross-section of nuclide leading to ground state of the daughter nuclide
        if (SNG != 0.0 && _col != -1) {
            int productID = (NuclID + 10) / 10 * 10;
            _row = ModecNuclideLibrary.GetNuclIndex(productID);
            TransMatrixCrossSection.AddElement(_row, _col, SNG);
        }

        //2: (n,2n) cross-section of nuclide leading to ground state of the daughter nuclide
        if (SN2N != 0.0 && _col != -1) {
            int productID = (NuclID - 10) / 10 * 10;
            _row = ModecNuclideLibrary.GetNuclIndex(productID);
            TransMatrixCrossSection.AddElement(_row, _col, SN2N);
        }

        //3: (n,3n)/(n,a) cross-section of nuclide leading to ground state of the daughter nuclide
        if (SN3N_SNA != 0.0 && _col != -1) {
            if (NuclID > 900000) {	// actinide  (n,3n)
                int productID = (NuclID - 20) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(productID);
                TransMatrixCrossSection.AddElement(_row, _col, SN3N_SNA);
            } else {	// activation product and fission product (n,a)
                int productID = (NuclID - 20030) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(productID);
                TransMatrixCrossSection.AddElement(_row, _col, SN3N_SNA);
                TransMatrixCrossSection.AddElement(alphaID, _col, SN3N_SNA); // He-4
            }
        }

        //4: (n,f)/(n,p) cross-section of nuclide leading to ground state of the daughter nuclide
        if (SNF_SNP != 0.0 && _col != -1) {
            if (NuclID < 900000) {	//  activation product and fission product (n,p)
                int productID = (NuclID - 10000) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(productID);
                TransMatrixCrossSection.AddElement(_row, _col, SNF_SNP);
                TransMatrixCrossSection.AddElement(protonID, _col, SNF_SNP); // H-1
            }
        }

        //5: (n,gx) cross-section of nuclide leading to excited state of the daughter nuclide
        if (SNGX != 0.0 && _col != -1) {
            int productID = (NuclID + 10) / 10 * 10 + 1;
            _row = ModecNuclideLibrary.GetNuclIndex(productID);
            TransMatrixCrossSection.AddElement(_row, _col, SNGX);

        }

        //6: (n,2nx) cross-section of nuclide leading to excited state of the daughter nuclide
        if (SN2NX != 0.0 && _col != -1) {
            int productID = (NuclID - 10) / 10 * 10 + 1;
            _row = ModecNuclideLibrary.GetNuclIndex(productID);
            TransMatrixCrossSection.AddElement(_row, _col, SN2NX);
        }

        // Absorption/Reduction Rate of the nuclide
        if (totXS != 0.0 && _col != -1) {
            TransMatrixCrossSection.AddElement(_col, _col, -totXS);
            ModecNuclideLibrary.nuclide_library_vector_[6][_col] = totXS;
        }



        if (NuclID > 900000 && SNF_SNP != 0.0 && _col != -1) { // Depth 程序中认为核素ID大于900000认为其为锕系核素，且锕系核素有SNF但没有SNP
            ModecNuclideLibrary.nuclide_library_vector_[7][_col] = SNF_SNP;
        }
        if (_col != -1) {
            ModecNuclideLibrary.nuclide_library_vector_[8][_col] = SNG + SNGX;
        }


        ReadDepth.getline(line, 200);
    }
    ////////////////////////////////////////////////

    // read fission yield data from library
    while (!ReadDepth.eof()) {
        if (ReadDepth.peek() == '*' || ReadDepth.peek() == '!' || ReadDepth.peek() == '/') {
            ReadDepth.getline(line, 200);
        } else {
            break;
        }
    }

    ReadDepth >> FisNucNum;
    ModecNuclideLibrary.fissionable_nuclide_id_vector_.resize(FisNucNum);
    ModecNuclideLibrary.fission_products_id_vector_.resize(FisNucNum);
    ModecNuclideLibrary.effective_fission_yields_vector_.resize(FisNucNum);
    for (int i = 0; i < FisNucNum; ++i) {
        ReadDepth >> ModecNuclideLibrary.fissionable_nuclide_id_vector_[i]; // 将拥有裂变场的核素ID存储起来
    }

    for (int i = 0; i < FisNucNum; ++i) {
        ReadDepth >> itemp;
        ReadDepth >> PrdctNum;



        ModecNuclideLibrary.fission_products_id_vector_[i].resize(PrdctNum);
        ModecNuclideLibrary.effective_fission_yields_vector_[i].resize(PrdctNum);


        for (int j = 0; j < PrdctNum; ++j) {
            ReadDepth >> FPid;
            ModecNuclideLibrary.fission_products_id_vector_[i][j] = FPid;

            ReadDepth >> ModecNuclideLibrary.effective_fission_yields_vector_[i][j];
        }
        /// tag line ///
        ReadDepth >> itemp;

        if (itemp != -1) {
            InfoMessage::ErrorMessage("Position: void ModecClass::ReadFromDepthLib; \n Error: error flag in DEPTH lib!" ,1);
        }
        ReadDepth.getline(line, 200);
    }
    ////////////////////////////////////////////////

    ReadDepth.close();
    return;
};

void ModecClass::ReadFromDepthLibForTta() {
    char line[200];
    int NuclID, ProdID, itemp;
    double dtemp;
    int _col, _row;
    double lamda;
    int alphaID = ModecNuclideLibrary.GetNuclIndex(20040);
    int protonID = ModecNuclideLibrary.GetNuclIndex(10010);
    int FisNucNum, PrdctNum; // 有裂变场的裂变核素个数及裂变产物个数
    int FPid; // 裂变产物ID
    double FPyield; // 裂变产物份额
    //vector<int> FisNucIndex;

    double halfl, fb1, fp, fp1, fa, ft;       //// first line
    double fsf, fbn, fb, fbb, fneut, fba;    //// second line
    double q, abund, ampc, wmpc;			//// third line (only for ORIGEN-s decay library)
    double SNG, SN2N, SN3N_SNA, SNF_SNP, SNGX, SN2NX;       //// fourth line cross-section
    double fisrate;

    //int mat_dim = ModecNuclideLibrary.nuclide_library_vector_[0].size(); // 矩阵维度

    ifstream ReadDepth;
    ReadDepth.open(depth_library_name_, ios::in);
    if (!ReadDepth) {
        InfoMessage::ErrorMessage("Position: void ModecClass::ReadFromDepthLibForTta; \n Error: cannot open DEPTH library file! go and check the existence of the lib file.", 1);
    }
    while (!ReadDepth.eof()) {
        if (ReadDepth.peek() == '*' || ReadDepth.peek() == '!' || ReadDepth.peek() == '/') {
            ReadDepth.getline(line, 200);
        } else {
            break;
        }
    }

    // read decay and cross-section data from library
    while (!ReadDepth.eof()) {
        ReadDepth >> NuclID;
        if (NuclID == -1) {
            ReadDepth.getline(line, 200);
            break;
        }
        _col = ModecNuclideLibrary.GetNuclIndex(NuclID);
        if (_col == -1) {
            stringstream ss;
            string id;
            ss << NuclID;
            ss >> id;
            info_message_ = "Position: void ModecClass::ReadFromDepthLibForTta; \n Warning: unknown nuclide id " + id + ".";
            InfoMessage::ErrorMessage(info_message_, 0);
        }
        if (!(ReadDepth >> line)) {
            info_message_ = "Position: void ModecClass::ReadFromDepthLibForTta; \n Warning: unknown flag " + string(line) + ".";
            InfoMessage::ErrorMessage(info_message_, 0);
        }
        /// first line ///
        ReadDepth >> lamda >> fb1 >> fp >> fp1 >> fa >> ft;

        /// second line ///
        ReadDepth >> fsf >> fbn >> fb >> fbb >> fneut >> fba;

        /// third line ///
        ReadDepth >> q >> abund >> ampc >> wmpc >> dtemp >> dtemp;

        /// fourth line ///
        ReadDepth >> SNG >> SN2N >> SN3N_SNA >> SNF_SNP >> SNGX >> SN2NX;
        double totXS = SNG + SN2N + SN3N_SNA + SNF_SNP + SNGX + SN2NX;

        /// tag line ///
        ReadDepth >> itemp;

        if (itemp != -1) {
            InfoMessage::ErrorMessage("Position: void ModecClass::ReadFromDepthLibForTta; \n Error: error flag in DEPTH lib!", 1);
        }

        // 衰变矩阵TransMatrixDecay建立
        if (lamda != 0.0 && _col != -1) {
            //lamda = log(2) / halfl;

            if (fb1 != 0.0) {
                ProdID = (NuclID + 10000) / 10 * 10 + 1;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fb1;
                TtaMatrixDecay.AddElementCCS(_row, _col, _val);
            }

            if (fp != 0.0) {
                ProdID = (NuclID - 10000) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fp;
                TtaMatrixDecay.AddElementCCS(_row, _col, _val);
            }

            if (fp1 != 0.0) {
                ProdID = (NuclID - 10000) / 10 * 10 + 1;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fp1;
                TtaMatrixDecay.AddElementCCS(_row, _col, _val);
            }

            if (fa != 0.0) {
                ProdID = (NuclID - 20040) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fa;
                TtaMatrixDecay.AddElementCCS(_row, _col, _val);

                TtaMatrixDecay.AddElementCCS(alphaID, _col, _val); // 阿尔法衰变产生的He-4元素
            }

            if (ft != 0.0) {
                ProdID = NuclID - 1;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*ft;
                TtaMatrixDecay.AddElementCCS(_row, _col, _val);
            }

            if (fbn != 0.0) {
                ProdID = (NuclID + 10000 - 10) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fbn;
                TtaMatrixDecay.AddElementCCS(_row, _col, _val);
            }

            if (fb != 0.0) {
                ProdID = (NuclID + 10000) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fb;
                TtaMatrixDecay.AddElementCCS(_row, _col, _val);
            }

            if (fbb != 0.0) {
                ProdID = (NuclID + 20000) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fbb;
                TtaMatrixDecay.AddElementCCS(_row, _col, _val);
            }

            if (fneut != 0.0) {
                ProdID = (NuclID - 10) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fneut;
                TtaMatrixDecay.AddElementCCS(_row, _col, _val);
            }

            if (fba != 0.0) {
                ProdID = (NuclID - 10040) / 10 * 10;
                if (ProdID > 992550)
                    _row = ModecNuclideLibrary.GetNuclIndex(ProdID);
                double _val = lamda*fba;
                TtaMatrixDecay.AddElementCCS(_row, _col, _val);
                TtaMatrixDecay.AddElementCCS(alphaID, _col, _val); // 阿尔法衰变产生的He-4元素
            }

            if (_col != -1) {
                TtaMatrixDecay.AddElementCCS(_col, _col, -lamda);
                ModecNuclideLibrary.nuclide_library_vector_[1][_col] = lamda;
                ModecNuclideLibrary.nuclide_library_vector_[2][_col] = fsf;
                ModecNuclideLibrary.nuclide_library_vector_[3][_col] = q;
                ModecNuclideLibrary.nuclide_library_vector_[4][_col] = ampc;
                ModecNuclideLibrary.nuclide_library_vector_[5][_col] = wmpc;
            }

        }

        /*if (if_tracking_stockage == true)
        {
        	TransMatrixStockage = TransMatrixDecay;
        }*/
        ////////////////////////

        // 截面矩阵matrix_XS建立

        //1: (n,g) cross-section of nuclide leading to ground state of the daughter nuclide
        if (SNG != 0.0 && _col != -1) {
            int productID = (NuclID + 10) / 10 * 10;
            _row = ModecNuclideLibrary.GetNuclIndex(productID);
            TtaMatrixCrossSection.AddElementCCS(_row, _col, SNG);
        }

        //2: (n,2n) cross-section of nuclide leading to ground state of the daughter nuclide
        if (SN2N != 0.0 && _col != -1) {
            int productID = (NuclID - 10) / 10 * 10;
            _row = ModecNuclideLibrary.GetNuclIndex(productID);
            TtaMatrixCrossSection.AddElementCCS(_row, _col, SN2N);
        }

        //3: (n,3n)/(n,a) cross-section of nuclide leading to ground state of the daughter nuclide
        if (SN3N_SNA != 0.0 && _col != -1) {
            if (NuclID > 900000) {	// actinide  (n,3n)
                int productID = (NuclID - 20) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(productID);
                TtaMatrixCrossSection.AddElementCCS(_row, _col, SN3N_SNA);
            } else {	// activation product and fission product (n,a)
                int productID = (NuclID - 20030) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(productID);
                TtaMatrixCrossSection.AddElementCCS(_row, _col, SN3N_SNA);
                TtaMatrixCrossSection.AddElementCCS(alphaID, _col, SN3N_SNA); // He-4
            }
        }

        //4: (n,f)/(n,p) cross-section of nuclide leading to ground state of the daughter nuclide
        if (SNF_SNP != 0.0 && _col != -1) {
            if (NuclID < 900000) {	//  activation product and fission product (n,p)
                int productID = (NuclID - 10000) / 10 * 10;
                _row = ModecNuclideLibrary.GetNuclIndex(productID);
                TtaMatrixCrossSection.AddElementCCS(_row, _col, SNF_SNP);
                TtaMatrixCrossSection.AddElementCCS(protonID, _col, SNF_SNP); // H-1
            }
        }

        //5: (n,gx) cross-section of nuclide leading to excited state of the daughter nuclide
        if (SNGX != 0.0 && _col != -1) {
            int productID = (NuclID + 10) / 10 * 10 + 1;
            _row = ModecNuclideLibrary.GetNuclIndex(productID);
            TtaMatrixCrossSection.AddElementCCS(_row, _col, SNGX);

        }

        //6: (n,2nx) cross-section of nuclide leading to excited state of the daughter nuclide
        if (SN2NX != 0.0 && _col != -1) {
            int productID = (NuclID - 10) / 10 * 10 + 1;
            _row = ModecNuclideLibrary.GetNuclIndex(productID);
            TtaMatrixCrossSection.AddElementCCS(_row, _col, SN2NX);
        }

        // Absorption/Reduction Rate of the nuclide
        if (totXS != 0.0 && _col != -1) {
            TtaMatrixCrossSection.AddElementCCS(_col, _col, -totXS);
            ModecNuclideLibrary.nuclide_library_vector_[6][_col] = totXS;
        }



        if (NuclID > 900000 && SNF_SNP != 0.0 && _col != -1) { // Depth 程序中认为核素ID大于900000认为其为锕系核素，且锕系核素有SNF但没有SNP
            ModecNuclideLibrary.nuclide_library_vector_[7][_col] = SNF_SNP;
        }
        if (_col != -1) {
            ModecNuclideLibrary.nuclide_library_vector_[8][_col] = SNG + SNGX;
        }


        ReadDepth.getline(line, 200);
    }
    ////////////////////////////////////////////////

    // read fission yield data from library
    while (!ReadDepth.eof()) {
        if (ReadDepth.peek() == '*' || ReadDepth.peek() == '!' || ReadDepth.peek() == '/') {
            ReadDepth.getline(line, 200);
        } else {
            break;
        }
    }

    ReadDepth >> FisNucNum;
    ModecNuclideLibrary.fissionable_nuclide_id_vector_.resize(FisNucNum);
    ModecNuclideLibrary.fission_products_id_vector_.resize(FisNucNum);
    ModecNuclideLibrary.effective_fission_yields_vector_.resize(FisNucNum);
    for (int i = 0; i < FisNucNum; ++i) {
        ReadDepth >> ModecNuclideLibrary.fissionable_nuclide_id_vector_[i]; // 将拥有裂变场的核素ID存储起来
    }

    for (int i = 0; i < FisNucNum; ++i) {
        ReadDepth >> itemp;
        ReadDepth >> PrdctNum;



        ModecNuclideLibrary.fission_products_id_vector_[i].resize(PrdctNum);
        ModecNuclideLibrary.effective_fission_yields_vector_[i].resize(PrdctNum);


        for (int j = 0; j < PrdctNum; ++j) {
            ReadDepth >> FPid;
            ModecNuclideLibrary.fission_products_id_vector_[i][j] = FPid;

            ReadDepth >> ModecNuclideLibrary.effective_fission_yields_vector_[i][j];
        }
        /// tag line ///
        ReadDepth >> itemp;

        if (itemp != -1) {
            InfoMessage::ErrorMessage("Position: void ModecClass::ReadFromDepthLibForTta; \n Error: error flag in DEPTH lib!", 1);
        }
        ReadDepth.getline(line, 200);
    }
    ////////////////////////////////////////////////

    ReadDepth.close();
    return;
};

// 在ORIGENS或者TRITON耦合中，不需要进行裂变场的修正，只有在读取ORIGEN2的截面文件中，才需要进行修正
void ModecClass::ConstructFissionYieldsSpMat() {
    TransMatrixFissionYields.Reset();

    int MaxFisNucIndex;

    double FisRateSum = 0, MaxFisRate = 0;
    int FistAc = ModecNuclideLibrary.GetNuclIndex(902260);
    int LastAc = ModecNuclideLibrary.GetNuclIndex(992550);
    for (int i = FistAc; i <= LastAc; ++i) {
        int FisNucID = ModecNuclideLibrary.nuclide_list_[i];
        bool HasYield = false;
        for (unsigned int j = 0; j < ModecNuclideLibrary.fissionable_nuclide_id_vector_.size(); ++j) {
            if (FisNucID == ModecNuclideLibrary.fissionable_nuclide_id_vector_[j]) {
                HasYield = true;
                break;
            }
        }
        if (HasYield) {
            continue;
        }
        double FisXs = ModecNuclideLibrary.nuclide_library_vector_[7][i];
        FisRateSum = FisRateSum + ModecNuclideLibrary.nuclide_library_vector_[0][i] * FisXs;
        if (ModecNuclideLibrary.nuclide_library_vector_[0][i] * FisXs > MaxFisRate) {
            MaxFisRate = ModecNuclideLibrary.nuclide_library_vector_[0][i] * FisXs;
            MaxFisNucIndex = FisNucID;
        }
    }
    if (FisRateSum == 0) {
        yield_factor_ = 1.0;
        nearest_neighbor_ = 0;
    } else {
        int i, Distant;
        Distant = abs(MaxFisNucIndex - ModecNuclideLibrary.fissionable_nuclide_id_vector_[0]);
        for (i = 1; i < ModecNuclideLibrary.fissionable_nuclide_id_vector_.size(); ++i) {
            if (abs(MaxFisNucIndex - ModecNuclideLibrary.fissionable_nuclide_id_vector_[i]) > Distant) {
                break;
            }
            Distant = abs(MaxFisNucIndex - ModecNuclideLibrary.fissionable_nuclide_id_vector_[i]);
        }

        nearest_neighbor_ = ModecNuclideLibrary.fissionable_nuclide_id_vector_[i - 1];

        int NeiIndex = ModecNuclideLibrary.GetNuclIndex(nearest_neighbor_);

        if (ModecNuclideLibrary.nuclide_library_vector_[0][NeiIndex] == 0.0) {
            yield_factor_ = 1.0;
            nearest_neighbor_ = 0;
        } else {
            //int NeighborID = ModecNuclideLibrary.GetNuclIndex(nearest_neighbor_);
            yield_factor_ = 1.0 + FisRateSum / (ModecNuclideLibrary.nuclide_library_vector_[0][NeiIndex] * ModecNuclideLibrary.nuclide_library_vector_[7][NeiIndex]);
        }
    }


    int _col, _row;
    int ReacID, ProdID;
    double FisRate, ReactRate;
    for (unsigned int i = 0; i < ModecNuclideLibrary.fissionable_nuclide_id_vector_.size(); ++i) {
        ReacID = ModecNuclideLibrary.fissionable_nuclide_id_vector_[i];
        _col = ModecNuclideLibrary.GetNuclIndex(ReacID);
        FisRate = ModecNuclideLibrary.nuclide_library_vector_[7][_col];
        if (nearest_neighbor_ == ReacID) {
            FisRate = FisRate*yield_factor_;
        }

        for (unsigned int j = 0; j < ModecNuclideLibrary.fission_products_id_vector_[i].size(); ++j) {
            int ProdId = ModecNuclideLibrary.fission_products_id_vector_[i][j];
            if (ProdId < 10010 || ProdId >= 1050000) {
                continue;
            }
            _row = ModecNuclideLibrary.GetNuclIndex(ProdId);
            ReactRate = FisRate * ModecNuclideLibrary.effective_fission_yields_vector_[i][j];
            if (ReactRate != 0) {
                TransMatrixFissionYields.AddElement(_row, _col, ReactRate);
            }

        }
    }
    return;
}

void ModecClass::ConstructFissionYieldsSpMatForTta() {
    TtaMatrixFissionYields.Reset();

    int MaxFisNucIndex;

    double FisRateSum = 0, MaxFisRate = 0;
    int FistAc = ModecNuclideLibrary.GetNuclIndex(902260);
    int LastAc = ModecNuclideLibrary.GetNuclIndex(992550);
    for (int i = FistAc; i <= LastAc; ++i) {
        int FisNucID = ModecNuclideLibrary.nuclide_list_[i];
        bool HasYield = false;
        for (unsigned int j = 0; j < ModecNuclideLibrary.fissionable_nuclide_id_vector_.size(); ++j) {
            if (FisNucID == ModecNuclideLibrary.fissionable_nuclide_id_vector_[j]) {
                HasYield = true;
                break;
            }
        }
        if (HasYield) {
            continue;
        }
        double FisXs = ModecNuclideLibrary.nuclide_library_vector_[7][i];
        FisRateSum = FisRateSum + ModecNuclideLibrary.nuclide_library_vector_[0][i] * FisXs;
        if (ModecNuclideLibrary.nuclide_library_vector_[0][i] * FisXs > MaxFisRate) {
            MaxFisRate = ModecNuclideLibrary.nuclide_library_vector_[0][i] * FisXs;
            MaxFisNucIndex = FisNucID;
        }
    }
    if (FisRateSum == 0) {
        yield_factor_ = 1.0;
        nearest_neighbor_ = 0;
    } else {
        int i, Distant;
        Distant = abs(MaxFisNucIndex - ModecNuclideLibrary.fissionable_nuclide_id_vector_[0]);
        for (i = 1; i < ModecNuclideLibrary.fissionable_nuclide_id_vector_.size(); ++i) {
            if (abs(MaxFisNucIndex - ModecNuclideLibrary.fissionable_nuclide_id_vector_[i]) > Distant) {
                break;
            }
            Distant = abs(MaxFisNucIndex - ModecNuclideLibrary.fissionable_nuclide_id_vector_[i]);
        }

        nearest_neighbor_ = ModecNuclideLibrary.fissionable_nuclide_id_vector_[i - 1];

        int NeiIndex = ModecNuclideLibrary.GetNuclIndex(nearest_neighbor_);

        if (ModecNuclideLibrary.nuclide_library_vector_[0][NeiIndex] == 0.0) {
            yield_factor_ = 1.0;
            nearest_neighbor_ = 0;
        } else {
            //int NeighborID = ModecNuclideLibrary.GetNuclIndex(nearest_neighbor_);
            yield_factor_ = 1.0 + FisRateSum / (ModecNuclideLibrary.nuclide_library_vector_[0][NeiIndex] * ModecNuclideLibrary.nuclide_library_vector_[7][NeiIndex]);
        }
    }


    int _col, _row;
    int ReacID, ProdID;
    double FisRate, ReactRate;
    for (unsigned int i = 0; i < ModecNuclideLibrary.fissionable_nuclide_id_vector_.size(); ++i) {
        ReacID = ModecNuclideLibrary.fissionable_nuclide_id_vector_[i];
        _col = ModecNuclideLibrary.GetNuclIndex(ReacID);
        FisRate = ModecNuclideLibrary.nuclide_library_vector_[7][_col];
        if (nearest_neighbor_ == ReacID) {
            FisRate = FisRate*yield_factor_;
        }

        for (unsigned int j = 0; j < ModecNuclideLibrary.fission_products_id_vector_[i].size(); ++j) {
            int ProdId = ModecNuclideLibrary.fission_products_id_vector_[i][j];
            if (ProdId < 10010 || ProdId >= 1050000) {
                continue;
            }
            _row = ModecNuclideLibrary.GetNuclIndex(ProdId);
            ReactRate = FisRate * ModecNuclideLibrary.effective_fission_yields_vector_[i][j];
            if (ReactRate != 0) {
                TtaMatrixFissionYields.AddElementCCS(_row, _col, ReactRate);
            }

        }
    }
    return;
}
//////////////////////////////////////////////////////////////////////////////////////////////////////
