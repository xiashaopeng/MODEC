//#include <iostream>
#include "SolveTrans.h"

void SolveTrans::PfdCramSolver(SpMat &matrix, vector <double> &N,const double &time) {
    const double alpha0 = 2.1248537104952237488E-16;
    const Complex alpha[8] = {
        Complex(-5.0901521865224915650E-07, -2.4220017652852287970E-05),
        Complex(+2.1151742182466030907E-04, +4.3892969647380673918E-03),
        Complex(+1.1339775178483930527E+02, +1.0194721704215856450E+02),
        Complex(+1.5059585270023467528E+01, -5.7514052776421819979E+00),
        Complex(-6.4500878025539646595E+01, -2.2459440762652096056E+02),
        Complex(-1.4793007113557999718E+00, +1.7686588323782937906E+00),
        Complex(-6.2518392463207918892E+01, -1.1190391094283228480E+01),
        Complex(+4.1023136835410021273E-02, -1.5743466173455468191E-01)
    };

    const Complex theta[8] = {
        Complex(-1.0843917078696988026E+01, +1.9277446167181652284E+01),
        Complex(-5.2649713434426468895E+00, +1.6220221473167927305E+01),
        Complex(+5.9481522689511774808E+00, +3.5874573620183222829E+00),
        Complex(+3.5091036084149180974E+00, +8.4361989858843750826E+00),
        Complex(+6.4161776990994341923E+00, +1.1941223933701386874E+00),
        Complex(+1.4193758971856659786E+00, +1.0925363484496722585E+01),
        Complex(+4.9931747377179963991E+00, +5.9968817136039422260E+00),
        Complex(-1.4139284624888862114E+00, +1.3497725698892745389E+01)
    };

    int dim = N.size();
    vector < Complex > N_temp(dim);
    vector < Complex > N_eol;
    //vector<double> N_eol_real;
    N_eol.resize(dim);
    //N_eol_real.resize(dim);

    for (int j = 0; j < 8; ++j) {
        for (unsigned int i = 0; i < dim; ++i) {
            N_temp[i] = N[i] * alpha[j];
        }
        matrix.LUEliminationForCram(theta[j], time, N_temp, dim);
        for (unsigned int i = 0; i < dim; ++i) {
            N_eol[i] += 2.0 * N_temp[i];
        }
    }
    for (int i = 0; i < dim; ++i) {
        //N_eol[i] += N[i] * alpha0;
        N[i] = N_eol[i].real_ + N[i] * alpha0;
    }

    //return N_eol_real;
}

void SolveTrans::PfdCramSolver(SpMat &matrix, const SpMat &TransMatrixReprocess, SpMat &TransMatrixStockage, vector < double > &N, const double &time) {
    const double alpha0 = 2.1248537104952237488E-16;
    const Complex alpha[8] = {
        Complex(-5.0901521865224915650E-07, -2.4220017652852287970E-05),
        Complex(+2.1151742182466030907E-04, +4.3892969647380673918E-03),
        Complex(+1.1339775178483930527E+02, +1.0194721704215856450E+02),
        Complex(+1.5059585270023467528E+01, -5.7514052776421819979E+00),
        Complex(-6.4500878025539646595E+01, -2.2459440762652096056E+02),
        Complex(-1.4793007113557999718E+00, +1.7686588323782937906E+00),
        Complex(-6.2518392463207918892E+01, -1.1190391094283228480E+01),
        Complex(+4.1023136835410021273E-02, -1.5743466173455468191E-01)
    };

    const Complex theta[8] = {
        Complex(-1.0843917078696988026E+01, +1.9277446167181652284E+01),
        Complex(-5.2649713434426468895E+00, +1.6220221473167927305E+01),
        Complex(+5.9481522689511774808E+00, +3.5874573620183222829E+00),
        Complex(+3.5091036084149180974E+00, +8.4361989858843750826E+00),
        Complex(+6.4161776990994341923E+00, +1.1941223933701386874E+00),
        Complex(+1.4193758971856659786E+00, +1.0925363484496722585E+01),
        Complex(+4.9931747377179963991E+00, +5.9968817136039422260E+00),
        Complex(-1.4139284624888862114E+00, +1.3497725698892745389E+01)
    };

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

        for (int j = 0; j < 8; ++j) {
            for (unsigned int i = 0; i < dim / 2; ++i) {
                N_core[i] = N[i] * alpha[j];
                N_stockage[i] = N[i + dim / 2] * alpha[j];
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
            N[i] = N_eol[i].real_ + N[i] * alpha0;
        }

        //return N_eol_real;
    } else {
        vector<Complex > N_core((dim + 1)/2); // 堆芯核素浓度向量多一个表示添料率的伪核素，其数值等于1
        //vector <Complex > N_temp(dim / 2);
        vector <Complex > N_stockage((dim - 1) / 2);
        vector <Complex > N_eol;
        N_eol.resize(dim);

        //vector <double> N_eol_real;
        //N_eol_real.resize(dim);


        for (int j = 0; j < 8; ++j) {
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
            N[i] = N_eol[i].real_ + N[i] * alpha0;
        }

        //return N_eol_real;
    }

}

void SolveTrans::IpfCramSolver(const int &order, SpMat &matrix, vector <double> &N, const double &time) {

    const double alpha0 = 2.1248537104952237488E-16;
    const Complex alpha[8] = {
        Complex(-5.0901521865224915650E-07, -2.4220017652852287970E-05),
        Complex(+2.1151742182466030907E-04, +4.3892969647380673918E-03),
        Complex(+1.1339775178483930527E+02, +1.0194721704215856450E+02),
        Complex(+1.5059585270023467528E+01, -5.7514052776421819979E+00),
        Complex(-6.4500878025539646595E+01, -2.2459440762652096056E+02),
        Complex(-1.4793007113557999718E+00, +1.7686588323782937906E+00),
        Complex(-6.2518392463207918892E+01, -1.1190391094283228480E+01),
        Complex(+4.1023136835410021273E-02, -1.5743466173455468191E-01)
    };

    const Complex theta[8] = {
        Complex(-1.0843917078696988026E+01, +1.9277446167181652284E+01),
        Complex(-5.2649713434426468895E+00, +1.6220221473167927305E+01),
        Complex(+5.9481522689511774808E+00, +3.5874573620183222829E+00),
        Complex(+3.5091036084149180974E+00, +8.4361989858843750826E+00),
        Complex(+6.4161776990994341923E+00, +1.1941223933701386874E+00),
        Complex(+1.4193758971856659786E+00, +1.0925363484496722585E+01),
        Complex(+4.9931747377179963991E+00, +5.9968817136039422260E+00),
        Complex(-1.4139284624888862114E+00, +1.3497725698892745389E+01)
    };

    int dim = N.size();
    vector <Complex > N_temp(dim);
    vector <Complex > N_eol;
    //vector<double> N_eol_real;
    N_eol.resize(dim);
    //N_eol_real.resize(dim);

    for (int j = 0; j < 8; ++j) {
        for (unsigned int i = 0; i < dim; ++i) {
            N_temp[i] = N[i] * alpha[j];
        }
        matrix.LUEliminationForCram(theta[j], time, N_temp, dim);
        for (unsigned int i = 0; i < dim; ++i) {
            N_eol[i] += 2.0 * N_temp[i];
        }
    }
    for (unsigned int i = 0; i < dim; ++i) {
        //N_eol[i] += N[i] * alpha0;
        N[i] = N_eol[i].real_ + N[i] * alpha0;
    }

    //return N_eol_real;
}

void SolveTrans::IpfCramSolver(const int &order, SpMat &matrix, const SpMat &TransMatrixReprocess, SpMat &TransMatrixStockage, vector < double > &N, const double &time) {
    const double alpha0 = 2.1248537104952237488E-16;
    const Complex alpha[8] = {
        Complex(-5.0901521865224915650E-07, -2.4220017652852287970E-05),
        Complex(+2.1151742182466030907E-04, +4.3892969647380673918E-03),
        Complex(+1.1339775178483930527E+02, +1.0194721704215856450E+02),
        Complex(+1.5059585270023467528E+01, -5.7514052776421819979E+00),
        Complex(-6.4500878025539646595E+01, -2.2459440762652096056E+02),
        Complex(-1.4793007113557999718E+00, +1.7686588323782937906E+00),
        Complex(-6.2518392463207918892E+01, -1.1190391094283228480E+01),
        Complex(+4.1023136835410021273E-02, -1.5743466173455468191E-01)
    };

    const Complex theta[8] = {
        Complex(-1.0843917078696988026E+01, +1.9277446167181652284E+01),
        Complex(-5.2649713434426468895E+00, +1.6220221473167927305E+01),
        Complex(+5.9481522689511774808E+00, +3.5874573620183222829E+00),
        Complex(+3.5091036084149180974E+00, +8.4361989858843750826E+00),
        Complex(+6.4161776990994341923E+00, +1.1941223933701386874E+00),
        Complex(+1.4193758971856659786E+00, +1.0925363484496722585E+01),
        Complex(+4.9931747377179963991E+00, +5.9968817136039422260E+00),
        Complex(-1.4139284624888862114E+00, +1.3497725698892745389E+01)
    };

    int dim = N.size();

    // 判断dim是否为奇数，若为奇数，意味着采用了利用增广矩阵方法求解常添料率的方法
    // 若为偶数，则采用数值积分方法求解常添料率燃耗方程
    if (dim % 2 == 0) {
        vector<Complex > N_core(dim / 2);
        //vector <Complex > N_temp(dim / 2);
        vector <Complex > N_stockage(dim / 2);
        vector <Complex > N_eol;
        N_eol.resize(dim);
        //N_eol_real.resize(dim / 2);
        //vector <double> N_eol_real;
        //N_eol_real.resize(dim);


        for (int j = 0; j < 8; ++j) {
            for (unsigned int i = 0; i < dim / 2; ++i) {
                N_core[i] = N[i] * alpha[j];
                N_stockage[i] = N[i + dim / 2] * alpha[j];
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
            N[i] = N_eol[i].real_ + N[i] * alpha0;
        }

        //return N_eol_real;
    } else {
        vector< Complex > N_core((dim + 1) / 2); // 堆芯核素浓度向量多一个表示添料率的伪核素，其数值等于1
        //vector <Complex > N_temp(dim / 2);
        vector < Complex > N_stockage((dim - 1) / 2);
        vector < Complex > N_eol;
        N_eol.resize(dim);
        //N_eol_real.resize(dim / 2);

        //vector <double> N_eol_real;
        //N_eol_real.resize(dim);


        for (int j = 0; j < 8; ++j) {
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
            N[i] = N_eol[i].real_ + N[i] * alpha0;
        }

        //return N_eol_real;
    }

}
