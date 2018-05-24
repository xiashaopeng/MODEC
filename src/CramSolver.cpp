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

// 采用拉普拉斯变换方法，来计算添料率的燃耗方程
void SolveTrans::PfdCramSolver(SpMat &matrix, vector <double> &N, vector <double> &F, const double &time) {
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
    N_eol.resize(dim);

	
    for (int j = 0; j < 8; ++j) {
		Complex temp = time / theta[j];
        for (unsigned int i = 0; i < dim; ++i) {
            N_temp[i] = (N[i] + F[i] * temp)* alpha[j];
        }
        matrix.LUEliminationForCram(theta[j], time, N_temp, dim);
        for (unsigned int i = 0; i < dim; ++i) {
            N_eol[i] += 2.0 * N_temp[i];
        }
    }
    for (int i = 0; i < dim; ++i) {
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

void SolveTrans::PfdCramSolver(SpMat &matrix, const SpMat &TransMatrixReprocess, SpMat &TransMatrixStockage, vector < double > &N, vector <double> &F, const double &time) {
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
	vector<Complex > N_core(dim / 2);
	//vector <Complex > N_temp(dim / 2);
	vector <Complex > N_stockage(dim / 2);
	vector <Complex > N_eol;
	N_eol.resize(dim);

	//vector <double> N_eol_real;
	//N_eol_real.resize(dim);

	for (int j = 0; j < 8; ++j) {
		Complex temp = time / theta[j];
		
		for (unsigned int i = 0; i < dim / 2; ++i) {
			N_core[i] = ( N[i] + F[i] * temp) * alpha[j];
			N_stockage[i] = ( N[i + dim / 2] + F[i + dim / 2] * temp ) * alpha[j];
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

}


void SolveTrans::IpfCramSolver32(SpMat &matrix, vector <double> &N, const double &time) {
	const double alpha0 = 6.932444346272945E-32;
	const Complex alpha[16] = {
		Complex(+2.531738087291248E+03, -2.175335564251554E+05),
		Complex(+1.124951925994460E+02, -2.293141996821969E+02),
		Complex(+1.928222545500035E+02, -5.268292754604315E+02),
		Complex(+4.159057536149641E+04, -3.509499779111824E+05),
		Complex(+4.981965659174993E+03, -5.519940772045004E+05),
		Complex(+1.498320382271818E+02, -6.815792364464349E+02),
		Complex(+6.462582753425457E+02, -6.286003508366936E+03),
		Complex(+3.108473530705140E+02, -2.227478875866931E+04),
		Complex(+2.885604705807475E+02, -1.308904072042900E+03),
		Complex(+7.047430374921731E+01, -3.234996756134937E+02),
		Complex(+7.344336306595115E+01, -2.581272076901578E+02),
		Complex(+4.519923142224831E+01, -4.759858266475396E+00),
		Complex(+5.271307338359870E+01, -2.583918268027449E+01),
		Complex(+4.983259420279240E+01, -1.542812979380035E+01),
		Complex(+6.057642473653785E+01, -6.941674585178426E+02),
		Complex(+4.534413157137819E+01, -3.114070953092643E+01)
	};

	const Complex theta[16] = {
		Complex(+9.093582328296485E+00, +1.321250259877626E+01),
		Complex(-1.727168130828052E+00, +2.572302624565582E+01),
		Complex(+1.166498099877330E+00, +2.314228200408637E+01),
		Complex(+5.780372991371440E+00, +1.811790953587012E+01),
		Complex(+1.201313618215791E+01, +5.977830703724750E+00),
		Complex(+7.585235004412852E+00, +1.565385209991484E+01),
		Complex(+3.652646888760184E+00, +2.061118862125454E+01),
		Complex(+1.272430588935007E+01, +1.194256446871873E+00),
		Complex(-5.098669890503975E+00, +2.837060507799692E+01),
		Complex(+1.129472112530614E+01, +8.378353071561875E+00),
		Complex(+1.032515517992622E+01, +1.078873149419210E+01),
		Complex(-2.734798356101249E+01, +4.069231229671756E+01),
		Complex(-1.377054550962031E+01, +3.399214933428828E+01),
		Complex(-1.958248015332671E+01, +3.710329224494451E+01),
		Complex(+1.248804035352348E+01, +3.584068696364335E+00),
		Complex(-9.054506283982660E+00, +3.111173004087709E+01)
	};

    int dim = N.size();
    vector <Complex > N_temp(dim);
    vector <Complex > N_eol;
    //vector<double> N_eol_real;
    N_eol.resize(dim);
    //N_eol_real.resize(dim);
	
	for (unsigned int i = 0; i < dim; ++i) {
		N_eol[i] = N[i];
	}
    for (int j = 0; j < 16; ++j) {
        N_temp = matrix.LUEliminationForIpfCram(theta[j], time, N_eol);
        for (unsigned int i = 0; i < dim; ++i) {
            N_eol[i] += 2.0 * (alpha[j] * N_temp[i]).real_;
        }
    }
    for (unsigned int i = 0; i < dim; ++i) {
        N[i] = N_eol[i].real_ * alpha0;
    }
}

void SolveTrans::IpfCramSolver48(SpMat &matrix, vector <double> &N, const double &time) {
	const double alpha0 = 2.258038182743983E-47;
	const Complex alpha[24] = {
		Complex(+6.387380733878774E+02, -6.743912502859256E+02),
		Complex(+1.909896179065730E+02, -3.973203432721332E+02),
		Complex(+4.236195226571914E+02, -2.041233768918671E+03),
		Complex(+4.645770595258726E+02, -1.652917287299683E+03),
		Complex(+7.765163276752433E+02, -1.783617639907328E+04),
		Complex(+1.907115136768522E+03, -5.887068595142284E+04),
		Complex(+2.909892685603256E+03, -9.953255345514560E+03),
		Complex(+1.944772206620450E+02, -1.427131226068449E+03),
		Complex(+1.382799786972332E+05, -3.256885197214938E+06),
		Complex(+5.628442079602433E+03, -2.924284515884309E+04),
		Complex(+2.151681283794220E+02, -1.121774011188224E+03),
		Complex(+1.324720240514420E+03, -6.370088443140973E+04),
		Complex(+1.617548476343347E+04, -1.008798413156542E+06),
		Complex(+1.112729040439685E+02, -8.837109731680418E+01),
		Complex(+1.074624783191125E+02, -1.457246116408180E+02),
		Complex(+8.835727765158191E+01, -6.388286188419360E+01),
		Complex(+9.354078136054179E+01, -2.195424319460237E+02),
		Complex(+9.418142823531573E+01, -6.719055740098035E+02),
		Complex(+1.040012390717851E+02, -1.693747595553868E+02),
		Complex(+6.861882624343235E+01, -1.177598523430493E+01),
		Complex(+8.766654491283722E+01, -4.596464999363902E+03),
		Complex(+1.056007619389650E+02, -1.738294585524067E+03),
		Complex(+7.738987569039419E+01, -4.311715386228984E+01),
		Complex(+1.041366366475571E+02, -2.777743732451969E+02)
	};

	const Complex theta[24] = {
		Complex(-4.465731934165702E+01, +6.233225190695437E+01),
		Complex(-5.284616241568964E+00, +4.057499381311059E+01),
		Complex(-8.867715667624458E+00, +4.325515754166724E+01),
		Complex(+3.493013124279215E+00, +3.281615453173585E+01),
		Complex(+1.564102508858634E+01, +1.558061616372237E+01),
		Complex(+1.742097597385893E+01, +1.076629305714420E+01),
		Complex(-2.834466755180654E+01, +5.492841024648724E+01),
		Complex(+1.661569367939544E+01, +1.316994930024688E+01),
		Complex(+8.011836167974721E+00, +2.780232111309410E+01),
		Complex(-2.056267541998229E+00, +3.794824788914354E+01),
		Complex(+1.449208170441839E+01, +1.799988210051809E+01),
		Complex(+1.853807176907916E+01, +5.974332563100539E+00),
		Complex(+9.932562704505182E+00, +2.532823409972962E+01),
		Complex(-2.244223871767187E+01, +5.179633600312162E+01),
		Complex(+8.590014121680897E-01, +3.536456194294350E+01),
		Complex(-1.286192925744479E+01, +4.600304902833652E+01),
		Complex(+1.164596909542055E+01, +2.287153304140217E+01),
		Complex(+1.806076684783089E+01, +8.368200580099821E+00),
		Complex(+5.870672154659249E+00, +3.029700159040121E+01),
		Complex(-3.542938819659747E+01, +5.834381701800013E+01),
		Complex(+1.901323489060250E+01, +1.194282058271408E+00),
		Complex(+1.885508331552577E+01, +3.583428564427879E+00),
		Complex(-1.734689708174982E+01, +4.883941101108207E+01),
		Complex(+1.316284237125190E+01, +2.042951874827759E+01)
	};

    int dim = N.size();
    vector <Complex > N_temp(dim);
    vector <Complex > N_eol;
    //vector<double> N_eol_real;
    N_eol.resize(dim);
    //N_eol_real.resize(dim);

	for (unsigned int i = 0; i < dim; ++i) {
		N_eol[i] = N[i];
	}
    for (int j = 0; j < 24; ++j) {
        N_temp = matrix.LUEliminationForIpfCram(theta[j], time, N_eol);
        for (unsigned int i = 0; i < dim; ++i) {
            N_eol[i] += 2.0 * (alpha[j] * N_temp[i]).real_;
        }
    }
    for (unsigned int i = 0; i < dim; ++i) {
        //N_eol[i] += N[i] * alpha0;
        N[i] = N_eol[i].real_ * alpha0;
    }

}


void SolveTrans::IpfCramSolver32(SpMat &matrix, const SpMat &TransMatrixReprocess, SpMat &TransMatrixStockage, vector < double > &N, const double &time) {
	const double alpha0 = 6.932444346272945E-32;
	const Complex alpha[16] = {
		Complex(+2.531738087291248E+03, -2.175335564251554E+05),
		Complex(+1.124951925994460E+02, -2.293141996821969E+02),
		Complex(+1.928222545500035E+02, -5.268292754604315E+02),
		Complex(+4.159057536149641E+04, -3.509499779111824E+05),
		Complex(+4.981965659174993E+03, -5.519940772045004E+05),
		Complex(+1.498320382271818E+02, -6.815792364464349E+02),
		Complex(+6.462582753425457E+02, -6.286003508366936E+03),
		Complex(+3.108473530705140E+02, -2.227478875866931E+04),
		Complex(+2.885604705807475E+02, -1.308904072042900E+03),
		Complex(+7.047430374921731E+01, -3.234996756134937E+02),
		Complex(+7.344336306595115E+01, -2.581272076901578E+02),
		Complex(+4.519923142224831E+01, -4.759858266475396E+00),
		Complex(+5.271307338359870E+01, -2.583918268027449E+01),
		Complex(+4.983259420279240E+01, -1.542812979380035E+01),
		Complex(+6.057642473653785E+01, -6.941674585178426E+02),
		Complex(+4.534413157137819E+01, -3.114070953092643E+01)
	};

	const Complex theta[16] = {
		Complex(+9.093582328296485E+00, +1.321250259877626E+01),
		Complex(-1.727168130828052E+00, +2.572302624565582E+01),
		Complex(+1.166498099877330E+00, +2.314228200408637E+01),
		Complex(+5.780372991371440E+00, +1.811790953587012E+01),
		Complex(+1.201313618215791E+01, +5.977830703724750E+00),
		Complex(+7.585235004412852E+00, +1.565385209991484E+01),
		Complex(+3.652646888760184E+00, +2.061118862125454E+01),
		Complex(+1.272430588935007E+01, +1.194256446871873E+00),
		Complex(-5.098669890503975E+00, +2.837060507799692E+01),
		Complex(+1.129472112530614E+01, +8.378353071561875E+00),
		Complex(+1.032515517992622E+01, +1.078873149419210E+01),
		Complex(-2.734798356101249E+01, +4.069231229671756E+01),
		Complex(-1.377054550962031E+01, +3.399214933428828E+01),
		Complex(-1.958248015332671E+01, +3.710329224494451E+01),
		Complex(+1.248804035352348E+01, +3.584068696364335E+00),
		Complex(-9.054506283982660E+00, +3.111173004087709E+01)
	};

    int dim = N.size();

    // 判断dim是否为奇数，若为奇数，意味着采用了利用增广矩阵方法求解常添料率的方法
    // 若为偶数，则采用数值积分方法求解常添料率燃耗方程
    if (dim % 2 == 0) {
        vector <Complex> N_core(dim / 2);
        vector <Complex> N_stockage(dim / 2);
        vector <Complex> N_eol;
        N_eol.resize(dim);
		
		vector <Complex> N_temp_core(dim / 2);
		vector <Complex> N_temp_stockage(dim / 2);

		vector <Complex> N_stockage_c(dim / 2);
		
		for (unsigned int i = 0; i < dim / 2; ++i) {
			N_core[i] = N[i];
			N_stockage[i] = N[i + dim / 2];
		}
        for (int j = 0; j < 16; ++j) {

            N_temp_core = matrix.LUEliminationForIpfCram(theta[j], time, N_core);

            for (unsigned int i = 0; i < dim / 2; ++i) {
                N_core[i] += 2.0 * (alpha[j] * N_temp_core[i]).real_;
            }

            for (int ii = 0; ii < dim / 2; ++ii) {
				N_stockage_c[ii] = N_stockage[ii] - TransMatrixReprocess.diagonal_val_[ii] * time * N_temp_core[ii];
            }
            N_temp_stockage = TransMatrixStockage.LUEliminationForIpfCram(theta[j], time, N_stockage_c);

            for (unsigned int i = 0; i < dim / 2; ++i) {
                N_stockage[i] += 2.0 * (alpha[j] * N_temp_stockage[i]).real_;
            }

        }

		for (unsigned int i = 0; i < dim; ++i) {
			if (i < dim / 2) {
				N[i] = N_core[i].real_ * alpha0;
			}
			else {
				N[i] = N_stockage[i - dim/2].real_ * alpha0;
			}
		}

        //return N_eol_real;
    } else {
        vector< Complex > N_core((dim + 1) / 2); // 堆芯核素浓度向量多一个表示添料率的伪核素，其数值等于1
        //vector <Complex > N_temp(dim / 2);
        vector < Complex > N_stockage((dim - 1) / 2);
        vector < Complex > N_eol;
        N_eol.resize(dim);
        //N_eol_real.resize(dim / 2);

        vector <Complex > N_temp_core((dim + 1) / 2);
		vector <Complex > N_temp_stockage((dim - 1) / 2);

		vector <Complex > N_stockage_c((dim - 1) / 2);
		
		for (int i = 0; i < (dim - 1) / 2; ++i) {
			N_core[i] = N[i];
			N_stockage[i] = N[i + (dim + 1) / 2];
		}
		N_core[(dim - 1) / 2] = N[(dim - 1) / 2];

        for (int j = 0; j < 16; ++j) {
            
            N_temp_core = matrix.LUEliminationForIpfCram(theta[j], time, N_core);

            for (int i = 0; i < (dim + 1) / 2; ++i) {
                N_core[i] += 2.0 * (alpha[j] * N_temp_core[i]).real_;
            }

            for (int ii = 0; ii < (dim - 1) / 2; ++ii) {
				N_stockage_c[ii] = N_stockage[ii] - TransMatrixReprocess.diagonal_val_[ii] * time * N_temp_core[ii];
            }
            N_temp_stockage = TransMatrixStockage.LUEliminationForIpfCram(theta[j], time, N_stockage_c);

            for (unsigned int i = 0; i < (dim - 1) / 2; ++i) {
                N_stockage[i] += 2.0 * (alpha[j] * N_temp_stockage[i]).real_;
            }

        }
        for (unsigned int i = 0; i < dim; ++i) {
            if (i < (dim + 1)/2) {
				N[i] = N_core[i].real_ * alpha0;
			}
			else {
				N[i] = N_stockage[i - (dim + 1) / 2].real_ * alpha0;
			}           
        }

        //return N_eol_real;
    }

}
