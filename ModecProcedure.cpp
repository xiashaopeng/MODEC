#include "ModecClass.h"


void ModecClass::ModecProcedure()
{
	int size = evolution_mode_.size();
	for (int i = 0; i < size; ++i)
	{
		int mode = evolution_mode_[i];
		if (mode == 1)
		{
			ModecNuclideLibrary.flux_ = evolution_value_[i];
		}
		if (mode == 2)
		{
			ModecNuclideLibrary.specified_power_ = evolution_value_[i];
		}
		double time = burnup_time_[i];
		int subtime = substep_[i];
		Evolution(mode + 3 * if_flow_mode_, time, subtime);
	}
};


//void ModecClass::CalEquilibrium(int mode)
//{
//	double time = 100000.0 * 24 * 3600;
//	double error_cutoff = 1.0e-7;
//
//	vector<double> error_;
//	error_.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size());
//
//	double error_max;
//
//	if (mode == 0) // ��˥�����
//	{
//		TransMatrixDecay.SymbolLUElimination();
//		flux_vector_.push_back(0.0);
//		power_vector_.push_back(0.0);
//		while(1)
//		{
//			ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrixDecay, ModecNuclideLibrary.nuclide_library_vector_[0], time);
//			if (if_constant_online_feeding_ == true)
//			{
//				int size_F = constant_feeding_vector_.size();
//
//				vector<double > F_mol;
//				F_mol.resize(size_F);
//
//				int size_GL = gauss_legendre_abscissa_.size();
//				for (int GL_i = 0; GL_i < size_GL; ++GL_i)
//				{
//					vector<double > F_temp;
//					F_temp.resize(size_F);
//					for (int F_i = 0; F_i < size_F; ++F_i)
//					{
//						F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
//					}
//					double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);
//
//					F_temp = Solver.PfdCramSolver(TransMatrixDecay, F_temp, time_gl);
//
//					for (int F_i = 0; F_i < size_F; ++F_i)
//					{
//						F_mol[F_i] += F_temp[F_i];
//					}
//				}
//
//				for (int F_i = 0; F_i < size_F; ++F_i)
//				{
//					ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
//				}
//			}
//
//			error_max = 0;
//			for (unsigned int i = 0; i < n_vector_.size(); ++i)
//			{
//				error_[i] = abs((ModecNuclideLibrary.nuclide_library_vector_[0][i] - n_vector_[1][i]) / n_vector_[1][i]);
//				if (error_[i] > error_max)
//				{
//					error_max = error_[i];
//				}
//			}
//
//			if (error_max < error_cutoff)
//			{
//				break;
//			}
//			else
//			{
//				n_vector_[1] = ModecNuclideLibrary.nuclide_library_vector_[0];
//			}
//		}
//		flux_vector_.push_back(0.0);
//		power_vector_.push_back(0.0);
//	}
//
//	if (mode == 1) // ��ͨ�����
//	{
//		if (lib_tag_ == 1) //lib_tag_ = 1��ζ�Ŷ�ȡdepth_library_name_
//		{
//			SpMat TransMatrix(TransMatrixDecay);
//			ModecNuclideLibrary.CalculateFlux(mode);
//			flux_vector_.push_back(ModecNuclideLibrary.flux_);
//			power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//			while (1)
//			{
//
//				TransMatrix = TransMatrixDecay + (TransMatrixCrossSection + TransMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
//				TransMatrix.SymbolLUElimination();
//				
//				
//					ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);
//					if (if_constant_online_feeding_ == true)
//					{
//						int size_F = constant_feeding_vector_.size();
//
//						vector<double > F_mol;
//						F_mol.resize(size_F);
//
//						int size_GL = gauss_legendre_abscissa_.size();
//						for (int GL_i = 0; GL_i < size_GL; ++GL_i)
//						{
//							vector<double > F_temp;
//							F_temp.resize(size_F);
//							for (int F_i = 0; F_i < size_F; ++F_i)
//							{
//								F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
//							}
//
//							double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);
//
//							F_temp = Solver.PfdCramSolver(TransMatrix, F_temp, time_gl);
//
//							for (int F_i = 0; F_i < size_F; ++F_i)
//							{
//								F_mol[F_i] += F_temp[F_i];
//							}
//						}
//
//						for (int F_i = 0; F_i < size_F; ++F_i)
//						{
//							ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
//						}
//
//					}
//
//					error_max = 0;
//					for (unsigned int i = 0; i < n_vector_.size(); ++i)
//					{
//						error_[i] = abs((ModecNuclideLibrary.nuclide_library_vector_[0][i] - n_vector_[1][i]) / n_vector_[1][i]);
//						if (error_[i] > error_max)
//						{
//							error_max = error_[i];
//						}
//					}
//
//					if (error_max < error_cutoff)
//					{
//						break;
//					}
//					else
//					{
//						n_vector_[1] = ModecNuclideLibrary.nuclide_library_vector_[0];
//					}
//
//				ConstructFissionYieldsSpMat(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
//
//				ModecNuclideLibrary.CalculateFlux(mode);
//			}
//			//ModecNuclideLibrary.CalculateFlux(mode);
//			flux_vector_.push_back(ModecNuclideLibrary.flux_);
//			power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//		}
//		else // �ڶ�ȡcouple�ļ�����ȼ�ľ���ʱ������Ҫ�����ѱ����ݶ��������Ҳ����ҪTransMatrixFissionYields�����Ѿ�������xs��
//		{
//			SpMat TransMatrix;
//			TransMatrix = TransMatrixDecay + TransMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);
//			TransMatrix.SymbolLUElimination();
//
//			ModecNuclideLibrary.CalculateFlux(mode);
//			flux_vector_.push_back(ModecNuclideLibrary.flux_);
//			power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//			while(1)
//			{
//				
//		
//					ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);
//					if (if_constant_online_feeding_ == true)
//					{
//						int size_F = constant_feeding_vector_.size();
//
//						vector<double > F_mol;
//						F_mol.resize(size_F);
//
//						int size_GL = gauss_legendre_abscissa_.size();
//						for (int GL_i = 0; GL_i < size_GL; ++GL_i)
//						{
//							vector<double > F_temp;
//							F_temp.resize(size_F);
//							for (int F_i = 0; F_i < size_F; ++F_i)
//							{
//								F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
//							}
//
//							double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);
//
//							F_temp = Solver.PfdCramSolver(TransMatrix, F_temp, time_gl);
//
//							for (int F_i = 0; F_i < size_F; ++F_i)
//							{
//								F_mol[F_i] += F_temp[F_i];
//							}
//						}
//
//						for (int F_i = 0; F_i < size_F; ++F_i)
//						{
//							ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
//						}
//
//					}
//
//					error_max = 0;
//					for (unsigned int i = 0; i < n_vector_.size(); ++i)
//					{
//						error_[i] = abs((ModecNuclideLibrary.nuclide_library_vector_[0][i] - n_vector_[1][i]) / n_vector_[1][i]);
//						if (error_[i] > error_max)
//						{
//							error_max = error_[i];
//						}
//					}
//
//					if (error_max < error_cutoff)
//					{
//						break;
//					}
//					else
//					{
//						n_vector_[1] = ModecNuclideLibrary.nuclide_library_vector_[0];
//					}
//
//					ModecNuclideLibrary.CalculateFlux(mode);
//			}
//			flux_vector_.push_back(ModecNuclideLibrary.flux_);
//			power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//		}
//	}
//
//	if (mode == 2) //���������
//	{
//		if (lib_tag_ == 1) //lib_tag_ = 1��ζ�Ŷ�ȡdepth_library_name_
//		{
//			SpMat TransMatrix(TransMatrixDecay);
//			ModecNuclideLibrary.CalculateFlux(mode);
//			flux_vector_.push_back(ModecNuclideLibrary.flux_);
//			power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//			while (1)
//			{
//
//				TransMatrix = TransMatrixDecay + (TransMatrixCrossSection + TransMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
//				TransMatrix.SymbolLUElimination();
//
//
//				ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);
//				if (if_constant_online_feeding_ == true)
//				{
//					int size_F = constant_feeding_vector_.size();
//
//					vector<double > F_mol;
//					F_mol.resize(size_F);
//
//					int size_GL = gauss_legendre_abscissa_.size();
//					for (int GL_i = 0; GL_i < size_GL; ++GL_i)
//					{
//						vector<double > F_temp;
//						F_temp.resize(size_F);
//						for (int F_i = 0; F_i < size_F; ++F_i)
//						{
//							F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
//						}
//
//						double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);
//
//						F_temp = Solver.PfdCramSolver(TransMatrix, F_temp, time_gl);
//
//						for (int F_i = 0; F_i < size_F; ++F_i)
//						{
//							F_mol[F_i] += F_temp[F_i];
//						}
//					}
//
//					for (int F_i = 0; F_i < size_F; ++F_i)
//					{
//						ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
//					}
//
//				}
//
//				error_max = 0;
//				for (unsigned int i = 0; i < n_vector_.size(); ++i)
//				{
//					error_[i] = abs((ModecNuclideLibrary.nuclide_library_vector_[0][i] - n_vector_[1][i]) / n_vector_[1][i]);
//					if (error_[i] > error_max)
//					{
//						error_max = error_[i];
//					}
//				}
//
//				if (error_max < error_cutoff)
//				{
//					break;
//				}
//				else
//				{
//					n_vector_[1] = ModecNuclideLibrary.nuclide_library_vector_[0];
//				}
//
//				ConstructFissionYieldsSpMat(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
//
//				ModecNuclideLibrary.CalculateFlux(1);
//			}
//			//ModecNuclideLibrary.CalculateFlux(mode);
//			flux_vector_.push_back(ModecNuclideLibrary.flux_);
//			power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//		}
//		else // �ڶ�ȡcouple�ļ�����ȼ�ľ���ʱ������Ҫ�����ѱ����ݶ��������Ҳ����ҪTransMatrixFissionYields�����Ѿ�������xs��
//		{
//			SpMat TransMatrix;
//			ModecNuclideLibrary.CalculateFlux(mode);
//
//			TransMatrix = TransMatrixDecay + TransMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);
//			TransMatrix.SymbolLUElimination();
//
//			
//			flux_vector_.push_back(ModecNuclideLibrary.flux_);
//			power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//			while (1)
//			{
//
//
//				ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);
//				if (if_constant_online_feeding_ == true)
//				{
//					int size_F = constant_feeding_vector_.size();
//
//					vector<double > F_mol;
//					F_mol.resize(size_F);
//
//					int size_GL = gauss_legendre_abscissa_.size();
//					for (int GL_i = 0; GL_i < size_GL; ++GL_i)
//					{
//						vector<double > F_temp;
//						F_temp.resize(size_F);
//						for (int F_i = 0; F_i < size_F; ++F_i)
//						{
//							F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
//						}
//
//						double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);
//
//						F_temp = Solver.PfdCramSolver(TransMatrix, F_temp, time_gl);
//
//						for (int F_i = 0; F_i < size_F; ++F_i)
//						{
//							F_mol[F_i] += F_temp[F_i];
//						}
//					}
//
//					for (int F_i = 0; F_i < size_F; ++F_i)
//					{
//						ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
//					}
//
//				}
//
//				error_max = 0;
//				for (unsigned int i = 0; i < n_vector_.size(); ++i)
//				{
//					error_[i] = abs((ModecNuclideLibrary.nuclide_library_vector_[0][i] - n_vector_[1][i]) / n_vector_[1][i]);
//					if (error_[i] > error_max)
//					{
//						error_max = error_[i];
//					}
//				}
//
//				if (error_max < error_cutoff)
//				{
//					break;
//				}
//				else
//				{
//					n_vector_[1] = ModecNuclideLibrary.nuclide_library_vector_[0];
//				}
//
//				ModecNuclideLibrary.CalculateFlux(1);
//			}
//			flux_vector_.push_back(ModecNuclideLibrary.flux_);
//			power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//		}
//	}
//}


void ModecClass::CalEquilibrium(int mode)
{
	double error_cutoff = 1.0e-10;
	
	vector<double> error_;
	error_.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size());
	
	double error_max;


	SpMat TransMatrix;
	ModecNuclideLibrary.CalculateFlux(mode);
	flux_vector_.push_back(ModecNuclideLibrary.flux_);
	power_vector_.push_back(ModecNuclideLibrary.specified_power_); // ��ʼʱ�̵�ͨ���͹���ͳ��

	if (lib_tag_ == 1) //lib_tag_ = 1��ζ�Ŷ�ȡdepth_library_name_
	{
		TransMatrix = TransMatrixDecay + (TransMatrixCrossSection + TransMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
	}
	else
	{
		TransMatrix = TransMatrixDecay + TransMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);
	}
	TransMatrix.SymbolLUElimination();


	int size_F = constant_feeding_vector_.size();

	int size_GL = gauss_legendre_abscissa_.size();

	double time = 1000.0 * 24 * 3600;
	double total_time = 0;

	int count = 0;

	for(;;)
	{
		vector<double > F_mol;
		F_mol.resize(size_F);

		for (int GL_i = 0; GL_i < size_GL; ++GL_i)
		{
			vector<double > F_temp;
			F_temp.resize(size_F);
			for (int F_i = 0; F_i < size_F; ++F_i)
			{
				F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
			}

			double time_gl = time / 2.0 * gauss_legendre_abscissa_[GL_i] + (total_time + time/2.0);

			F_temp = Solver.PfdCramSolver(TransMatrix, F_temp, time_gl);

			for (int F_i = 0; F_i < size_F; ++F_i)
			{
				F_mol[F_i] += F_temp[F_i];
			}
		}

		if (count > 10)
		{
			error_max = 0;
			for (int F_i = 0; F_i < size_F; ++F_i)
			{
				error_[F_i] = abs(F_mol[F_i] / n_vector_[1][F_i]); // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
				if (error_[F_i] > error_max)
				{
					error_max = error_[F_i];
				}
			}
			if (error_max < error_cutoff)
			{
				cout << "error_max = " << error_max << '\n';
				break;
			}
			else
			{
				for (int F_i = 0; F_i < size_F; ++F_i)
				{
					n_vector_[1][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
				}
			}
		}
		else
		{
			for (int F_i = 0; F_i < size_F; ++F_i)
			{
				n_vector_[1][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
			}
		}
		
		total_time += time;
		if (count < 20)
		{
			time += time;
		}
		count++;
	}
	
	
	ModecNuclideLibrary.CalculateFlux(mode);
	flux_vector_.push_back(ModecNuclideLibrary.flux_);
	power_vector_.push_back(ModecNuclideLibrary.specified_power_);// ƽ��̬��ͨ���͹���ͳ��
}

void ModecClass::Evolution(int mode, double time, int subtime)
{
	if (if_tracking_stockage == true)
	{
		TransMatrixStockage.SymbolLUElimination();
	}

	switch (mode)
	{
		case 0:	// ��˥�����
		{
			{
				if (solver_selection_ == 1)
				{
					if (if_constant_online_feeding_ == true && constant_feeding_calculation_methods_ == 2)
					{
						int size_matrix = TransMatrixDecay.spmat_dimen_;
						TransMatrixDecay.Resize(size_matrix + 1); // ��������
						int size_nucl = constant_feeding_nuclide_id_vector_.size();
						for (int i = 0; i < size_nucl; ++i)
						{
							int index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
							TransMatrixDecay.AddElement(index, size_matrix, constant_feeding_rate_[i]);
						}
						TransMatrixDecay.SymbolLUElimination();



						if (if_tracking_stockage == false)
						{
							vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
							F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 1.0);
							for (int i = 1; i <= subtime; ++i)
							{
								F_mol = Solver.PfdCramSolver(TransMatrixDecay, F_mol, time);
								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
								}
								F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size()] = 1.0;

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
								flux_vector_.push_back(0.0);
								power_vector_.push_back(0.0);
							}
						}
						else
						{
							vector<double > F_mol;
							F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1);
							for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
							{
								if (j < ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
								{
									F_mol[j] = ModecNuclideLibrary.nuclide_library_vector_[0][j];
								}
								if (j >= ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
								{
									F_mol[j + 1] = ModecNuclideLibrary.nuclide_library_vector_[0][j];
								}
							}
							F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2] = 1.0;

							for (int i = 1; i <= subtime; ++i)
							{
								F_mol = Solver.PfdCramSolver(TransMatrixDecay, TransMatrixReprocess, TransMatrixStockage, F_mol, time);
								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									if (j < ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
									{
										ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
									}
									if (j >= ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
									{
										ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j + 1];
									}
								}
								F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2] = 1.0;

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
								flux_vector_.push_back(0.0);
								power_vector_.push_back(0.0);
							}
						}
					}
					else
					{
						TransMatrixDecay.SymbolLUElimination();
						for (int i = 1; i <= subtime; ++i)
						{
							if (if_tracking_stockage == true)
							{
								ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrixDecay, TransMatrixReprocess, TransMatrixStockage, ModecNuclideLibrary.nuclide_library_vector_[0], time);
								if (if_constant_online_feeding_ == true)
								{
									int size_F = constant_feeding_vector_.size();

									vector<double > F_mol;
									F_mol.resize(size_F);

									int size_GL = gauss_legendre_abscissa_.size();
									for (int GL_i = 0; GL_i < size_GL; ++GL_i)
									{
										vector<double > F_temp;
										F_temp.resize(size_F);
										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
										}

										double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);

										F_temp = Solver.PfdCramSolver(TransMatrixDecay, TransMatrixReprocess, TransMatrixStockage, F_temp, time_gl);

										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											F_mol[F_i] += F_temp[F_i];
										}
									}

									for (int F_i = 0; F_i < size_F; ++F_i)
									{
										ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
									}

								}
							}
							else
							{
								ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrixDecay, ModecNuclideLibrary.nuclide_library_vector_[0], time);
								if (if_constant_online_feeding_ == true)
								{

									int size_F = constant_feeding_vector_.size();

									vector<double > F_mol;
									F_mol.resize(size_F);

									int size_GL = gauss_legendre_abscissa_.size();
									for (int GL_i = 0; GL_i < size_GL; ++GL_i)
									{
										vector<double > F_temp;
										F_temp.resize(size_F);
										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
										}

										double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);

										F_temp = Solver.PfdCramSolver(TransMatrixDecay, F_temp, time_gl);

										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											F_mol[F_i] += F_temp[F_i];
										}
									}

									for (int F_i = 0; F_i < size_F; ++F_i)
									{
										ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
									}

								}
							}

							n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
							flux_vector_.push_back(0.0);
							power_vector_.push_back(0.0);
						}
					}
				}
				else if (solver_selection_ == 0)
				{
					if (if_constant_online_feeding_ == false) {
						Solver.TtaInitialize(ModecNuclideLibrary.nuclide_library_vector_[0]);
						double re_time = 0;
						for (int i = 1; i <= subtime; ++i)
						{
							re_time += time;
							ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.TtaSolver(TtaMatrixDecay, re_time);
							n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
							flux_vector_.push_back(0.0);
							power_vector_.push_back(0.0);
						}
					}
					else
					{
						// TTA��������������ʵ�� 
						// ����һ�����أ�����ȼ���ڽӾ�������һά
						int size_matrix(TtaMatrixDecay.spmat_dimen_);
						TtaMatrixDecay.Resize(size_matrix + 1);

						int size_nucl(constant_feeding_nuclide_id_vector_.size());
						double tot_feeding_rate(0.0);
						for (int i = 0; i < size_nucl; ++i)
						{
							int index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
							TtaMatrixDecay.AddElementCCS(index, size_matrix, constant_feeding_rate_[i]);
							tot_feeding_rate += constant_feeding_rate_[i];
						}
						TtaMatrixDecay.AddElementCCS(size_matrix, size_matrix, -tot_feeding_rate);
						Solver.tot_feeding_rate_ = tot_feeding_rate;

						vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
						F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 0.0);

						Solver.TtaInitialize(F_mol);
						double re_time = 0;
						for (int i = 1; i <= subtime; ++i)
						{
							re_time += time;
							F_mol = Solver.TtaSolver(TtaMatrixDecay, re_time);
							for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
							{
								ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
							}

							n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
							flux_vector_.push_back(0.0);
							power_vector_.push_back(0.0);
						}
					}
				}
			}
			break;
		}

		case 1:	// ��ͨ�����
		{
			{
				if (lib_tag_ == 1) //lib_tag_ = 1��ζ�Ŷ�ȡDepthLib
				{
					if(solver_selection_ == 1)
					{
						if (if_constant_online_feeding_ == true && constant_feeding_calculation_methods_ == 2)
						{
							int size_matrix = TransMatrixDecay.spmat_dimen_;
							SpMat TransMatrix(size_matrix + 1);
							TransMatrixDecay.Resize(size_matrix + 1);
							TransMatrixCrossSection.Resize(size_matrix + 1);
							TransMatrixFissionYields.Resize(size_matrix + 1);

							int size_nucl = constant_feeding_nuclide_id_vector_.size();
							for (int i = 0; i < size_nucl; ++i)
							{
								int index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
								TransMatrixDecay.AddElement(index, size_matrix, constant_feeding_rate_[i]);
							}

							if (if_tracking_stockage == false)
							{
								vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
								F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 1.0);

								ModecNuclideLibrary.CalculateFlux(mode);

								for (int i = 1; i <= subtime; ++i)
								{

									TransMatrix = TransMatrixDecay + (TransMatrixCrossSection + TransMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
									TransMatrix.SymbolLUElimination();

									F_mol = Solver.PfdCramSolver(TransMatrix, F_mol, time);
									for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
									{
										ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
									}
									F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size()] = 1.0;

									n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

									ModecNuclideLibrary.CalculateFlux(mode);
									flux_vector_.push_back(ModecNuclideLibrary.flux_);
									power_vector_.push_back(ModecNuclideLibrary.specified_power_);

									ConstructFissionYieldsSpMat(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
								}
							}
							else
							{
								vector<double > F_mol;
								F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1);
								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									if (j < ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
									{
										F_mol[j] = ModecNuclideLibrary.nuclide_library_vector_[0][j];
									}
									if (j >= ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
									{
										F_mol[j + 1] = ModecNuclideLibrary.nuclide_library_vector_[0][j];
									}
								}
								F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2] = 1.0;
								ModecNuclideLibrary.CalculateFlux(mode);

								for (int i = 1; i <= subtime; ++i)
								{
									TransMatrix = TransMatrixDecay + (TransMatrixCrossSection + TransMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
									TransMatrix.SymbolLUElimination();

									F_mol = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, F_mol, time);
									for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
									{
										if (j < ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
										}
										if (j >= ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j + 1];
										}
									}
									F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2] = 1.0;

									n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

									ModecNuclideLibrary.CalculateFlux(mode);
									flux_vector_.push_back(ModecNuclideLibrary.flux_);
									power_vector_.push_back(ModecNuclideLibrary.specified_power_);
									ConstructFissionYieldsSpMat(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
								}
							}

						}
						else
						{
							SpMat TransMatrix;// (TransMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);

							for (int i = 1; i <= subtime; ++i)
							{

								TransMatrix = TransMatrixDecay + (TransMatrixCrossSection + TransMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
								//SpMat TransMatrix(TransMatrixDecay + (TransMatrixCrossSection + TransMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24));
								TransMatrix.SymbolLUElimination();
								if (if_tracking_stockage == true)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, ModecNuclideLibrary.nuclide_library_vector_[0], time);
									if (if_constant_online_feeding_ == true)
									{
										int size_F = constant_feeding_vector_.size();

										vector<double > F_mol;
										F_mol.resize(size_F);

										int size_GL = gauss_legendre_abscissa_.size();
										for (int GL_i = 0; GL_i < size_GL; ++GL_i)
										{
											vector<double > F_temp;
											F_temp.resize(size_F);
											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
											}

											double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);

											F_temp = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, F_temp, time_gl);

											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_mol[F_i] += F_temp[F_i];
											}
										}

										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
										}

									}
								}
								else
								{
									ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);
									if (if_constant_online_feeding_ == true)
									{
										int size_F = constant_feeding_vector_.size();

										vector<double > F_mol;
										F_mol.resize(size_F);

										int size_GL = gauss_legendre_abscissa_.size();
										for (int GL_i = 0; GL_i < size_GL; ++GL_i)
										{
											vector<double > F_temp;
											F_temp.resize(size_F);
											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
											}

											double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);

											F_temp = Solver.PfdCramSolver(TransMatrix, F_temp, time_gl);

											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_mol[F_i] += F_temp[F_i];
											}
										}

										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
										}

									}
								}

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);

								ConstructFissionYieldsSpMat(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
							}

						}
					}
					else if (solver_selection_ == 0)
					{
						if (if_constant_online_feeding_ == false) {
							//Solver.TtaInitialize(ModecNuclideLibrary.nuclide_library_vector_[0]);
							Solver.TtaInitialize(ModecNuclideLibrary.nuclide_library_vector_[0].size());
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);

							//double re_time = 0;

							for (int i = 1; i <= subtime; ++i)
							{
								//re_time += time;
								TransMatrix = TtaMatrixDecay + (TtaMatrixCrossSection + TtaMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
								Solver.TtaSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);

								ConstructFissionYieldsSpMatForTta(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
							}
						}
						else
						{
							// TTA��������������ʵ�� 
							// ����һ�����أ�����ȼ���ڽӾ�������һά
							int size_matrix(TtaMatrixDecay.spmat_dimen_);
							TtaMatrixDecay.Resize(size_matrix + 1);
							TtaMatrixCrossSection.Resize(size_matrix + 1);
							TtaMatrixFissionYields.Resize(size_matrix + 1);

							int size_nucl(constant_feeding_nuclide_id_vector_.size());
							double tot_feeding_rate(0.0);
							for (int i = 0; i < size_nucl; ++i)
							{
								int index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
								//TtaMatrixDecay.AddElementCCS(index, size_matrix, constant_feeding_rate_[i]);
								Solver.feed_nuclide_id_.push_back(index);

								Solver.feed_rate_.push_back(constant_feeding_rate_[i]);
								tot_feeding_rate += abs(constant_feeding_rate_[i]);
							}


							//TtaMatrixDecay.AddElementCCS(size_matrix, size_matrix, -tot_feeding_rate);
							Solver.tot_feeding_rate_ = tot_feeding_rate;
							//Solver.feed_nuclide_id_ = constant_feeding_nuclide_id_vector_;

							vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
							F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 0.0);

							Solver.TtaInitialize(F_mol.size());
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);

							for (int i = 1; i <= subtime; ++i)
							{
								TransMatrix = TtaMatrixDecay + (TtaMatrixCrossSection + TtaMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
								Solver.TtaSolverForFeeding(TransMatrix, F_mol, time);

								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
								}

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);

								ConstructFissionYieldsSpMatForTta(); // modify fission-yields at the end of each substep 
							}
						}
					}

				}
				else // �ڶ�ȡcouple�ļ�����ȼ�ľ���ʱ������Ҫ�����ѱ����ݶ��������Ҳ����ҪTransMatrixFissionYields�����Ѿ�������xs��
				{
					if (solver_selection_ == 1)
					{
						if (if_constant_online_feeding_ == true && constant_feeding_calculation_methods_ == 2)
						{
							int size_matrix = TransMatrixDecay.spmat_dimen_;
							SpMat TransMatrix(size_matrix + 1);
							TransMatrixDecay.Resize(size_matrix + 1);
							TransMatrixCrossSection.Resize(size_matrix + 1);

							TransMatrix = TransMatrixDecay + TransMatrixCrossSection * (ModecNuclideLibrary.flux_ * 1.0e-24);

							int size_nucl = constant_feeding_nuclide_id_vector_.size();
							for (int i = 0; i < size_nucl; ++i)
							{
								int index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
								TransMatrix.AddElement(index, size_matrix, constant_feeding_rate_[i]);
							}

							TransMatrix.SymbolLUElimination();

							if (if_tracking_stockage == false)
							{
								vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
								F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 1.0);
								ModecNuclideLibrary.CalculateFlux(mode);

								for (int i = 1; i <= subtime; ++i)
								{

									F_mol = Solver.PfdCramSolver(TransMatrix, F_mol, time);
									for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
									{
										ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
									}
									F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size()] = 1.0;

									n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

									ModecNuclideLibrary.CalculateFlux(mode);
									flux_vector_.push_back(ModecNuclideLibrary.flux_);
									power_vector_.push_back(ModecNuclideLibrary.specified_power_);
								}
							}
							else
							{
								vector<double > F_mol;
								F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1);
								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									if (j < ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
									{
										F_mol[j] = ModecNuclideLibrary.nuclide_library_vector_[0][j];
									}
									if (j >= ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
									{
										F_mol[j + 1] = ModecNuclideLibrary.nuclide_library_vector_[0][j];
									}
								}
								F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2] = 1.0;
								ModecNuclideLibrary.CalculateFlux(mode);

								for (int i = 1; i <= subtime; ++i)
								{

									F_mol = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, F_mol, time);
									for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
									{
										if (j < ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
										}
										if (j >= ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j + 1];
										}
									}
									F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2] = 1.0;

									n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

									ModecNuclideLibrary.CalculateFlux(mode);
									flux_vector_.push_back(ModecNuclideLibrary.flux_);
									power_vector_.push_back(ModecNuclideLibrary.specified_power_);
								}
							}

						}
						else
						{
							SpMat TransMatrix;
							TransMatrix = TransMatrixDecay + TransMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);
							//TransitionMatrixOutput(TransMatrixDecay*time); // �������ĸ���Ԫ�ص��ļ���

							TransMatrix.SymbolLUElimination();
							ModecNuclideLibrary.CalculateFlux(mode);

							for (int i = 1; i <= subtime; ++i)
							{

								if (if_tracking_stockage == true)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, ModecNuclideLibrary.nuclide_library_vector_[0], time);
									if (if_constant_online_feeding_ == true)
									{
										int size_F = constant_feeding_vector_.size();

										vector<double > F_mol;
										F_mol.resize(size_F);

										int size_GL = gauss_legendre_abscissa_.size();
										for (int GL_i = 0; GL_i < size_GL; ++GL_i)
										{
											vector<double > F_temp;
											F_temp.resize(size_F);
											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
											}

											double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);

											F_temp = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, F_temp, time_gl);

											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_mol[F_i] += F_temp[F_i];
											}
										}

										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
										}

									}
								}
								else
								{
									ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);
									if (if_constant_online_feeding_ == true)
									{
										int size_F = constant_feeding_vector_.size();

										vector<double > F_mol;
										F_mol.resize(size_F);

										int size_GL = gauss_legendre_abscissa_.size();
										for (int GL_i = 0; GL_i < size_GL; ++GL_i)
										{
											vector<double > F_temp;
											F_temp.resize(size_F);
											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
											}

											double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);

											F_temp = Solver.PfdCramSolver(TransMatrix, F_temp, time_gl);

											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_mol[F_i] += F_temp[F_i];
											}
										}

										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
										}

									}
								}
								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
							}

						}
					}
					else if (solver_selection_ == 0)
					{
						if (if_constant_online_feeding_ == false) {
							Solver.TtaInitialize(ModecNuclideLibrary.nuclide_library_vector_[0].size());
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);
							TransMatrix = TtaMatrixDecay + TtaMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);

							for (int i = 1; i <= subtime; ++i)
							{
								Solver.TtaSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
							}
						}
						else
						{
							int size_matrix(TtaMatrixDecay.spmat_dimen_);
							TtaMatrixDecay.Resize(size_matrix + 1);
							TtaMatrixCrossSection.Resize(size_matrix + 1);
							
							int size_nucl(constant_feeding_nuclide_id_vector_.size());
							double tot_feeding_rate(0.0);
							for (int i = 0; i < size_nucl; ++i)
							{
								int index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
								//TtaMatrixDecay.AddElementCCS(index, size_matrix, constant_feeding_rate_[i]);
								Solver.feed_nuclide_id_.push_back(index);
								Solver.feed_rate_.push_back(constant_feeding_rate_[i]);	
								tot_feeding_rate += abs(constant_feeding_rate_[i]);
							}
							//TtaMatrixDecay.AddElementCCS(size_matrix, size_matrix, -tot_feeding_rate);
							Solver.tot_feeding_rate_ = tot_feeding_rate;

							vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
							F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 0.0);

							Solver.TtaInitialize(F_mol.size());
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);
							TransMatrix = TtaMatrixDecay + TtaMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);

							for (int i = 1; i <= subtime; ++i)
							{
								Solver.TtaSolver(TransMatrix, F_mol, time);

								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
								}

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
							}
						}
					}
				}
			}
			break;
		}
	
		case 2:	//���������
		{		
			{
				if (lib_tag_ == 1)
				{
					if (solver_selection_ == 1)
					{
						if (if_constant_online_feeding_ == true && constant_feeding_calculation_methods_ == 2)
						{
							int size_matrix = TransMatrixDecay.spmat_dimen_;
							SpMat TransMatrix(size_matrix + 1);
							TransMatrixDecay.Resize(size_matrix + 1);
							TransMatrixCrossSection.Resize(size_matrix + 1);
							TransMatrixFissionYields.Resize(size_matrix + 1);

							int size_nucl = constant_feeding_nuclide_id_vector_.size();
							for (int i = 0; i < size_nucl; ++i)
							{
								int index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
								TransMatrixDecay.AddElement(index, size_matrix, constant_feeding_rate_[i]);
							}

							if (if_tracking_stockage == false)
							{
								vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
								F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 1.0);
								ModecNuclideLibrary.CalculateFlux(mode);

								for (int i = 1; i <= subtime; ++i)
								{

									TransMatrix = TransMatrixDecay + (TransMatrixCrossSection + TransMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
									TransMatrix.SymbolLUElimination();

									F_mol = Solver.PfdCramSolver(TransMatrix, F_mol, time);
									for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
									{
										ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];

									}
									F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size()] = 1.0;

									n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

									ModecNuclideLibrary.CalculateFlux(mode);
									flux_vector_.push_back(ModecNuclideLibrary.flux_);
									power_vector_.push_back(ModecNuclideLibrary.specified_power_);
									ConstructFissionYieldsSpMat(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
								}
							}
							else
							{
								vector<double > F_mol;
								F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1);
								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									if (j < ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
									{
										F_mol[j] = ModecNuclideLibrary.nuclide_library_vector_[0][j];
									}
									if (j >= ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
									{
										F_mol[j + 1] = ModecNuclideLibrary.nuclide_library_vector_[0][j];
									}
								}
								F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2] = 1.0;

								ModecNuclideLibrary.CalculateFlux(mode);

								for (int i = 1; i <= subtime; ++i)
								{

									TransMatrix = TransMatrixDecay + (TransMatrixCrossSection + TransMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
									TransMatrix.SymbolLUElimination();

									F_mol = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, F_mol, time);
									for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
									{
										if (j < ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
										}
										if (j >= ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j + 1];
										}
									}
									F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2] = 1.0;

									n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

									ModecNuclideLibrary.CalculateFlux(mode);
									flux_vector_.push_back(ModecNuclideLibrary.flux_);
									power_vector_.push_back(ModecNuclideLibrary.specified_power_);
									ConstructFissionYieldsSpMat(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
								}
							}

						}
						else
						{
							SpMat TransMatrix;
							ModecNuclideLibrary.CalculateFlux(mode);

							for (int i = 1; i <= subtime; ++i)
							{

								TransMatrix = (TransMatrixCrossSection + TransMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24) + TransMatrixDecay;
								TransMatrix.SymbolLUElimination();
								if (if_tracking_stockage == true)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, ModecNuclideLibrary.nuclide_library_vector_[0], time);
									if (if_constant_online_feeding_ == true)
									{
										int size_F = constant_feeding_vector_.size();

										vector<double > F_mol;
										F_mol.resize(size_F);

										int size_GL = gauss_legendre_abscissa_.size();
										for (int GL_i = 0; GL_i < size_GL; ++GL_i)
										{
											vector<double > F_temp;
											F_temp.resize(size_F);
											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
											}

											double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);

											F_temp = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, F_temp, time_gl);

											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_mol[F_i] += F_temp[F_i];
											}
										}

										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
										}

									}
								}
								else
								{
									ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);
									if (if_constant_online_feeding_ == true)
									{
										int size_F = constant_feeding_vector_.size();

										vector<double > F_mol;
										F_mol.resize(size_F);

										int size_GL = gauss_legendre_abscissa_.size();
										for (int GL_i = 0; GL_i < size_GL; ++GL_i)
										{
											vector<double > F_temp;
											F_temp.resize(size_F);
											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
											}

											double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);

											F_temp = Solver.PfdCramSolver(TransMatrix, F_temp, time_gl);

											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_mol[F_i] += F_temp[F_i];
											}
										}

										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
										}

									}
								}
								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
								ConstructFissionYieldsSpMat(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
							}

						}
					}
					else if (solver_selection_ == 0)
					{
						if (if_constant_online_feeding_ == false) {
							Solver.TtaInitialize(ModecNuclideLibrary.nuclide_library_vector_[0].size());
							SparseMatrixMCS TransMatrix;
							ModecNuclideLibrary.CalculateFlux(mode);

							for (int i = 1; i <= subtime; ++i)
							{
								TransMatrix = TtaMatrixDecay + (TtaMatrixCrossSection + TtaMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);

								Solver.TtaSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
								ConstructFissionYieldsSpMatForTta(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
							}
						}
						else
						{
							int size_matrix(TtaMatrixDecay.spmat_dimen_);
							TtaMatrixDecay.Resize(size_matrix + 1);
							TtaMatrixCrossSection.Resize(size_matrix + 1);
							TtaMatrixFissionYields.Resize(size_matrix + 1);

							int size_nucl(constant_feeding_nuclide_id_vector_.size());
							double tot_feeding_rate(0.0);
							for (int i = 0; i < size_nucl; ++i)
							{
								int index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
								//TtaMatrixDecay.AddElementCCS(index, size_matrix, constant_feeding_rate_[i]);
								Solver.feed_nuclide_id_.push_back(index);
								Solver.feed_rate_.push_back(constant_feeding_rate_[i]);

								tot_feeding_rate += abs(constant_feeding_rate_[i]);
							}
							//TtaMatrixDecay.AddElementCCS(size_matrix, size_matrix, -tot_feeding_rate);
							Solver.tot_feeding_rate_ = tot_feeding_rate;

							vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
							F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 0.0);

							Solver.TtaInitialize(F_mol.size());
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);

							for (int i = 1; i <= subtime; ++i)
							{
								TransMatrix = (TtaMatrixCrossSection + TtaMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24) + TtaMatrixDecay;

								Solver.TtaSolverForFeeding(TransMatrix, F_mol, time);

								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
								}

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
								ConstructFissionYieldsSpMatForTta(); // ÿ��ȼ�Ĳ������ѱ����ݶ�
							}
						}
					}

				}
				else// �ڶ�ȡcouple�ļ�����ȼ�ľ���ʱ������Ҫ�����ѱ����ݶ��������Ҳ����ҪTransMatrixFissionYields�����Ѿ�������xs��
				{
					if (solver_selection_ == 1)
					{
						if (if_constant_online_feeding_ == true && constant_feeding_calculation_methods_ == 2)
						{
							int size_matrix = TransMatrixDecay.spmat_dimen_;
							SpMat TransMatrix(size_matrix + 1);
							vector<int> _IRC;
							vector<int> _ICFR;
							vector<int> _LUP;

							TransMatrixDecay.Resize(size_matrix + 1);
							TransMatrixCrossSection.Resize(size_matrix + 1);


							int size_nucl = constant_feeding_nuclide_id_vector_.size();
							for (int i = 0; i < size_nucl; ++i)
							{
								int index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
								TransMatrixDecay.AddElement(index, size_matrix, constant_feeding_rate_[i]);
							}

							if (if_tracking_stockage == false)
							{
								vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
								F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 1.0);
								ModecNuclideLibrary.CalculateFlux(mode);

								for (int i = 1; i <= subtime; ++i)
								{

									TransMatrix = TransMatrixCrossSection * (ModecNuclideLibrary.flux_ * 1.0e-24) + TransMatrixDecay;

									if (i == 1)
									{
										TransMatrix.SymbolLUElimination();
										_IRC = TransMatrix.IRC;
										_ICFR = TransMatrix.ICFR;
										_LUP = TransMatrix.LUP;
									}
									else
									{
										TransMatrix.IRC = _IRC;
										TransMatrix.ICFR = _ICFR;
										TransMatrix.LUP = _LUP;
									}
									F_mol = Solver.PfdCramSolver(TransMatrix, F_mol, time);
									for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
									{
										ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
									}
									F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size()] = 1.0;

									n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

									ModecNuclideLibrary.CalculateFlux(mode);
									flux_vector_.push_back(ModecNuclideLibrary.flux_);
									power_vector_.push_back(ModecNuclideLibrary.specified_power_);
								}
							}
							else
							{
								vector<double > F_mol;
								F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1);
								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									if (j < ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
									{
										F_mol[j] = ModecNuclideLibrary.nuclide_library_vector_[0][j];
									}
									if (j >= ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
									{
										F_mol[j + 1] = ModecNuclideLibrary.nuclide_library_vector_[0][j];
									}
								}
								F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2] = 1.0;

								ModecNuclideLibrary.CalculateFlux(mode);

								for (int i = 1; i <= subtime; ++i)
								{

									TransMatrix = TransMatrixCrossSection * (ModecNuclideLibrary.flux_ * 1.0e-24) + TransMatrixDecay;

									if (i == 1)
									{
										TransMatrix.SymbolLUElimination();
										_IRC = TransMatrix.IRC;
										_ICFR = TransMatrix.ICFR;
										_LUP = TransMatrix.LUP;
									}
									else
									{
										TransMatrix.IRC = _IRC;
										TransMatrix.ICFR = _ICFR;
										TransMatrix.LUP = _LUP;
									}
									F_mol = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, F_mol, time);
									for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
									{
										if (j < ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
										}
										if (j >= ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j + 1];
										}
									}
									F_mol[ModecNuclideLibrary.nuclide_library_vector_[0].size() / 2] = 1.0;

									n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

									ModecNuclideLibrary.CalculateFlux(mode);
									flux_vector_.push_back(ModecNuclideLibrary.flux_);
									power_vector_.push_back(ModecNuclideLibrary.specified_power_);
								}
							}

						}
						else
						{
							SpMat TransMatrix;
							vector<int> _IRC;
							vector<int> _ICFR;
							vector<int> _LUP;

							ModecNuclideLibrary.CalculateFlux(mode);
							for (int i = 1; i <= subtime; ++i)
							{

								TransMatrix = TransMatrixCrossSection * (ModecNuclideLibrary.flux_ * 1.0e-24) + TransMatrixDecay;

								if (i == 1)
								{
									TransMatrix.SymbolLUElimination();
									_IRC = TransMatrix.IRC;
									_ICFR = TransMatrix.ICFR;
									_LUP = TransMatrix.LUP;
								}
								else
								{
									TransMatrix.IRC = _IRC;
									TransMatrix.ICFR = _ICFR;
									TransMatrix.LUP = _LUP;
								}

								if (if_tracking_stockage == true)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, ModecNuclideLibrary.nuclide_library_vector_[0], time);
									if (if_constant_online_feeding_ == true)
									{
										int size_F = constant_feeding_vector_.size();

										vector<double > F_mol;
										F_mol.resize(size_F);

										int size_GL = gauss_legendre_abscissa_.size();
										for (int GL_i = 0; GL_i < size_GL; ++GL_i)
										{
											vector<double > F_temp;
											F_temp.resize(size_F);
											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
											}

											double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);

											F_temp = Solver.PfdCramSolver(TransMatrix, TransMatrixReprocess, TransMatrixStockage, F_temp, time_gl);

											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_mol[F_i] += F_temp[F_i];
											}
										}

										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
										}

									}
								}
								else
								{
									ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.PfdCramSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);
									if (if_constant_online_feeding_ == true)
									{
										int size_F = constant_feeding_vector_.size();

										vector<double > F_mol;
										F_mol.resize(size_F);

										int size_GL = gauss_legendre_abscissa_.size();
										for (int GL_i = 0; GL_i < size_GL; ++GL_i)
										{
											vector<double > F_temp;
											F_temp.resize(size_F);
											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_temp[F_i] = constant_feeding_vector_[F_i] * gauss_legendre_weight_[GL_i] * time / 2.0;
											}

											double time_gl = time / 2.0*(1 - gauss_legendre_abscissa_[GL_i]);

											F_temp = Solver.PfdCramSolver(TransMatrix, F_temp, time_gl);

											for (int F_i = 0; F_i < size_F; ++F_i)
											{
												F_mol[F_i] += F_temp[F_i];
											}
										}

										for (int F_i = 0; F_i < size_F; ++F_i)
										{
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // �������ʳ����Ĺ��׼����ܵĺ���Ũ����ȥ
										}

									}
								}
								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);

								//ConstructFissionYieldsSpMat(TransMatrixFissionYields, ModecNuclideLibrary); // ÿ��ȼ�Ĳ������ѱ����ݶ�
							}

						}
					}
					else if (solver_selection_ == 0)
					{
						if (if_constant_online_feeding_ == false) {
							Solver.TtaInitialize(ModecNuclideLibrary.nuclide_library_vector_[0].size());
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);

							for (int i = 1; i <= subtime; ++i)
							{
								TransMatrix = TtaMatrixDecay + TtaMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);
								Solver.TtaSolver(TransMatrix, ModecNuclideLibrary.nuclide_library_vector_[0], time);

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
							}
						}
						else
						{
							int size_matrix(TtaMatrixDecay.spmat_dimen_);
							TtaMatrixDecay.Resize(size_matrix + 1);
							TtaMatrixCrossSection.Resize(size_matrix + 1);

							int size_nucl(constant_feeding_nuclide_id_vector_.size());
							double tot_feeding_rate(0.0);
							for (int i = 0; i < size_nucl; ++i)
							{
								int index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
								//TtaMatrixDecay.AddElementCCS(index, size_matrix, constant_feeding_rate_[i]);
								Solver.feed_nuclide_id_.push_back(index);
								Solver.feed_rate_.push_back(constant_feeding_rate_[i]);
								tot_feeding_rate += abs(constant_feeding_rate_[i]);
							}
							Solver.tot_feeding_rate_ = tot_feeding_rate;

							vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
							F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 0.0);

							Solver.TtaInitialize(F_mol.size());
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);

							for (int i = 1; i <= subtime; ++i)
							{
								TransMatrix = TtaMatrixDecay + TtaMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);
								Solver.TtaSolverForFeeding(TransMatrix, F_mol, time);

								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
								}

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
							}
						}
					}
				}
			}
			break;
		}
		case 3: // ��˥��+����
		{
			int num_depletion_zone( residue_time_.size() ); // ȼ��������
			int spmat_dimen( TransMatrixDecay.spmat_dimen_ );
			SpMat MatrixFlow( num_depletion_zone * spmat_dimen );
			vector< complex<double> > temp_mol;
			temp_mol.resize( num_depletion_zone * spmat_dimen);
			for ( int i = 0; i < spmat_dimen; ++i )
			{
				temp_mol[2 * i] = ModecNuclideLibrary.nuclide_library_vector_[0][i];
				temp_mol[2 * i + 1] = ModecNuclideLibrary.nuclide_library_vector_[0][i];  
			}
			
			n_vector_.resize(0);
			n_vector_.push_back(temp_mol);

			for (int _row = 0; _row < spmat_dimen; ++ _row)
			{
				for( int _col = 0; _col < spmat_dimen; ++ _col)
				{
					double element = TransMatrixDecay.Element(_row,_col);
					if (element != 0.0)
					{
						MatrixFlow.AddElement( 2 * _row, 2 * _col, element );
						MatrixFlow.AddElement( 2 * _row + 1, 2 * _col + 1, element );
					}
				}
			}
			for ( int _row = 0; _row < spmat_dimen; ++_row)
			{
				MatrixFlow.AddElement( 2 * _row, 2 * _row, -1 / residue_time[0]);
				MatrixFlow.AddElement( 2 * _row + 1, 2 * _row + 1, -1 / residue_time_[1]);
				MatrixFlow.AddElement( 2 * _row, 2 * _row + 1, 1 / residue_time[0]);
				MatrixFlow.AddElement( 2 * _row + 1, 2 * _row, 1 / residue_time[1]);
			}
			MatrixFlow.SymbolLUElimination();
			
			for ( int i = 1; i < subtime; ++ i)
			{
				temp_mol = Solver.PfdCramSolver(MatrixFlow, temp_mol, time);
				n_vector_.push_back(temp_mol);
				flux_vector_.push_back(0);
				power_vector_.push_back(0);
			}
			break;
		}
		case 4: // ��ͨ��+����
		{
			if( lib_tag_ == 1) // ��ȡDEPTH���ݿ�
			{
				int num_depletion_zone( residue_time_size() ); // ȼ��������
				int spmat_dimen( TransMatrixDecay.spmat_dimen_ );
				SpMat MatrixFlow( num_depletion_zone * spmat_dimen );
				vector<int> _IRC;
				vector<int> _ICFR;
				vector<int> _LUP;
				
				vector< complex<double> > temp_mol;
				temp_mol.resize(num_depletion_zone * spmat_dimen);
				for( int i = 0; i < spmat_dimen; ++i)
				{
					temp_mol[2 * i] = ModecNuclideLibrary.nuclide_library_vector_[0][i];
					temp_mol[2 * i + 1] = ModecNuclideLibrary.nuclide_library_vector_[0][i];
				}
			
				n_vector_.resize(0);
				n_vector_.push_back(temp_mol);	
			
				for (int _row = 0; _row < spmat_dimen; ++ _row)
				{	
					MatrixFlow.AddElement( 2 * _row, 2 * _row, -1 / residue_time[0]);
					MatrixFlow.AddElement( 2 * _row + 1, 2 * _row + 1, -1 / residue_time_[1]);
					MatrixFlow.AddElement( 2 * _row, 2 * _row + 1, 1 / residue_time[0]);
					MatrixFlow.AddElement( 2 * _row + 1, 2 * _row, 1 / residue_time[1]);
				}

				SpMat TransMatrix;
				
				ModecNuclideLibrary.CalculateFlux( mode - 3 ); // ��������ȼ�ĵı�־
				for (int i = 1; i <= subtime; ++i)
				{
					TransMatrix = TransMatrixDecay + (TransMatrixCrossSection + TransMatrixFissionYields) * (ModecNuclideLibrary.flux_ * 1.0e-24);
					for ( int _row = 0; _row < spmat_dimen; ++ _row )
					{
						for( int _col = 0; _col < spmat_dimen; ++ _col)
						{
							double element1 = TransMatrix.Element(_row, _col);
							double element2 = TransDecayMatrix.Element(_row, _col);
							if ( element1 != 0.0)
							{
								MatrixFlow.AddElement(2 * _row, 2 * _col, element1);
							}
							if ( element2 != 0.0)
							{
								MatrixFlow.AddElement(2 * _row + 1, 2 * _col + 1, element2);
							}
						}
					}
					if ( i == 1)
					{
						MatrixFlow.SymbolLUElimination();
						_IRC = MatrixFlow.IRC;
						_ICFR = MatrixFlow.ICFR;
						_LUP = MatrixFlow.LUP;
					}
					else
					{
						MatrixFlow.IRC = _IRC;
						MatrixFlow.ICFR = _ICFR;
						MatrixFlow.LUP = _LUP;
					}
					
					temp_mol = Solver.PfdCramSolver(MatrixFlow, temp_mol, time);
					
					n_vector_.push_back(temp_mol);
					
					for (int i = 0; i < spmat_dimen; ++i)
					{
						ModecNuclideLibrary.nuclide_library_vector_[0][i] = temp_mol[2 * i];
					}
					ModecNuclideLibrary.CalculateFlux( mode - 3);
					flux_vector_.push_back(ModecNuclideLibrary.flux_);
					power_vector_.push_back(ModecNuclideLibrary.specified_power_);
				}
		}
		case 5: // ������+����
		{
			
		}
	}
}
