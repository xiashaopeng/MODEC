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
		Evolution(mode, time, subtime);
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
//	if (mode == 0) // 纯衰变情况
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
//					ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
//	if (mode == 1) // 定通量情况
//	{
//		if (lib_tag_ == 1) //lib_tag_ = 1意味着读取depth_library_name_
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
//							ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
//				ConstructFissionYieldsSpMat(); // 每个燃耗步调整裂变产物份额
//
//				ModecNuclideLibrary.CalculateFlux(mode);
//			}
//			//ModecNuclideLibrary.CalculateFlux(mode);
//			flux_vector_.push_back(ModecNuclideLibrary.flux_);
//			power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//		}
//		else // 在读取couple文件建立燃耗矩阵时，不需要进行裂变产物份额的修正，也不需要TransMatrixFissionYields，其已经包含在xs中
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
//							ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
//	if (mode == 2) //定功率情况
//	{
//		if (lib_tag_ == 1) //lib_tag_ = 1意味着读取depth_library_name_
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
//						ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
//				ConstructFissionYieldsSpMat(); // 每个燃耗步调整裂变产物份额
//
//				ModecNuclideLibrary.CalculateFlux(1);
//			}
//			//ModecNuclideLibrary.CalculateFlux(mode);
//			flux_vector_.push_back(ModecNuclideLibrary.flux_);
//			power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//		}
//		else // 在读取couple文件建立燃耗矩阵时，不需要进行裂变产物份额的修正，也不需要TransMatrixFissionYields，其已经包含在xs中
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
//						ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
	power_vector_.push_back(ModecNuclideLibrary.specified_power_); // 初始时刻的通量和功率统计

	if (lib_tag_ == 1) //lib_tag_ = 1意味着读取depth_library_name_
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
				error_[F_i] = abs(F_mol[F_i] / n_vector_[1][F_i]); // 将添料率常数的贡献加入总的核素浓度中去
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
					n_vector_[1][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
				}
			}
		}
		else
		{
			for (int F_i = 0; F_i < size_F; ++F_i)
			{
				n_vector_[1][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
	power_vector_.push_back(ModecNuclideLibrary.specified_power_);// 平衡态的通量和功率统计
}

void ModecClass::Evolution(int mode, double time, int subtime)
{
	if (if_tracking_stockage == true)
	{
		TransMatrixStockage.SymbolLUElimination();
	}

	switch (mode)
	{
		case 0:	// 纯衰变情况
		{
			{
				if (solver_selection_ == 1)
				{
					if (if_constant_online_feeding_ == true && constant_feeding_calculation_methods_ == 2)
					{
						int size_matrix = TransMatrixDecay.spmat_dimen_;
						TransMatrixDecay.Resize(size_matrix + 1); // 矩阵增广
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
										ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
										ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
						// TTA方法的在线添料实现 
						// 增加一个核素，并给燃耗邻接矩阵增加一维
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

		case 1:	// 定通量情况
		{
			{
				if (lib_tag_ == 1) //lib_tag_ = 1意味着读取DepthLib
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

									ConstructFissionYieldsSpMat(); // 每个燃耗步调整裂变产物份额
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
									ConstructFissionYieldsSpMat(); // 每个燃耗步调整裂变产物份额
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
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
										}

									}
								}

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);

								ConstructFissionYieldsSpMat(); // 每个燃耗步调整裂变产物份额
							}

						}
					}
					else if (solver_selection_ == 0)
					{
						if (if_constant_online_feeding_ == false) {
							Solver.TtaInitialize(ModecNuclideLibrary.nuclide_library_vector_[0]);
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);

							double re_time = 0;
							for (int i = 1; i <= subtime; ++i)
							{
								re_time += time;
								TransMatrix = TtaMatrixDecay + (TtaMatrixCrossSection + TtaMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
								ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.TtaSolver(TransMatrix, re_time);

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);

								ConstructFissionYieldsSpMatForTta(); // 每个燃耗步调整裂变产物份额
							}
						}
						else
						{
							// TTA方法的在线添料实现 
							// 增加一个核素，并给燃耗邻接矩阵增加一维
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

							Solver.TtaInitialize(F_mol);
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);

							double re_time = 0;
							for (int i = 1; i <= subtime; ++i)
							{
								re_time += time;
								TransMatrix = TtaMatrixDecay + (TtaMatrixCrossSection + TtaMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24);
								F_mol = Solver.TtaSolverForFeeding(TransMatrix, re_time);

								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
								}

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);

								ConstructFissionYieldsSpMatForTta(); // 每个燃耗步调整裂变产物份额
							}
						}
					}

				}
				else // 在读取couple文件建立燃耗矩阵时，不需要进行裂变产物份额的修正，也不需要TransMatrixFissionYields，其已经包含在xs中
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
							//TransitionMatrixOutput(TransMatrixDecay*time); // 输出矩阵的各个元素到文件中

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
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
							Solver.TtaInitialize(ModecNuclideLibrary.nuclide_library_vector_[0]);
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);
							TransMatrix = TtaMatrixDecay + TtaMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);

							double re_time = 0;
							for (int i = 1; i <= subtime; ++i)
							{
								re_time += time;
								ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.TtaSolver(TransMatrix, re_time);

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
								TtaMatrixDecay.AddElementCCS(index, size_matrix, constant_feeding_rate_[i]);
								tot_feeding_rate += constant_feeding_rate_[i];
							}
							TtaMatrixDecay.AddElementCCS(size_matrix, size_matrix, -tot_feeding_rate);
							Solver.tot_feeding_rate_ = tot_feeding_rate;

							vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
							F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 0.0);

							Solver.TtaInitialize(F_mol);
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);
							TransMatrix = TtaMatrixDecay + TtaMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);

							double re_time = 0;
							for (int i = 1; i <= subtime; ++i)
							{
								re_time += time;
								F_mol = Solver.TtaSolver(TransMatrix, re_time);

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
	
		case 2:	//定功率情况
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
									ConstructFissionYieldsSpMat(); // 每个燃耗步调整裂变产物份额
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
									ConstructFissionYieldsSpMat(); // 每个燃耗步调整裂变产物份额
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
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
										}

									}
								}
								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
								ConstructFissionYieldsSpMat(); // 每个燃耗步调整裂变产物份额
							}

						}
					}
					else if (solver_selection_ == 0)
					{
						if (if_constant_online_feeding_ == false) {
							Solver.TtaInitialize(ModecNuclideLibrary.nuclide_library_vector_[0]);
							SparseMatrixMCS TransMatrix;
							ModecNuclideLibrary.CalculateFlux(mode);

							double re_time = 0;
							for (int i = 1; i <= subtime; ++i)
							{
								re_time += time;
								TransMatrix = (TtaMatrixCrossSection + TtaMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24) + TtaMatrixDecay;

								ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.TtaSolver(TransMatrix, re_time);

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
								ConstructFissionYieldsSpMatForTta(); // 每个燃耗步调整裂变产物份额
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
								TtaMatrixDecay.AddElementCCS(index, size_matrix, constant_feeding_rate_[i]);
								tot_feeding_rate += constant_feeding_rate_[i];
							}
							TtaMatrixDecay.AddElementCCS(size_matrix, size_matrix, -tot_feeding_rate);
							Solver.tot_feeding_rate_ = tot_feeding_rate;

							vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
							F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 0.0);

							Solver.TtaInitialize(F_mol);
							SparseMatrixMCS TransMatrix(size_matrix + 1);
							ModecNuclideLibrary.CalculateFlux(mode);

							double re_time = 0;
							for (int i = 1; i <= subtime; ++i)
							{
								re_time += time;
								TransMatrix = (TtaMatrixCrossSection + TtaMatrixFissionYields)*(ModecNuclideLibrary.flux_ * 1.0e-24) + TtaMatrixDecay;

								F_mol = Solver.TtaSolver(TransMatrix, re_time);

								for (unsigned int j = 0; j < ModecNuclideLibrary.nuclide_library_vector_[0].size(); ++j)
								{
									ModecNuclideLibrary.nuclide_library_vector_[0][j] = F_mol[j];
								}

								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);
								ConstructFissionYieldsSpMatForTta(); // 每个燃耗步调整裂变产物份额
							}
						}
					}

				}
				else// 在读取couple文件建立燃耗矩阵时，不需要进行裂变产物份额的修正，也不需要TransMatrixFissionYields，其已经包含在xs中
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
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
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
											ModecNuclideLibrary.nuclide_library_vector_[0][F_i] += F_mol[F_i]; // 将添料率常数的贡献加入总的核素浓度中去
										}

									}
								}
								n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

								ModecNuclideLibrary.CalculateFlux(mode);
								flux_vector_.push_back(ModecNuclideLibrary.flux_);
								power_vector_.push_back(ModecNuclideLibrary.specified_power_);

								//ConstructFissionYieldsSpMat(TransMatrixFissionYields, ModecNuclideLibrary); // 每个燃耗步调整裂变产物份额
							}

						}
					}
					else if (solver_selection_ == 0)
					{
						if (if_constant_online_feeding_ == false) {
							Solver.TtaInitialize(ModecNuclideLibrary.nuclide_library_vector_[0]);
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);


							double re_time = 0;
							for (int i = 1; i <= subtime; ++i)
							{
								re_time += time;

								TransMatrix = TtaMatrixDecay + TtaMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);
								ModecNuclideLibrary.nuclide_library_vector_[0] = Solver.TtaSolver(TransMatrix, re_time);

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
								TtaMatrixDecay.AddElementCCS(index, size_matrix, constant_feeding_rate_[i]);
								tot_feeding_rate += constant_feeding_rate_[i];
							}
							TtaMatrixDecay.AddElementCCS(size_matrix, size_matrix, -tot_feeding_rate);
							Solver.tot_feeding_rate_ = tot_feeding_rate;

							vector<double > F_mol(ModecNuclideLibrary.nuclide_library_vector_[0]);
							F_mol.resize(ModecNuclideLibrary.nuclide_library_vector_[0].size() + 1, 0.0);

							Solver.TtaInitialize(F_mol);
							SparseMatrixMCS TransMatrix(TtaMatrixDecay);
							ModecNuclideLibrary.CalculateFlux(mode);

							double re_time = 0;
							for (int i = 1; i <= subtime; ++i)
							{
								re_time += time;

								TransMatrix = TtaMatrixDecay + TtaMatrixCrossSection*(ModecNuclideLibrary.flux_ * 1.0e-24);
								F_mol = Solver.TtaSolver(TransMatrix, re_time);

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
	}
}