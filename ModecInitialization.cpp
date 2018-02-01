#include "ModecClass.h"

void ModecClass::ModecInitial(int argc, char *argv[])
{
	ifstream GetInputFileName;
	GetInputFileName.open("GetInputFile.burn");
	while (GetInputFileName && !GetInputFileName.eof())
	{
		getline(GetInputFileName, work_direc_);
		getline(GetInputFileName, input_filename_);
		if (work_direc_.length()>0 && work_direc_[work_direc_.length() - 1] != '/')
		{
			work_direc_ = work_direc_ + "/";
		}
		break;
	}
	GetInputFileName.close();

	int num_nucl;								// 初始核素个数 

	if (argc == 2)
	{
		input_filename_ = argv[1];
	}
	else if (argc > 2)
	{
		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: extra message in console.", 1);
	}

	input_file_ = work_direc_ + input_filename_;

input_mark:
	modec_inp_.open(input_file_);
	if (!modec_inp_.is_open())
	{
		cerr << "No Default MODEC input file (modec.input) !! " << '\n';
		cerr << "Please Enter User-Defined MODEC Input Filename or Enter '0' to exit: ";
		cin >> input_filename_;
		if (input_filename_ == "0")
		{
			exit(0);
		}
		input_file_ = work_direc_ + input_filename_;
		goto input_mark;
	}
	InfoMessage::start = clock();
	InfoMessage::InputName = input_filename_ + ".log";
	string initial_tag;
	modec_inp_ >> initial_tag;
	if (initial_tag != "=MODEC")
	{
		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no '=MODEC' tag in the top of input file.", 1);
	}

readfile:
	while (!modec_inp_.eof())
	{
		string line;
		getline(modec_inp_, line);
		istringstream record(line);
		string tag;
		record >> tag;

		/////////////////// Block 1 //////////////////////
		if (tag == "DepthLib")
		{
			record >> depth_library_name_;
			lib_tag_ = 1;
		}
		if (tag == "DecayLib")
		{
			record >> decay_library_name_;
			lib_tag_ = 0;
		}
		if (tag == "FYlib")
		{
			record >> fission_yields_library_name_;
			lib_tag_ = 0;
		}
		if (tag == "CoupleLib")
		{
			record >> couple_library_name_;
			lib_tag_ = 2;
		}
		if (tag == "CalculateEquilibrium")
		{
			record >> if_calculate_equilibrium_;
		}
		if (tag == "GL_order")
		{
			record >> GaussLegendreQuadrature::GL_order;
		}

		/*if (tag == "GL_order_eql")
		{
			int temp;
			record >> temp;
			GL_weight_eql = GaussLaguerreQuadrature::GLa_weight[temp];
			GL_abscissa_eql = GaussLaguerreQuadrature::GLa_abscissa[temp];
		}*/

		if (tag == "Solver")
		{
			record >> solver_selection_;
		}
		if (tag == "TTA_cutoff")
		{
			record >> Solver.cutoff_std_;
		}
		if (tag == "Solver_feed")
		{
			record >> constant_feeding_calculation_methods_;
		}
		if (tag == "Print")
		{
			record >> print_mode_;
		}
		if (tag == "Radioactivity")
		{
			record >> if_print_activity_;
		}
		if (tag == "DecayEnergy")
		{
			record >> if_print_decayenergy_;
		}
		if (tag == "AMPCtoxicity")
		{
			record >> if_print_ampc_;
		}
		if (tag == "WMPCtoxicity")
		{
			record >> if_print_wmpc_;
		}
		if (tag == "Svtoxicity")
		{
			record >> if_print_toxicity_;
		}
		if (tag == "Print_Kinf")
		{
			record >> if_print_kinf_;
		}
		if (tag == "Print_FissionRate")
		{
			record >> if_print_fission_rate_;
		}
		if (tag == "Print_AbsorptionRate")
		{
			record >> if_print_absorption_rate_;
		}
		//////////////////////////////////////////////////

		/////////////////// Block 2 //////////////////////
		if (tag == "Decay")
		{
			if (if_calculate_equilibrium_ == true)
			{
				goto readfile;
			}
			if_read_mode_tag_ = 1;
			double time;
			string ctemp;
			record >> time;
			record >> ctemp;
			time_unit_.push_back(ctemp);
			if (ctemp == "d" || ctemp == "day" || ctemp == "days")
			{
				time = time * 24 * 3600;
			}
			if (ctemp == "m" || ctemp == "month" || ctemp == "months")
			{
				time = time * 24 * 3600 * 30;
			}
			if (ctemp == "y" || ctemp == "year" || ctemp == "years")
			{
				time = time * 365.25 * 24 * 3600;
			}
			burnup_time_.push_back(time);
			record >> ctemp;
			int sub;
			record >> sub;
			substep_.push_back(sub);
			evolution_value_.push_back(0.0);
			evolution_mode_.push_back(0);
		}

		if (tag == "Flux")
		{
			if_read_mode_tag_ = 1;
			double flux;
			record >> flux;
			evolution_value_.push_back(flux);
			evolution_mode_.push_back(1);
			if (if_calculate_equilibrium_ == true)
			{
				goto readfile;
			}
			double time;
			string ctemp;
			record >> time;
			record >> ctemp;
			time_unit_.push_back(ctemp);
			if (ctemp == "d" || ctemp == "day" || ctemp == "days")
			{
				time = time * 24 * 3600;
			}
			if (ctemp == "m" || ctemp == "month" || ctemp == "months")
			{
				time = time * 24 * 3600 * 30;
			}
			if (ctemp == "y" || ctemp == "year" || ctemp == "years")
			{
				time = time * 365.25 * 24 * 3600;
			}
			burnup_time_.push_back(time);
			record >> ctemp;
			int sub;
			record >> sub;
			substep_.push_back(sub);		
		}
		
		if (tag == "Power")
		{
			if_read_mode_tag_ = 1;
			double power;
			record >> power;
			//power_vector_.push_back(ModecNuclideLibrary.specified_power_);
			evolution_value_.push_back(power);
			evolution_mode_.push_back(2);

			if (if_calculate_equilibrium_ == true)
			{
				goto readfile;
			}
			double time;
			record >> time;
			string ctemp;
			record >> ctemp;
			time_unit_.push_back(ctemp);
			if (ctemp == "d" || ctemp == "day" || ctemp == "days")
			{
				time = time * 24 * 3600;
			}
			if (ctemp == "m" || ctemp == "month" || ctemp == "months")
			{
				time = time * 24 * 3600 * 30;
			}
			if (ctemp == "y" || ctemp == "year" || ctemp == "years")
			{
				time = time * 365.25 * 24 * 3600;
			}
			burnup_time_.push_back(time);
			record >> ctemp;
			int sub;
			record >> sub;
			substep_.push_back(sub);
			
		}

		if (tag == "Calc_Flow_Depletion")
		{
			record >> if_flow_mode_;
		}
		if (tag == "Residue_Time")
		{
			double temp;
			record >> temp;
			residue_time_.push_back(temp);
			record >> temp;
			residue_time_.push_back(temp);
		}
		//////////////////////////////////////////////////

		/////////////////// Block 3 //////////////////////
		if (tag == "Num_Nucl")
		{
			record >> num_nucl;
		}
		if (tag == "Dens_unit")
		{
			record >> dens_unit_;
		}
		if (tag == "Nuclide")
		{
			if_read_density_tag_ = 1;
			//ModecNuclideLibrary.heavy_metal_mass_ = 0; // gram
			if (dens_unit_ == "mol") { goto mol; }
			else if (dens_unit_ == "g") { goto g; }
			else if (dens_unit_ == "kg") { goto kg; }
			else if (dens_unit_ == "atom") { goto atom; }
			else if (dens_unit_ == "atom/(barn-cm)") { goto atombarn; }
			else
			{
				InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: error density unit! only the fllowing units can be used: mol, g, kg, atom, atom/(barn-cm) !!", 1);
			}

		mol:
			for (int i = 0; i < num_nucl; ++i)
			{
				string nuclstr;
				getline(modec_inp_, nuclstr);
				istringstream record(nuclstr);
				int nucl_id;
				double n_mol;
				record >> nucl_id >> n_mol;
				int id = ModecNuclideLibrary.GetNuclIndex(nucl_id);
				ModecNuclideLibrary.nuclide_library_vector_[0][id] = n_mol;
				double awt = nucl_id - nucl_id / 10000 * 10000;
				if (awt >= 220)
				{
					ModecNuclideLibrary.heavy_metal_mass_ += n_mol * awt;
				}
			}
			goto readfile;

		g:
			for (int i = 0; i < num_nucl; ++i)
			{
				string nuclstr;
				getline(modec_inp_, nuclstr);
				istringstream record(nuclstr);
				int nucl_id;
				double mass;
				record >> nucl_id >> mass;
				int id = ModecNuclideLibrary.GetNuclIndex(nucl_id);
				double awt = double((nucl_id - nucl_id / 10000 * 10000) / 10);

				ModecNuclideLibrary.nuclide_library_vector_[0][id] = mass / awt; // 转换成mol进行存储
				
				if (awt >= 220)
				{
					ModecNuclideLibrary.heavy_metal_mass_ += mass;
				}
			}
			goto readfile;
		
		kg:
			for (int i = 0; i < num_nucl; ++i)
			{
				string nuclstr;
				getline(modec_inp_, nuclstr);
				istringstream record(nuclstr);
				int nucl_id;
				double mass;
				record >> nucl_id >> mass;
				int id = ModecNuclideLibrary.GetNuclIndex(nucl_id);
				double awt = double((nucl_id - nucl_id / 10000 * 10000) / 10);

				ModecNuclideLibrary.nuclide_library_vector_[0][id] = mass * 1000.0 / awt; // 转换成mol进行存储

				if (awt >= 220)
				{
					ModecNuclideLibrary.heavy_metal_mass_ += mass * 1000.0; // 转换成g来存储
				}
			}
			goto readfile;
			
		atom:
			for (int i = 0; i < num_nucl; ++i)
			{
				string nuclstr;
				getline(modec_inp_, nuclstr);
				istringstream record(nuclstr);
				int nucl_id;
				double atom;
				record >> nucl_id >> atom;
				int id = ModecNuclideLibrary.GetNuclIndex(nucl_id);
				double awt = nucl_id - nucl_id / 10000 * 10000;

				ModecNuclideLibrary.nuclide_library_vector_[0][id] = atom / (6.022140857E+23); // 转换成mol进行存储

				if (awt >= 220)
				{
					ModecNuclideLibrary.heavy_metal_mass_ += atom / (6.022140857E+23) * awt; // 转换成g来存储
				}
			}
			goto readfile;

		atombarn:
			for (int i = 0; i < num_nucl; ++i)
			{
				string nuclstr;
				getline(modec_inp_, nuclstr);
				istringstream record(nuclstr);
				int nucl_id;
				double atom;
				record >> nucl_id >> atom;
				int id = ModecNuclideLibrary.GetNuclIndex(nucl_id);
				double awt = nucl_id - nucl_id / 10000 * 10000;

				ModecNuclideLibrary.nuclide_library_vector_[0][id] = atom / (6.022140857E-1); // 转换成mol进行存储

				if (awt >= 220)
				{
					ModecNuclideLibrary.heavy_metal_mass_ += atom / (6.022140857E-1) * awt; // 转换成g来存储
				}
			}
			goto readfile;
		}
		//////////////////////////////////////////////////

		/////////////////// Block 4 //////////////////////
		if (tag == "OnlineReprocessing")
		{
			record >> if_continously_remove_;
			if (if_continously_remove_ == true)
			{
				string removeline;
				getline(modec_inp_, removeline);
				istringstream record(removeline);
				string removetag;
				record >> removetag;
				if (removetag != "Ele_GroupNum")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Ele_GroupNum' tag in the input file.", 1);
				}
				record >> remove_group_number_;
				remove_group_vector_.resize(remove_group_number_);
				remove_rate_vector_.resize(remove_group_number_);
				remove_element_vector_.resize(remove_group_number_);

				removeline.clear();
				record.clear();

				getline(modec_inp_, removeline);
				record.str(removeline);
				record >> removetag;
				if (removetag != "Ele_Group")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Ele_Group' tag in the input file.", 1);
				}
				for (int i = 0; i < remove_group_number_; ++i)
				{
					record >> remove_group_vector_[i];
				}
				removeline.clear();
				record.clear();

				getline(modec_inp_, removeline);
				record.str(removeline);
				record >> removetag;
				if (removetag != "Ele_RemoveRate")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Ele_RemoveRate' tag in the input file.", 1);
				}
				for (int i = 0; i < remove_group_number_; ++i)
				{
					record >> remove_rate_vector_[i];
				}
				removeline.clear();
				record.clear();

				for (int i = 0; i < remove_group_number_; ++i)
				{
					remove_element_vector_[i].resize(remove_group_vector_[i]);
				}
				modec_inp_ >> removetag;
				if (removetag != "Ele_ID")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Ele_ID' tag in the input file.", 1);
				}
				for (int i = 0; i < remove_group_number_; ++i)
				{
					int size = remove_group_vector_[i];
					for (int j = 0; j < size; ++j)
					{
						modec_inp_ >> remove_element_vector_[i][j];
					}
				}
				modec_inp_ >> removetag;
				if (removetag != "TrackingStokage")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'TrackingStokage' tag in the input file.", 1);
				}
				modec_inp_ >> if_tracking_stockage;
				if (if_tracking_stockage == true)
				{
					int size = ModecNuclideLibrary.nuclide_library_vector_[0].size();
					ModecNuclideLibrary.nuclide_library_vector_[0].resize(2 * size); // 若要追踪堆外核素演化，则矩阵中核素总数应为原来的两倍
				}
				modec_inp_ >> removetag;
				if (removetag != "StokageRadioactivity")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'StokageRadioactivity' tag in the input file.", 1);
				}
				modec_inp_ >> if_print_stockage_activity_;

				modec_inp_ >> removetag;
				if (removetag != "StokageDecayEnergy")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'StokageDecayEnergy' tag in the input file.", 1);
				}
				modec_inp_ >> if_print_stockage_decayenergy_;

				modec_inp_ >> removetag;
				if (removetag != "StokageAMPCtoxicity")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'StokageAMPCtoxicity' tag in the input file.", 1);
				}
				modec_inp_ >> if_print_stockage_ampc_;

				modec_inp_ >> removetag;
				if (removetag != "StokageWMPCtoxicity")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'StokageWMPCtoxicity' tag in the input file.", 1);
				}
				modec_inp_ >> if_print_stockage_wmpc_;

				modec_inp_ >> removetag;
				if (removetag != "StokageSvtoxicity")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'StokageSvtoxicity' tag in the input file.", 1);
				}
				modec_inp_ >> if_print_stockage_toxicity_;
			}
		}
		//////////////////////////////////////////////////

		if (tag == "KeepingEutecticStable")
		{
			record >> if_keeping_eutectic_stable_;
		}

		/////////////////// Block 4 //////////////////////
		if (tag == "ContinuouslyFeeding")
		{
			record >> if_variable_feeding_;
			if (if_variable_feeding_ == true)
			{
				if (if_constant_online_feeding_ == true)
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: the variables of 'if_variable_feeding_' and 'if_constant_online_feeding_' cannot be both true.", 1);
				}
				string feedingline;
				getline(modec_inp_, feedingline);
				istringstream record(feedingline);
				string feedingtag;
				record >> feedingtag;
				if (feedingtag != "Nucl_GroupNum")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_GroupNum' tag in the input file.", 1);
				}
				record >> variable_feeding_group_num_;
				variable_feeding_group_vector_.resize(variable_feeding_group_num_);
				variable_feeding_nuclide_id_vector_.resize(variable_feeding_group_num_);
				variable_feeding_nuclide_ratio_vector_.resize(variable_feeding_group_num_);
				feedingline.clear();
				record.clear();

				getline(modec_inp_, feedingline);
				record.str(feedingline);
				record >> feedingtag;
				if (feedingtag != "Nucl_Ratio")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_Ratio' tag in the input file.", 1);
				}
				record >> variable_feeding_ratio_;
				feedingline.clear();
				record.clear();

				getline(modec_inp_, feedingline);
				record.str(feedingline);
				record >> feedingtag;
				if (feedingtag != "Nucl_Group")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_Group' tag in the input file.", 1);
				}
				for (int i = 0; i < variable_feeding_group_num_; ++i)
				{
					record >> variable_feeding_group_vector_[i];
				}
				feedingline.clear();
				record.clear();

				for (int i = 0; i < variable_feeding_group_num_; ++i)
				{
					variable_feeding_nuclide_id_vector_[i].resize(variable_feeding_group_vector_[i]);
					variable_feeding_nuclide_ratio_vector_[i].resize(variable_feeding_group_vector_[i]);
				}
				modec_inp_ >> feedingtag;
				if (feedingtag != "Nucl_ID")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_ID' tag in the input file.", 1);
				}
				for (int i = 0; i < variable_feeding_group_num_; ++i)
				{
					int size = variable_feeding_group_vector_[i];
					for (int j = 0; j < size; ++j)
					{
						modec_inp_ >> variable_feeding_nuclide_id_vector_[i][j];
						modec_inp_ >> variable_feeding_nuclide_ratio_vector_[i][j];
					}
				}
			}
		}
		//////////////////////////////////////////////////
		if (tag == "ConstantContinuouslyFeeding")
		{
			record >> if_constant_online_feeding_;
			if (if_constant_online_feeding_ == true)
			{
				gauss_legendre_weight_ = GaussLegendreQuadrature::gauss_legendre_weight_[GaussLegendreQuadrature::GL_order];
				gauss_legendre_abscissa_ = GaussLegendreQuadrature::gauss_legendre_abscissa_[GaussLegendreQuadrature::GL_order];
				if (if_variable_feeding_ == true)
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: the variables of 'if_variable_feeding_' and 'if_constant_online_feeding_' cannot be both true.", 1);
				}
				string feedingline;
				getline(modec_inp_, feedingline);
				istringstream record(feedingline);
				string feedingtag;
				record >> feedingtag;
				if (feedingtag != "Nucl_Num")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_Num' tag in the input file.", 1);
				}
				record >> constant_feeding_nuclide_num_;
				constant_feeding_nuclide_id_vector_.resize(constant_feeding_nuclide_num_);
				constant_feeding_rate_.resize(constant_feeding_nuclide_num_);
				feedingline.clear();
				record.clear();

				getline(modec_inp_, feedingline);
				record.str(feedingline);
				record >> feedingtag;
				if (feedingtag != "Nucl_ID")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_ID' tag in the input file.", 1);
				}
				for (int i = 0; i < constant_feeding_nuclide_num_; ++i)
				{
					record >> constant_feeding_nuclide_id_vector_[i];
				}
				feedingline.clear();
				record.clear();

				getline(modec_inp_, feedingline);
				record.str(feedingline);
				record >> feedingtag;
				if (feedingtag != "Nucl_Rate")
				{
					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_Rate' tag in the input file.", 1);
				}
				for (int i = 0; i < constant_feeding_nuclide_num_; ++i)
				{
					record >> constant_feeding_rate_[i];
				}
				feedingline.clear();
				record.clear();
			}
		}
	}
	modec_inp_.close();

	if (if_flow_mode_ == 1)
	{
		if(residue_time_.size() ==0)
			InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no residue time data are read.", 1);

		if (residue_time_[0] == 0 || residue_time_[1] == 0)
		{
			if_flow_mode_ = 0;
			InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Warning: zero residue time means totally well-mixed.",0);
			for (int i = 0; i < evolution_value_.size(); ++i)
			{
				evolution_value_[i] /= 2.0; // 完全均匀则意味着减半
			}
		}
		
		if(solver_selection_ == 0)
		{
			InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Warning: CRAM method must be used when flow depletion is calculated.",0);
			solver_selection_ == 1;
		}
	}
	
	if (decay_library_name_.length() == 0 && depth_library_name_.length() == 0 && couple_library_name_.length() == 0)
	{
		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no library name in the input file.", 1);
	}
	if (if_read_mode_tag_ == 0)
	{
		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no depletion model setting in the input file. Depletion model can be set as one or combined of the followings: Decay, Flux, Power.", 1);
	}
	if (if_read_density_tag_ == 0)
	{
		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: error of getting nuclide densities.", 1);
	}

	int size_mode = evolution_mode_.size();
	if (if_calculate_equilibrium_ == true && size_mode != 1)
	{
		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: only one depletion model can be set in equilibrium calculation. Depletion model can be set as one of the followings: Decay, Flux, Power.", 1);
	}
	if (if_calculate_equilibrium_ == true && if_tracking_stockage == true)
	{
		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: cannot track nuclide evolution outside the core in equilibrium calculation, the tag 'TrackingStokage' must be set false.", 0);
		if_tracking_stockage = false;
	}

	n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);

	if (if_calculate_equilibrium_ == true)
	{
		n_vector_.resize(2);
		n_vector_[1].resize(ModecNuclideLibrary.nuclide_library_vector_[0].size());
	}

	if (if_constant_online_feeding_ == true && if_tracking_stockage == false && constant_feeding_calculation_methods_ == 1)
	{
		constant_feeding_vector_.resize(ModecNuclideLibrary.nuclide_number_);
		int size = constant_feeding_nuclide_id_vector_.size();
		for (int i = 0; i < size; ++i)
		{
			int Index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
			constant_feeding_vector_[Index] = constant_feeding_rate_[i];
		}
	}
	if (if_constant_online_feeding_ == true && if_tracking_stockage == true && constant_feeding_calculation_methods_ == 1)
	{
		constant_feeding_vector_.resize(ModecNuclideLibrary.nuclide_number_*2);
		int size = constant_feeding_nuclide_id_vector_.size();
		for (int i = 0; i < size; ++i)
		{
			int Index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
			constant_feeding_vector_[Index] = constant_feeding_rate_[i];
		}
	}
};

void ModecClass::BuildSpMat()
{
	if (lib_tag_ == 0)
	{
		int size = evolution_mode_.size();
		int i;
		for (i = 0; i < size; ++i)
		{
			if (evolution_mode_[i] > 0)
			{
				break;
			}
		}
		if (solver_selection_ == 1)
		{
			DecayToSpMat();
		}
		else if(solver_selection_ == 0)
		{
			DecayToSpMatForTta();
		}
		if (i < size) // 耦合TRITON
		{
			XSfromTriton();
			CalculateEffectiveFissionYields();
		}
	}
	else if (lib_tag_ == 1) // 读取depth_library_name_
	{
		int size = evolution_mode_.size();
		int i;
		for (i = 0; i < size; ++i)
		{
			if (evolution_mode_[i] > 0)
			{
				break;
			}
		}
		if (i >= size) // 纯衰变情形
		{
			if (solver_selection_ == 1)
			{
				ReadFromDepthLib();
			}
			else if (solver_selection_ == 0)
			{
				ReadFromDepthLibForTta();
			}
		}
		else
		{
			if (solver_selection_ == 1)
			{
				ReadFromDepthLib();
				ConstructFissionYieldsSpMat();
			}
			else if (solver_selection_ == 0)
			{
				ReadFromDepthLibForTta();
				ConstructFissionYieldsSpMatForTta();
			}
		}
	}
	else if (lib_tag_ == 2)// 读取Couple
	{
		int size = evolution_mode_.size();
		int i;
		for (i = 0; i < size; ++i)
		{
			if (evolution_mode_[i] > 0)
			{
				break;
			}
		}
		if (solver_selection_ == 1)
		{
			ReadFromCouple();
		}
		else if( solver_selection_ == 0 )
		{
			ReadFromCoupleForTta();
		}
	}

	if (if_flow_mode_ == 1)
	{
		TransMatrixPureDecay = TransMatrixDecay; // 定义纯衰变矩阵
	}

 	if (if_continously_remove_ == true)
	{
		AddOnlineReprocessingCoeffi();
	}
	if (if_variable_feeding_ == true)
	{
		ContinuouslyFeeding();
	}
}
