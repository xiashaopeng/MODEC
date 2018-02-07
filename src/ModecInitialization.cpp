#include "ModecClass.h"

using namespace tinyxml2;

void ModecClass::ModecInitial(int argc, char **argv) {
    if (argc == 2) {
        input_filename_ = argv[1];
    } else if (argc > 2) {
        InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: extra message in console.", 1);
    }
    input_file_ = work_direc_ + input_filename_;

    XMLDocument InputFile;

    while (InputFile.LoadFile(input_file_.c_str())) {
        int error_id = InputFile.ErrorID();
        if(error_id == 0)
            break;

        if(error_id != 3) {
            info_message_ ="ERROR OCCURRED WHEN INPUT FILE WAS BEING READ!\nError MESSAGE: " + string(InputFile.ErrorName());
            InfoMessage::ErrorMessage(info_message_ , 1);
        }

        cerr << "NO DEFAULT MODEC INPUT FILE!! " << '\n';
        cerr << "ENTER USER-DEFINED MODEC INPUT FILENAME OR ENTER '0' TO EXIT: ";
        cin >> input_filename_;
        if (input_filename_ == "0") {
            exit(0);
        }
        input_file_ = work_direc_ + input_filename_;
    }


//input_mark:
//	InputFile.LoadFile(input_file_.c_str());
//	int error_id = InputFile.ErrorID();
//	if(error_id != 0)
//	{
//		if(error_id != 3)
//		{
//			info_message_ ="ERROR OCCURRED WHEN INPUT FILE WAS BEING READ!\nError MESSAGE: " + string(InputFile.ErrorName());
//			InfoMessage::ErrorMessage(info_message_ , 1);
//		}
//
//		cerr << "NO DEFAULT MODEC INPUT FILE!! " << '\n';
//		cerr << "ENTER USER-DEFINED MODEC INPUT FILENAME OR ENTER '0' TO EXIT: ";
//		cin >> input_filename_;
//		if (input_filename_ == "0")
//		{
//			exit(0);
//		}
//		input_file_ = work_direc_ + input_filename_;
//		goto input_mark;
//	}

    InfoMessage::start = clock();
    InfoMessage::InputName = input_filename_ + ".log";

    XMLElement *modec = InputFile.RootElement();
    if( string(modec->Name()) != "MODEC") {
        InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'MODEC' tag in the top of input file.", 1);
    }

    /* ------------------------------------------------------------------------------------*
    				Mother Node "Settings"
       ------------------------------------------------------------------------------------
       |	Daughter Node	|			Attributes
       ------------------------------------------------------------------------------------
       |	"DataLib"	|	"type" "file"
       |	"SolverSet"	|	"solver_basic" "solver_nonhom"
       |	"FlowSet"	|	"tag" "residue_time"
       |	"PrintSet"	|	"print_den" "print_act" "print_q" "print_ampc"
       |			|	"print_wmpc" "print_tox" "print_kinf" "print_prod" "print_abs"
       ------------------------------------------------------------------------------------ */
    XMLElement * settings = modec->FirstChildElement();
    if (string(settings->Name()) != "Settings") {
        InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Settings' tag in the input file" ,1);
    }

    //	燃耗数据库设置
    XMLElement * data_lib = settings->FirstChildElement();
    string data_lib_type = data_lib->Attribute("type");
    string lib_file = data_lib->Attribute("file");
    if (data_lib_type == "DepthLib") {
        depth_library_name_ = lib_file;
        lib_tag_ = 1;
    }
    if (data_lib_type == "DecayLib") {
        decay_library_name_ = lib_file;
        lib_tag_ = 0;
    }
    if (data_lib_type == "FYlib") {
        fission_yields_library_name_ = lib_file;
        lib_tag_ = 0;
    }
    if (data_lib_type == "CoupleLib") {
        couple_library_name_ = lib_file;
        lib_tag_ = 2;
    }

    //	求解器设置
    XMLElement *solver_set = data_lib->NextSiblingElement();
    string solver_basic = solver_set->Attribute("solver_basic");
    if (solver_basic == "CRAM") {
        solver_selection_ = 1;
    } else if (solver_basic.substr(0,3) == "TTA") {
        solver_selection_ = 0;
        Solver.cutoff_std_ = atof(solver_basic.substr(3).c_str());
    }

    //	流动燃耗设置
    XMLElement *flow_set = solver_set->NextSiblingElement();
    if_flow_mode_ = atoi(flow_set->Attribute("tag"));
    if (if_flow_mode_ == 1) {
        if (solver_basic != "CRAM") {
            InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Warning: CRAM method must be used when flow depletion is calculated.",0);
            solver_selection_ = 1;
        }
        string residue_time_s = flow_set->Attribute("residue_time");
        stringstream residue_time_ss(residue_time_s);

        double temp1,temp2;
        residue_time_ss >> temp1 >> temp2;
        residue_time_.push_back(temp1);
        residue_time_.push_back(temp2);
    }

    //	打印设置
    XMLElement *print_set = flow_set->NextSiblingElement();
    print_mode_ = atoi(print_set->Attribute("print_den"));
    if_print_activity_ = atoi(print_set->Attribute("print_act"));
    if_print_decayenergy_ = atoi(print_set->Attribute("print_q"));
    if_print_ampc_ = atoi(print_set->Attribute("print_ampc"));
    if_print_wmpc_ = atoi(print_set->Attribute("print_wmpc"));
    if_print_toxicity_ = atoi(print_set->Attribute("print_tox"));
    if_print_kinf_ = atoi(print_set->Attribute("print_kinf"));
    if_print_fission_rate_ = atoi(print_set->Attribute("print_prod"));
    if_print_absorption_rate_ = atoi(print_set->Attribute("print_abs"));


    /* ------------------------------------------------------------------------------------*
    				Mother Node "Depletions"
       ------------------------------------------------------------------------------------
       |	Daughter Node	|			Attributes
       ------------------------------------------------------------------------------------
       |	"Depletion"	|	"mode" "value" "time" "step"
       ------------------------------------------------------------------------------------ */
    XMLElement *depletions = settings->NextSiblingElement();
    if (string(depletions->Name()) != "Depletions") {
        InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Depletions' tag in the input file" ,1);
    }
    XMLElement *depletion = depletions->FirstChildElement();
    while(depletion) {
        string mode = depletion->Attribute("mode");
        if (mode == "decay") {
            evolution_mode_.push_back(0);
            evolution_value_.push_back(0.0);
        } else if (mode == "flux") {
            evolution_mode_.push_back(1);
            evolution_value_.push_back(atof(depletion->Attribute("value")));
        } else if (mode == "power") {
            evolution_mode_.push_back(2);
            evolution_value_.push_back(atof(depletion->Attribute("value")));
        }

        stringstream time_ss(depletion->Attribute("time"));
        double time;
        string time_unit;
        time_ss >> time;
        time_ss >> time_unit;

        time_unit_.push_back(time_unit);
        if (time_unit == "d" || time_unit == "day" || time_unit == "days")
            time = time * 24.0 * 3600.0;
        if (time_unit == "m" || time_unit == "month" || time_unit == "months")
            time = time * 24.0 * 3600.0 * 30.0;
        if (time_unit == "y" || time_unit == "year" || time_unit == "years")
            time = time * 365.25 * 24.0 * 3600.0;
        burnup_time_.push_back(time);

        substep_.push_back(atoi(depletion->Attribute("step")));

        depletion = depletion->NextSiblingElement();
    }


    if( (if_flow_mode_ == 1) && (residue_time_[0] == 0 || residue_time_[1] == 0)) {
        if_flow_mode_ = 0;
        InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Warning: zero residue time means totally well-mixed.",0);
        for (int i = 0; i < evolution_value_.size(); ++i)
            evolution_value_[i] /= 2.0;
    }

    /* ------------------------------------------------------------------------------------*
    			Mother Node "Nuclides", attributes "units"
       ------------------------------------------------------------------------------------
       |	Daughter Node	|			Attributes
       ------------------------------------------------------------------------------------
       |	"Nuclide"	|		     "zai" "density"
       ------------------------------------------------------------------------------------ */
    XMLElement *nuclides = depletions->NextSiblingElement();
    if (string(nuclides->Name()) != "Nuclides") {
        InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nuclides' tag in the input file" ,1);
    }
    dens_unit_ = nuclides->Attribute("units");

    XMLElement *nuclide = nuclides->FirstChildElement();
    while (nuclide) {
        int nucl_id = atoi(nuclide->Attribute("zai"));
        double awt = double((nucl_id - nucl_id / 10000 * 10000) / 10);

        double coeff;
        if (dens_unit_ == "mol") coeff = 1.0;
        if (dens_unit_ == "g") coeff = 1.0 / awt;
        if (dens_unit_ == "kg") coeff = 1000.0/awt;
        if (dens_unit_ == "atom") coeff = 1.0 / (6.022140857E+23);
        if (dens_unit_ == "atom/(barn-cm)") coeff = 1.0 / 0.6022140857;
        ModecNuclideLibrary.nuclide_library_vector_[0][ModecNuclideLibrary.GetNuclIndex(nucl_id)]
            = atof(nuclide->Attribute("density")) * coeff;

        nuclide = nuclide->NextSiblingElement();
    }


    /* ------------------------------------------------------------------------------------*
    		Mother Node "OnlineReprocessing", attributes "tag" "track_storage"
       ------------------------------------------------------------------------------------
       |	Daughter Node	|			Attributes
       ------------------------------------------------------------------------------------
       |	"RemoveElement"	|		     	  "coeff"
       ------------------------------------------------------------------------------------ */
    XMLElement * online_reprocessing = nuclides->NextSiblingElement();
    if (string(online_reprocessing->Name()) != "OnlineReprocessing") {
        InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'OnlineReprocessing' tag in the input file" ,1);
    }
    if_continously_remove_ = atoi(online_reprocessing->Attribute("tag"));
    if (if_continously_remove_ == true) {
        if_tracking_stockage = atoi(online_reprocessing->Attribute("track_storage"));
        if (if_tracking_stockage == true) {
            int size = ModecNuclideLibrary.nuclide_library_vector_[0].size();
            ModecNuclideLibrary.nuclide_library_vector_[0].resize(2 * size);

            if_print_stockage_activity_ = if_print_activity_;
            if_print_stockage_decayenergy_ = if_print_decayenergy_;
            if_print_stockage_ampc_ = if_print_ampc_;
            if_print_stockage_wmpc_ = if_print_wmpc_;
            if_print_stockage_toxicity_ = if_print_toxicity_;
        }

        XMLElement * remove_element = online_reprocessing->FirstChildElement();
        while (remove_element) {
            remove_rate_vector_.push_back(atof(remove_element->Attribute("coeff")));

            vector< int > remove_element_vector;
            stringstream remove_ele_ss(remove_element->GetText());
            int temp_ele;
            while (remove_ele_ss >> temp_ele)
                remove_element_vector.push_back(temp_ele);
            remove_element_vector_.push_back(remove_element_vector);
            remove_group_vector_.push_back(remove_element_vector.size());

            remove_element = remove_element->NextSiblingElement();
        }
        remove_group_number_ = remove_group_vector_.size();
    }


    /* ------------------------------------------------------------------------------------*
    		Mother Node "ContinouslyFeeding", attributes "tag"
       ------------------------------------------------------------------------------------
       |	Daughter Node	|			Attributes
       ------------------------------------------------------------------------------------
       |	"FeedNuclide"	|		     "zai" "feed_rate"
       ------------------------------------------------------------------------------------ */
    XMLElement * continously_feeding = online_reprocessing->NextSiblingElement();
    if (string(continously_feeding->Name()) != "ContinouslyFeeding") {
        InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'ContinouslyFeeding' tag in the input file" ,1);
    }
    if_constant_online_feeding_ = atoi(continously_feeding->Attribute("tag"));
    if (if_constant_online_feeding_ == true) {
        stringstream method(continously_feeding->Attribute("method"));
        string method_name;
        method >> method_name;
        if( method_name == "Gauss") {
            constant_feeding_calculation_methods_ = 1;
            method >> GaussLegendreQuadrature::GL_order;
            gauss_legendre_weight_ = GaussLegendreQuadrature::gauss_legendre_weight_[GaussLegendreQuadrature::GL_order];
            gauss_legendre_abscissa_ = GaussLegendreQuadrature::gauss_legendre_abscissa_[GaussLegendreQuadrature::GL_order];
        }
        XMLElement * feed_nuclide = continously_feeding->FirstChildElement();
        while (feed_nuclide) {
            constant_feeding_nuclide_id_vector_.push_back(atoi(feed_nuclide->Attribute("zai")));
            constant_feeding_rate_.push_back(atof(feed_nuclide->Attribute("feed_rate")));
            feed_nuclide = feed_nuclide->NextSiblingElement();
        }
        constant_feeding_nuclide_num_ = constant_feeding_rate_.size();
    }


    /*---------------------------------------------------------------------------------------------------------------------------------------------*/

    if (if_constant_online_feeding_ == true && if_tracking_stockage == false && constant_feeding_calculation_methods_ == 1) {
        constant_feeding_vector_.resize(ModecNuclideLibrary.nuclide_number_);
        int size = constant_feeding_nuclide_id_vector_.size();
        for (int i = 0; i < size; ++i) {
            int Index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
            constant_feeding_vector_[Index] = constant_feeding_rate_[i];
        }
    }
    if (if_constant_online_feeding_ == true && if_tracking_stockage == true && constant_feeding_calculation_methods_ == 1) {
        constant_feeding_vector_.resize(ModecNuclideLibrary.nuclide_number_*2);
        int size = constant_feeding_nuclide_id_vector_.size();
        for (int i = 0; i < size; ++i) {
            int Index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
            constant_feeding_vector_[Index] = constant_feeding_rate_[i];
        }
    }

    n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
}

//void ModecClass::ModecInitial(int argc, char *argv[])
//{
//	ifstream GetInputFileName;
//	GetInputFileName.open("GetInputFile.burn");
//	while (GetInputFileName && !GetInputFileName.eof())
//	{
//		getline(GetInputFileName, work_direc_);
//		getline(GetInputFileName, input_filename_);
//		if (work_direc_.length()>0 && work_direc_[work_direc_.length() - 1] != '/')
//		{
//			work_direc_ = work_direc_ + "/";
//		}
//		break;
//	}
//	GetInputFileName.close();
//
//	int num_nucl;								// 初始核素个数
//
//	if (argc == 2)
//	{
//		input_filename_ = argv[1];
//	}
//	else if (argc > 2)
//	{
//		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: extra message in console.", 1);
//	}
//
//	input_file_ = work_direc_ + input_filename_;
//
//input_mark:
//	modec_inp_.open(input_file_);
//	if (!modec_inp_.is_open())
//	{
//		cerr << "No Default MODEC input file (modec.input) !! " << '\n';
//		cerr << "Please Enter User-Defined MODEC Input Filename or Enter '0' to exit: ";
//		cin >> input_filename_;
//		if (input_filename_ == "0")
//		{
//			exit(0);
//		}
//		input_file_ = work_direc_ + input_filename_;
//		goto input_mark;
//	}
//	InfoMessage::start = clock();
//	InfoMessage::InputName = input_filename_ + ".log";
//	string initial_tag;
//	modec_inp_ >> initial_tag;
//	if (initial_tag != "=MODEC")
//	{
//		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no '=MODEC' tag in the top of input file.", 1);
//	}
//
//readfile:
//	while (!modec_inp_.eof())
//	{
//		string line;
//		getline(modec_inp_, line);
//		istringstream record(line);
//		string tag;
//		record >> tag;
//
//		/////////////////// Block 1 //////////////////////
//		if (tag == "DepthLib")
//		{
//			record >> depth_library_name_;
//			lib_tag_ = 1;
//		}
//		if (tag == "DecayLib")
//		{
//			record >> decay_library_name_;
//			lib_tag_ = 0;
//		}
//		if (tag == "FYlib")
//		{
//			record >> fission_yields_library_name_;
//			lib_tag_ = 0;
//		}
//		if (tag == "CoupleLib")
//		{
//			record >> couple_library_name_;
//			lib_tag_ = 2;
//		}
//		if (tag == "CalculateEquilibrium")
//		{
//			record >> if_calculate_equilibrium_;
//		}
//		if (tag == "GL_order")
//		{
//			record >> GaussLegendreQuadrature::GL_order;
//		}
//
//		/*if (tag == "GL_order_eql")
//		{
//			int temp;
//			record >> temp;
//			GL_weight_eql = GaussLaguerreQuadrature::GLa_weight[temp];
//			GL_abscissa_eql = GaussLaguerreQuadrature::GLa_abscissa[temp];
//		}*/
//
//		if (tag == "Solver")
//		{
//			record >> solver_selection_;
//		}
//		if (tag == "TTA_cutoff")
//		{
//			record >> Solver.cutoff_std_;
//		}
//		if (tag == "Solver_feed")
//		{
//			record >> constant_feeding_calculation_methods_;
//		}
//		if (tag == "Print")
//		{
//			record >> print_mode_;
//		}
//		if (tag == "Radioactivity")
//		{
//			record >> if_print_activity_;
//		}
//		if (tag == "DecayEnergy")
//		{
//			record >> if_print_decayenergy_;
//		}
//		if (tag == "AMPCtoxicity")
//		{
//			record >> if_print_ampc_;
//		}
//		if (tag == "WMPCtoxicity")
//		{
//			record >> if_print_wmpc_;
//		}
//		if (tag == "Svtoxicity")
//		{
//			record >> if_print_toxicity_;
//		}
//		if (tag == "Print_Kinf")
//		{
//			record >> if_print_kinf_;
//		}
//		if (tag == "Print_FissionRate")
//		{
//			record >> if_print_fission_rate_;
//		}
//		if (tag == "Print_AbsorptionRate")
//		{
//			record >> if_print_absorption_rate_;
//		}
//		//////////////////////////////////////////////////
//
//		/////////////////// Block 2 //////////////////////
//		if (tag == "Decay")
//		{
//			if (if_calculate_equilibrium_ == true)
//			{
//				goto readfile;
//			}
//			if_read_mode_tag_ = 1;
//			double time;
//			string ctemp;
//			record >> time;
//			record >> ctemp;
//			time_unit_.push_back(ctemp);
//			if (ctemp == "d" || ctemp == "day" || ctemp == "days")
//			{
//				time = time * 24 * 3600;
//			}
//			if (ctemp == "m" || ctemp == "month" || ctemp == "months")
//			{
//				time = time * 24 * 3600 * 30;
//			}
//			if (ctemp == "y" || ctemp == "year" || ctemp == "years")
//			{
//				time = time * 365.25 * 24 * 3600;
//			}
//			burnup_time_.push_back(time);
//			record >> ctemp;
//			int sub;
//			record >> sub;
//			substep_.push_back(sub);
//			evolution_value_.push_back(0.0);
//			evolution_mode_.push_back(0);
//		}
//
//		if (tag == "Flux")
//		{
//			if_read_mode_tag_ = 1;
//			double flux;
//			record >> flux;
//			evolution_value_.push_back(flux);
//			evolution_mode_.push_back(1);
//			if (if_calculate_equilibrium_ == true)
//			{
//				goto readfile;
//			}
//			double time;
//			string ctemp;
//			record >> time;
//			record >> ctemp;
//			time_unit_.push_back(ctemp);
//			if (ctemp == "d" || ctemp == "day" || ctemp == "days")
//			{
//				time = time * 24 * 3600;
//			}
//			if (ctemp == "m" || ctemp == "month" || ctemp == "months")
//			{
//				time = time * 24 * 3600 * 30;
//			}
//			if (ctemp == "y" || ctemp == "year" || ctemp == "years")
//			{
//				time = time * 365.25 * 24 * 3600;
//			}
//			burnup_time_.push_back(time);
//			record >> ctemp;
//			int sub;
//			record >> sub;
//			substep_.push_back(sub);
//		}
//
//		if (tag == "Power")
//		{
//			if_read_mode_tag_ = 1;
//			double power;
//			record >> power;
//			//power_vector_.push_back(ModecNuclideLibrary.specified_power_);
//			evolution_value_.push_back(power);
//			evolution_mode_.push_back(2);
//
//			if (if_calculate_equilibrium_ == true)
//			{
//				goto readfile;
//			}
//			double time;
//			record >> time;
//			string ctemp;
//			record >> ctemp;
//			time_unit_.push_back(ctemp);
//			if (ctemp == "d" || ctemp == "day" || ctemp == "days")
//			{
//				time = time * 24 * 3600;
//			}
//			if (ctemp == "m" || ctemp == "month" || ctemp == "months")
//			{
//				time = time * 24 * 3600 * 30;
//			}
//			if (ctemp == "y" || ctemp == "year" || ctemp == "years")
//			{
//				time = time * 365.25 * 24 * 3600;
//			}
//			burnup_time_.push_back(time);
//			record >> ctemp;
//			int sub;
//			record >> sub;
//			substep_.push_back(sub);
//
//		}
//
//		if (tag == "Calc_Flow_Depletion")
//		{
//			record >> if_flow_mode_;
//		}
//		if (tag == "Residue_Time")
//		{
//			double temp;
//			record >> temp;
//			residue_time_.push_back(temp);
//			record >> temp;
//			residue_time_.push_back(temp);
//		}
//		//////////////////////////////////////////////////
//
//		/////////////////// Block 3 //////////////////////
//		if (tag == "Num_Nucl")
//		{
//			record >> num_nucl;
//		}
//		if (tag == "Dens_unit")
//		{
//			record >> dens_unit_;
//		}
//		if (tag == "Nuclide")
//		{
//			if_read_density_tag_ = 1;
//			//ModecNuclideLibrary.heavy_metal_mass_ = 0; // gram
//			if (dens_unit_ == "mol") { goto mol; }
//			else if (dens_unit_ == "g") { goto g; }
//			else if (dens_unit_ == "kg") { goto kg; }
//			else if (dens_unit_ == "atom") { goto atom; }
//			else if (dens_unit_ == "atom/(barn-cm)") { goto atombarn; }
//			else
//			{
//				InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: error density unit! only the fllowing units can be used: mol, g, kg, atom, atom/(barn-cm) !!", 1);
//			}
//
//		mol:
//			for (int i = 0; i < num_nucl; ++i)
//			{
//				string nuclstr;
//				getline(modec_inp_, nuclstr);
//				istringstream record(nuclstr);
//				int nucl_id;
//				double n_mol;
//				record >> nucl_id >> n_mol;
//				int id = ModecNuclideLibrary.GetNuclIndex(nucl_id);
//				ModecNuclideLibrary.nuclide_library_vector_[0][id] = n_mol;
//				double awt = nucl_id - nucl_id / 10000 * 10000;
//				if (awt >= 220)
//				{
//					ModecNuclideLibrary.heavy_metal_mass_ += n_mol * awt;
//				}
//			}
//			goto readfile;
//
//		g:
//			for (int i = 0; i < num_nucl; ++i)
//			{
//				string nuclstr;
//				getline(modec_inp_, nuclstr);
//				istringstream record(nuclstr);
//				int nucl_id;
//				double mass;
//				record >> nucl_id >> mass;
//				int id = ModecNuclideLibrary.GetNuclIndex(nucl_id);
//				double awt = double((nucl_id - nucl_id / 10000 * 10000) / 10);
//
//				ModecNuclideLibrary.nuclide_library_vector_[0][id] = mass / awt; // 转换成mol进行存储
//
//				if (awt >= 220)
//				{
//					ModecNuclideLibrary.heavy_metal_mass_ += mass;
//				}
//			}
//			goto readfile;
//
//		kg:
//			for (int i = 0; i < num_nucl; ++i)
//			{
//				string nuclstr;
//				getline(modec_inp_, nuclstr);
//				istringstream record(nuclstr);
//				int nucl_id;
//				double mass;
//				record >> nucl_id >> mass;
//				int id = ModecNuclideLibrary.GetNuclIndex(nucl_id);
//				double awt = double((nucl_id - nucl_id / 10000 * 10000) / 10);
//
//				ModecNuclideLibrary.nuclide_library_vector_[0][id] = mass * 1000.0 / awt; // 转换成mol进行存储
//
//				if (awt >= 220)
//				{
//					ModecNuclideLibrary.heavy_metal_mass_ += mass * 1000.0; // 转换成g来存储
//				}
//			}
//			goto readfile;
//
//		atom:
//			for (int i = 0; i < num_nucl; ++i)
//			{
//				string nuclstr;
//				getline(modec_inp_, nuclstr);
//				istringstream record(nuclstr);
//				int nucl_id;
//				double atom;
//				record >> nucl_id >> atom;
//				int id = ModecNuclideLibrary.GetNuclIndex(nucl_id);
//				double awt = nucl_id - nucl_id / 10000 * 10000;
//
//				ModecNuclideLibrary.nuclide_library_vector_[0][id] = atom / (6.022140857E+23); // 转换成mol进行存储
//
//				if (awt >= 220)
//				{
//					ModecNuclideLibrary.heavy_metal_mass_ += atom / (6.022140857E+23) * awt; // 转换成g来存储
//				}
//			}
//			goto readfile;
//
//		atombarn:
//			for (int i = 0; i < num_nucl; ++i)
//			{
//				string nuclstr;
//				getline(modec_inp_, nuclstr);
//				istringstream record(nuclstr);
//				int nucl_id;
//				double atom;
//				record >> nucl_id >> atom;
//				int id = ModecNuclideLibrary.GetNuclIndex(nucl_id);
//				double awt = nucl_id - nucl_id / 10000 * 10000;
//
//				ModecNuclideLibrary.nuclide_library_vector_[0][id] = atom / (6.022140857E-1); // 转换成mol进行存储
//
//				if (awt >= 220)
//				{
//					ModecNuclideLibrary.heavy_metal_mass_ += atom / (6.022140857E-1) * awt; // 转换成g来存储
//				}
//			}
//			goto readfile;
//		}
//		//////////////////////////////////////////////////
//
//		/////////////////// Block 4 //////////////////////
//		if (tag == "OnlineReprocessing")
//		{
//			record >> if_continously_remove_;
//			if (if_continously_remove_ == true)
//			{
//				string removeline;
//				getline(modec_inp_, removeline);
//				istringstream record(removeline);
//				string removetag;
//				record >> removetag;
//				if (removetag != "Ele_GroupNum")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Ele_GroupNum' tag in the input file.", 1);
//				}
//				record >> remove_group_number_;
//				remove_group_vector_.resize(remove_group_number_);
//				remove_rate_vector_.resize(remove_group_number_);
//				remove_element_vector_.resize(remove_group_number_);
//
//				removeline.clear();
//				record.clear();
//
//				getline(modec_inp_, removeline);
//				record.str(removeline);
//				record >> removetag;
//				if (removetag != "Ele_Group")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Ele_Group' tag in the input file.", 1);
//				}
//				for (int i = 0; i < remove_group_number_; ++i)
//				{
//					record >> remove_group_vector_[i];
//				}
//				removeline.clear();
//				record.clear();
//
//				getline(modec_inp_, removeline);
//				record.str(removeline);
//				record >> removetag;
//				if (removetag != "Ele_RemoveRate")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Ele_RemoveRate' tag in the input file.", 1);
//				}
//				for (int i = 0; i < remove_group_number_; ++i)
//				{
//					record >> remove_rate_vector_[i];
//				}
//				removeline.clear();
//				record.clear();
//
//				for (int i = 0; i < remove_group_number_; ++i)
//				{
//					remove_element_vector_[i].resize(remove_group_vector_[i]);
//				}
//				modec_inp_ >> removetag;
//				if (removetag != "Ele_ID")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Ele_ID' tag in the input file.", 1);
//				}
//				for (int i = 0; i < remove_group_number_; ++i)
//				{
//					int size = remove_group_vector_[i];
//					for (int j = 0; j < size; ++j)
//					{
//						modec_inp_ >> remove_element_vector_[i][j];
//					}
//				}
//				modec_inp_ >> removetag;
//				if (removetag != "TrackingStokage")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'TrackingStokage' tag in the input file.", 1);
//				}
//				modec_inp_ >> if_tracking_stockage;
//				if (if_tracking_stockage == true)
//				{
//					int size = ModecNuclideLibrary.nuclide_library_vector_[0].size();
//					ModecNuclideLibrary.nuclide_library_vector_[0].resize(2 * size); // 若要追踪堆外核素演化，则矩阵中核素总数应为原来的两倍
//				}
//				modec_inp_ >> removetag;
//				if (removetag != "StokageRadioactivity")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'StokageRadioactivity' tag in the input file.", 1);
//				}
//				modec_inp_ >> if_print_stockage_activity_;
//
//				modec_inp_ >> removetag;
//				if (removetag != "StokageDecayEnergy")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'StokageDecayEnergy' tag in the input file.", 1);
//				}
//				modec_inp_ >> if_print_stockage_decayenergy_;
//
//				modec_inp_ >> removetag;
//				if (removetag != "StokageAMPCtoxicity")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'StokageAMPCtoxicity' tag in the input file.", 1);
//				}
//				modec_inp_ >> if_print_stockage_ampc_;
//
//				modec_inp_ >> removetag;
//				if (removetag != "StokageWMPCtoxicity")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'StokageWMPCtoxicity' tag in the input file.", 1);
//				}
//				modec_inp_ >> if_print_stockage_wmpc_;
//
//				modec_inp_ >> removetag;
//				if (removetag != "StokageSvtoxicity")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'StokageSvtoxicity' tag in the input file.", 1);
//				}
//				modec_inp_ >> if_print_stockage_toxicity_;
//			}
//		}
//		//////////////////////////////////////////////////
//
//		if (tag == "KeepingEutecticStable")
//		{
//			record >> if_keeping_eutectic_stable_;
//		}
//
//		/////////////////// Block 4 //////////////////////
//		if (tag == "ContinuouslyFeeding")
//		{
//			record >> if_variable_feeding_;
//			if (if_variable_feeding_ == true)
//			{
//				if (if_constant_online_feeding_ == true)
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: the variables of 'if_variable_feeding_' and 'if_constant_online_feeding_' cannot be both true.", 1);
//				}
//				string feedingline;
//				getline(modec_inp_, feedingline);
//				istringstream record(feedingline);
//				string feedingtag;
//				record >> feedingtag;
//				if (feedingtag != "Nucl_GroupNum")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_GroupNum' tag in the input file.", 1);
//				}
//				record >> variable_feeding_group_num_;
//				variable_feeding_group_vector_.resize(variable_feeding_group_num_);
//				variable_feeding_nuclide_id_vector_.resize(variable_feeding_group_num_);
//				variable_feeding_nuclide_ratio_vector_.resize(variable_feeding_group_num_);
//				feedingline.clear();
//				record.clear();
//
//				getline(modec_inp_, feedingline);
//				record.str(feedingline);
//				record >> feedingtag;
//				if (feedingtag != "Nucl_Ratio")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_Ratio' tag in the input file.", 1);
//				}
//				record >> variable_feeding_ratio_;
//				feedingline.clear();
//				record.clear();
//
//				getline(modec_inp_, feedingline);
//				record.str(feedingline);
//				record >> feedingtag;
//				if (feedingtag != "Nucl_Group")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_Group' tag in the input file.", 1);
//				}
//				for (int i = 0; i < variable_feeding_group_num_; ++i)
//				{
//					record >> variable_feeding_group_vector_[i];
//				}
//				feedingline.clear();
//				record.clear();
//
//				for (int i = 0; i < variable_feeding_group_num_; ++i)
//				{
//					variable_feeding_nuclide_id_vector_[i].resize(variable_feeding_group_vector_[i]);
//					variable_feeding_nuclide_ratio_vector_[i].resize(variable_feeding_group_vector_[i]);
//				}
//				modec_inp_ >> feedingtag;
//				if (feedingtag != "Nucl_ID")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_ID' tag in the input file.", 1);
//				}
//				for (int i = 0; i < variable_feeding_group_num_; ++i)
//				{
//					int size = variable_feeding_group_vector_[i];
//					for (int j = 0; j < size; ++j)
//					{
//						modec_inp_ >> variable_feeding_nuclide_id_vector_[i][j];
//						modec_inp_ >> variable_feeding_nuclide_ratio_vector_[i][j];
//					}
//				}
//			}
//		}
//		//////////////////////////////////////////////////
//		if (tag == "ConstantContinuouslyFeeding")
//		{
//			record >> if_constant_online_feeding_;
//			if (if_constant_online_feeding_ == true)
//			{
//				gauss_legendre_weight_ = GaussLegendreQuadrature::gauss_legendre_weight_[GaussLegendreQuadrature::GL_order];
//				gauss_legendre_abscissa_ = GaussLegendreQuadrature::gauss_legendre_abscissa_[GaussLegendreQuadrature::GL_order];
//				if (if_variable_feeding_ == true)
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: the variables of 'if_variable_feeding_' and 'if_constant_online_feeding_' cannot be both true.", 1);
//				}
//				string feedingline;
//				getline(modec_inp_, feedingline);
//				istringstream record(feedingline);
//				string feedingtag;
//				record >> feedingtag;
//				if (feedingtag != "Nucl_Num")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_Num' tag in the input file.", 1);
//				}
//				record >> constant_feeding_nuclide_num_;
//				constant_feeding_nuclide_id_vector_.resize(constant_feeding_nuclide_num_);
//				constant_feeding_rate_.resize(constant_feeding_nuclide_num_);
//				feedingline.clear();
//				record.clear();
//
//				getline(modec_inp_, feedingline);
//				record.str(feedingline);
//				record >> feedingtag;
//				if (feedingtag != "Nucl_ID")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_ID' tag in the input file.", 1);
//				}
//				for (int i = 0; i < constant_feeding_nuclide_num_; ++i)
//				{
//					record >> constant_feeding_nuclide_id_vector_[i];
//				}
//				feedingline.clear();
//				record.clear();
//
//				getline(modec_inp_, feedingline);
//				record.str(feedingline);
//				record >> feedingtag;
//				if (feedingtag != "Nucl_Rate")
//				{
//					InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no 'Nucl_Rate' tag in the input file.", 1);
//				}
//				for (int i = 0; i < constant_feeding_nuclide_num_; ++i)
//				{
//					record >> constant_feeding_rate_[i];
//				}
//				feedingline.clear();
//				record.clear();
//			}
//		}
//	}
//	modec_inp_.close();
//
//	if (if_flow_mode_ == 1)
//	{
//		if(residue_time_.size() ==0)
//			InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no residue time data are read.", 1);
//
//		if (residue_time_[0] == 0 || residue_time_[1] == 0)
//		{
//			if_flow_mode_ = 0;
//			InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Warning: zero residue time means totally well-mixed.",0);
//			for (int i = 0; i < evolution_value_.size(); ++i)
//			{
//				evolution_value_[i] /= 2.0; // 完全均匀则意味着减半
//			}
//		}
//
//		if(solver_selection_ == 0)
//		{
//			InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Warning: CRAM method must be used when flow depletion is calculated.",0);
//			solver_selection_ == 1;
//		}
//	}
//
//	if (decay_library_name_.length() == 0 && depth_library_name_.length() == 0 && couple_library_name_.length() == 0)
//	{
//		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no library name in the input file.", 1);
//	}
//	if (if_read_mode_tag_ == 0)
//	{
//		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: no depletion model setting in the input file. Depletion model can be set as one or combined of the followings: Decay, Flux, Power.", 1);
//	}
//	if (if_read_density_tag_ == 0)
//	{
//		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: error of getting nuclide densities.", 1);
//	}
//
//	int size_mode = evolution_mode_.size();
//	if (if_calculate_equilibrium_ == true && size_mode != 1)
//	{
//		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: only one depletion model can be set in equilibrium calculation. Depletion model can be set as one of the followings: Decay, Flux, Power.", 1);
//	}
//	if (if_calculate_equilibrium_ == true && if_tracking_stockage == true)
//	{
//		InfoMessage::ErrorMessage("Position: void ModecClass::ModecInitial; \n Error: cannot track nuclide evolution outside the core in equilibrium calculation, the tag 'TrackingStokage' must be set false.", 0);
//		if_tracking_stockage = false;
//	}
//
//	n_vector_.push_back(ModecNuclideLibrary.nuclide_library_vector_[0]);
//
//	if (if_calculate_equilibrium_ == true)
//	{
//		n_vector_.resize(2);
//		n_vector_[1].resize(ModecNuclideLibrary.nuclide_library_vector_[0].size());
//	}
//
//	if (if_constant_online_feeding_ == true && if_tracking_stockage == false && constant_feeding_calculation_methods_ == 1)
//	{
//		constant_feeding_vector_.resize(ModecNuclideLibrary.nuclide_number_);
//		int size = constant_feeding_nuclide_id_vector_.size();
//		for (int i = 0; i < size; ++i)
//		{
//			int Index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
//			constant_feeding_vector_[Index] = constant_feeding_rate_[i];
//		}
//	}
//	if (if_constant_online_feeding_ == true && if_tracking_stockage == true && constant_feeding_calculation_methods_ == 1)
//	{
//		constant_feeding_vector_.resize(ModecNuclideLibrary.nuclide_number_*2);
//		int size = constant_feeding_nuclide_id_vector_.size();
//		for (int i = 0; i < size; ++i)
//		{
//			int Index = ModecNuclideLibrary.GetNuclIndex(constant_feeding_nuclide_id_vector_[i]);
//			constant_feeding_vector_[Index] = constant_feeding_rate_[i];
//		}
//	}
//};

void ModecClass::BuildSpMat() {
    if (lib_tag_ == 0) {
        int size = evolution_mode_.size();
        int i;
        for (i = 0; i < size; ++i) {
            if (evolution_mode_[i] > 0) {
                break;
            }
        }
        if (solver_selection_ == 1) {
            DecayToSpMat();
        } else if(solver_selection_ == 0) {
            DecayToSpMatForTta();
        }
        if (i < size) { // 耦合TRITON
            XSfromTriton();
            CalculateEffectiveFissionYields();
        }
    } else if (lib_tag_ == 1) { // 读取depth_library_name_
        int size = evolution_mode_.size();
        int i;
        for (i = 0; i < size; ++i) {
            if (evolution_mode_[i] > 0) {
                break;
            }
        }
        if (i >= size) { // 纯衰变情形
            if (solver_selection_ == 1) {
                ReadFromDepthLib();
            } else if (solver_selection_ == 0) {
                ReadFromDepthLibForTta();
            }
        } else {
            if (solver_selection_ == 1) {
                ReadFromDepthLib();
                ConstructFissionYieldsSpMat();
            } else if (solver_selection_ == 0) {
                ReadFromDepthLibForTta();
                ConstructFissionYieldsSpMatForTta();
            }
        }
    } else if (lib_tag_ == 2) { // 读取Couple
        int size = evolution_mode_.size();
        int i;
        for (i = 0; i < size; ++i) {
            if (evolution_mode_[i] > 0) {
                break;
            }
        }
        if (solver_selection_ == 1) {
            ReadFromCouple();
        } else if( solver_selection_ == 0 ) {
            ReadFromCoupleForTta();
        }
    }

    if (if_flow_mode_ == 1) {
        TransMatrixPureDecay = TransMatrixDecay; // 定义纯衰变矩阵
    }

    if (if_continously_remove_ == true) {
        AddOnlineReprocessingCoeffi();
    }
    if (if_variable_feeding_ == true) {
        ContinuouslyFeeding();
    }
}
