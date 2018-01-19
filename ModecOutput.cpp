#include "ModecClass.h"

const double Avogadro_Constant = 6.022140857E+23; // 阿伏加德罗常数   Source: 2014 CODATA  
const double Electron_Coulomb = 1.6021766208E-19;//1.6021766208E-19; // 电子电量        Source: 2014 CODATA  

const double cutoff = 1e-50; // 浓度的截断误差

void ModecClass::ModecOutput()
{
	int Nucl_size = ModecNuclideLibrary.nuclide_number_;

	if (if_calculate_equilibrium_ == false)
	{
		if (input_filename_ != "modec.inp")
		{
			output_filename_ = input_filename_ + ".concentration";
		}

		output_file_ = work_direc_ + output_filename_;
		modec_out_.open(output_file_);
		modec_out_ << "*********************************** Nuclide Concentrations, " << dens_unit_ << " ***********************************";
		modec_out_ << '\n';
		modec_out_ << '\n';
		modec_out_.setf(ios::left);
		int size = evolution_mode_.size();
		int size_tot = 0;
		vector<string> time_unit;
		for (int i = 0; i < size; ++i)
		{
			size_tot += substep_[i];		
		}

		vector<double> Cum_time;
		Cum_time.resize(size_tot + 1);

		int count = 0;
		for (int i = 0; i < size; ++i)
		{
			for (int j = 0; j < substep_[i]; ++j)
			{
				Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
				time_unit.push_back(time_unit_[i]);
				count++;
			}
		}


		modec_out_.width(8);
		modec_out_ << "";
		modec_out_.width(8);
		modec_out_ << "Time: ";
		modec_out_.width(18);
		modec_out_ << "ModecInitial";

		if (print_mode_ == 0)
		{
			modec_out_.width(10);
			modec_out_.setf(ios::right);
			if (time_unit[size_tot - 1] == "d" || time_unit[size_tot - 1] == "day" || time_unit[size_tot - 1] == "days")
			{
				Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0;
			}
			if (time_unit[size_tot - 1] == "m" || time_unit[size_tot - 1] == "month" || time_unit[size_tot - 1] == "months")
			{
				Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0 / 30.0;
			}
			else if (time_unit[size_tot - 1] == "y" || time_unit[size_tot - 1] == "year" || time_unit[size_tot - 1] == "years")
			{
				Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0 / 365.25;
				time_unit[size_tot - 1] = "yr";
			}

			if (Cum_time[size_tot] >= 1.0e+3)
			{
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
			}
			else
			{
				modec_out_.precision(5);
			}
			
			modec_out_ << Cum_time[size_tot];
			modec_out_.unsetf(ios::right);
			modec_out_ << " ";
			modec_out_.width(7);
			modec_out_ << time_unit[size_tot - 1];
		}
		else
		{
			for (int j = 1; j <= size_tot; ++j)
			{
				modec_out_.setf(ios::right);
				if (time_unit[j - 1] == "d" || time_unit[j - 1] == "day" || time_unit[j - 1] == "days")
				{
					Cum_time[j] = Cum_time[j] / 24.0 / 3600.0;
				}
				if (time_unit[j - 1] == "m" || time_unit[j - 1] == "month" || time_unit[j - 1] == "months")
				{
					Cum_time[j] = Cum_time[j] / 24.0 / 3600.0 / 30.0;
				}
				else if (time_unit[j - 1] == "y" || time_unit[j - 1] == "year" || time_unit[j - 1] == "years")
				{
					Cum_time[j] = Cum_time[j] / 24.0 / 3600.0 / 365.25;
					time_unit[j - 1] = "yr";
				}
				modec_out_.width(10);
				if (Cum_time[j] >= 1.0e+3)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
				}
				else
				{
					modec_out_.precision(5);
				}
				modec_out_ << Cum_time[j];
				modec_out_.unsetf(ios::right);
				modec_out_ << " ";
				modec_out_.width(7);
				modec_out_ << time_unit[j - 1];
			}
		}

		modec_out_ << '\n';

		modec_out_.width(8);
		modec_out_ << "";
		modec_out_.width(8);
		modec_out_ << "Power: ";

		if (print_mode_ == 0)
		{
			modec_out_.width(10);
			modec_out_.setf(ios::scientific | ios::uppercase);
			modec_out_.precision(4);
			modec_out_ << power_vector_[0];
			modec_out_.width(8);
			modec_out_ << "MW";

			modec_out_.width(10);
			modec_out_.setf(ios::scientific | ios::uppercase);
			modec_out_.precision(4);
			modec_out_ << power_vector_[size_tot];
			modec_out_.width(8);
			modec_out_ << "MW";
		}
		else
		{
			/*modec_out_.width(10);
			modec_out_.setf(ios::scientific | ios::uppercase);
			modec_out_.precision(4);
			modec_out_ << 0.0;
			modec_out_.width(8);
			modec_out_ << "MW";*/

			for (int j = 0; j <= size_tot; ++j)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << power_vector_[j];
				modec_out_.width(8);
				modec_out_ << " MW";
			}
		}

		modec_out_ << '\n';

		modec_out_.width(8);
		modec_out_ << "";
		modec_out_.width(8);
		modec_out_ << "Flux: ";

		if (print_mode_ == 0)
		{
			modec_out_.width(10);
			modec_out_.setf(ios::scientific | ios::uppercase);
			modec_out_.precision(4);
			modec_out_ << flux_vector_[0];
			modec_out_.width(8);
			modec_out_ << "";

			modec_out_.width(10);
			modec_out_.setf(ios::scientific | ios::uppercase);
			modec_out_.precision(4);
			modec_out_ << flux_vector_[size_tot];
			modec_out_.width(8);
			modec_out_ << "";
		}
		else
		{
			//modec_out_.width(10);
			//modec_out_.setf(ios::scientific | ios::uppercase);
			//modec_out_.precision(4);
			//modec_out_ << 0.0;
			//modec_out_.width(8);
			//modec_out_ << "";
			for (int j = 0; j <= size_tot; ++j)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << flux_vector_[j];
				modec_out_.width(8);
				modec_out_ << "";
			}
		}
		modec_out_ << '\n';
		
//////////////    Calculate kinf_vector_, Neutron_Production_Rate and Neutron_Absorption_Rate      /////////
		vector<double> kinf;
		vector<double> prod_neu;
		vector<double> absorption_neu;
		if (if_print_kinf_ != 0 || if_print_fission_rate_ != 0 || if_print_absorption_rate_ != 0)
		{		
			kinf.resize(size_tot + 1);
			prod_neu.resize(size_tot + 1);
			absorption_neu.resize(size_tot + 1);

			for (int j = 0; j <= size_tot; ++j)
			{
				for (int i = 0; i < Nucl_size; ++i)
				{
					prod_neu[j] += n_vector_[j][i] * Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24;
					absorption_neu[j] += n_vector_[j][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24;
				}
				prod_neu[j] = prod_neu[j] * flux_vector_[j];
				absorption_neu[j] = absorption_neu[j] * flux_vector_[j];

				if (absorption_neu[j] != 0.0)
				{
					kinf[j] = prod_neu[j] / absorption_neu[j];
				}
				
			}
			
			if (if_print_kinf_ == true)
			{

				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(8);
				modec_out_ << "kinf: ";

				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << kinf[0];
					modec_out_.width(18);
					modec_out_ << kinf[size_tot];
					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						modec_out_ << kinf[j];
					}
					modec_out_ << '\n';
				}
			}
			
			if (if_print_fission_rate_ != 0)
			{

				modec_out_.width(2);
				modec_out_ << "";
				modec_out_.width(14);
				modec_out_ << "NeuProdRate: ";

				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << prod_neu[0];
					modec_out_.width(18);
					modec_out_ << prod_neu[size_tot];
					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						modec_out_ << prod_neu[j];
					}
					modec_out_ << '\n';
				}
			}

			if (if_print_absorption_rate_ != 0)
			{

				modec_out_.width(2);
				modec_out_ << "";
				modec_out_.width(14);
				modec_out_ << "NeuAbsRate: ";

				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << absorption_neu[0];
					modec_out_.width(18);
					modec_out_ << absorption_neu[size_tot];
					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						modec_out_ << absorption_neu[j];
					}
					modec_out_ << '\n';
				}
			}
			modec_out_ << '\n';
		}
//////////////////////////////////////////////////////////////////////////////////////////////////


		modec_out_ << '\n';

		for (unsigned int i = 0; i < Nucl_size; ++i)
		{
			modec_out_.width(8);
			int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
			modec_out_ << nucl_id;
			int NZ = nucl_id / 10000;
			int NA = (nucl_id - NZ * 10000) / 10;
			stringstream ss;
			ss << NA;
			int NG = nucl_id - nucl_id / 10 * 10;
			string Nucl_Name;
			Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
			if (NG == 1) { Nucl_Name += "m"; }
			modec_out_.width(8);
			modec_out_ << Nucl_Name;

			double convert_coeff;

			if (dens_unit_ == "mol") { convert_coeff = 1.0; }
			else if (dens_unit_ == "g") { convert_coeff = double(NA); }
			else if (dens_unit_ == "kg") { convert_coeff = double(NA) / 1000.0; }
			else if (dens_unit_ == "atom") { convert_coeff = Avogadro_Constant; }
			else if (dens_unit_ == "atom/(barn-cm)") { convert_coeff = Avogadro_Constant * (1.0e-24); }

			if (print_mode_ == 0)
			{
				modec_out_.width(18);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(9);
				modec_out_ << n_vector_[0][i] * convert_coeff;
				modec_out_.width(18);
				if (abs(n_vector_[size_tot][i]) < cutoff)
				{
					modec_out_ << 0.0;
				}
				else
				{
					modec_out_ << n_vector_[size_tot][i] * convert_coeff;
				}

				modec_out_ << '\n';
			}
			if (print_mode_ == 1)
			{
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(9);
				for (int j = 0; j <= size_tot; ++j)
				{
					modec_out_.width(18);
					if (abs(n_vector_[j][i]) < cutoff || n_vector_[j][i] < 0.0)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[j][i] * convert_coeff;
					}
				}
				modec_out_ << '\n';
			}
		}

		modec_out_.close();

		if (if_print_activity_ == true)
		{

			vector<double > lamda;
			lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

			output_filename_ = input_filename_ + ".radioactivity";

			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Radioactivity, Ci ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';

			modec_out_.setf(ios::left);
			int size = evolution_mode_.size();
			int size_tot = 0;
			for (int i = 0; i < size; ++i)
			{
				size_tot += substep_[i];
			}

			vector<double> Cum_time;
			Cum_time.resize(size_tot + 1);

			int count = 0;
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < substep_[i]; ++j)
				{
					Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
					count++;
				}
			}


			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			if (print_mode_ == 0)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << Cum_time[size_tot];
				modec_out_.width(8);
				modec_out_ << "s";
			}
			else
			{
				for (int j = 1; j <= size_tot; ++j)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[j];
					modec_out_.width(8);
					modec_out_ << "s";
				}
			}

			modec_out_ << '\n';
			modec_out_ << '\n';

			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;
				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / (3.7e10);
					modec_out_.width(18);
					if (abs(n_vector_[size_tot][i]) < cutoff || n_vector_[size_tot][i] < 0.0)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
					}

					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						if (abs(n_vector_[j][i]) < cutoff || n_vector_[j][i] < 0.0)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[j][i]*Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
						}
					}
					modec_out_ << '\n';
				}
			}

			modec_out_.close();

		}

		if (if_print_decayenergy_ == true)
		{
			vector<double > lamda;
			lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

			vector<double > Q;
			Q = ModecNuclideLibrary.nuclide_library_vector_[3];

			output_filename_ = input_filename_ + ".decayheat";

			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Decay Heat, Watts ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';

			modec_out_.setf(ios::left);
			int size = evolution_mode_.size();
			int size_tot = 0;
			for (int i = 0; i < size; ++i)
			{
				size_tot += substep_[i];
			}

			vector<double> Cum_time;
			Cum_time.resize(size_tot + 1);

			int count = 0;
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < substep_[i]; ++j)
				{
					Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
					count++;
				}
			}


			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			if (print_mode_ == 0)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << Cum_time[size_tot];
				modec_out_.width(8);
				modec_out_ << "s";
			}
			else
			{
				for (int j = 1; j <= size_tot; ++j)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[j];
					modec_out_.width(8);
					modec_out_ << "s";
				}
			}

			modec_out_ << '\n';
			modec_out_ << '\n';

			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;
				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6;
					modec_out_.width(18);
					if (abs(n_vector_[size_tot][i]) < cutoff || n_vector_[size_tot][i] < 0.0)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6; // 转换成单位 W
					}

					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						if (abs(n_vector_[j][i]) < cutoff || n_vector_[j][i] < 0.0)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[j][i]*Avogadro_Constant*lamda[i]*Q[i]* Electron_Coulomb*1.0e6;
						}
					}
					modec_out_ << '\n';
				}
			}

			modec_out_.close();

		}

		if (if_print_ampc_ == true)
		{
			vector<double > lamda;
			lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

			vector<double > AMPC;
			AMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

			output_filename_ = input_filename_ + ".ampc";

			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Radioactivity, m3-air/Bq ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';

			modec_out_.setf(ios::left);
			int size = evolution_mode_.size();
			int size_tot = 0;
			for (int i = 0; i < size; ++i)
			{
				size_tot += substep_[i];
			}

			vector<double> Cum_time;
			Cum_time.resize(size_tot + 1);

			int count = 0;
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < substep_[i]; ++j)
				{
					Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
					count++;
				}
			}


			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			if (print_mode_ == 0)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << Cum_time[size_tot];
				modec_out_.width(8);
				modec_out_ << "s";
			}
			else
			{
				for (int j = 1; j <= size_tot; ++j)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[j];
					modec_out_.width(8);
					modec_out_ << "s";
				}
			}

			modec_out_ << '\n';
			modec_out_ << '\n';

			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;
				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
					modec_out_.width(18);
					if (abs(n_vector_[size_tot][i]) < cutoff || n_vector_[size_tot][i] < 0.0)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
					}

					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						if (abs(n_vector_[j][i]) < cutoff || n_vector_[j][i] < 0.0)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[j][i]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
						}
					}
					modec_out_ << '\n';
				}
			}

			modec_out_.close();

		}

		if (if_print_wmpc_ == true)
		{
			vector<double > lamda;
			lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

			vector<double > WMPC;
			WMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

			output_filename_ = input_filename_ + ".wmpc";

			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Radioactivity, m3-air/Bq ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';

			modec_out_.setf(ios::left);
			int size = evolution_mode_.size();
			int size_tot = 0;
			for (int i = 0; i < size; ++i)
			{
				size_tot += substep_[i];
			}

			vector<double> Cum_time;
			Cum_time.resize(size_tot + 1);

			int count = 0;
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < substep_[i]; ++j)
				{
					Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
					count++;
				}
			}


			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			if (print_mode_ == 0)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << Cum_time[size_tot];
				modec_out_.width(8);
				modec_out_ << "s";
			}
			else
			{
				for (int j = 1; j <= size_tot; ++j)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[j];
					modec_out_.width(8);
					modec_out_ << "s";
				}
			}

			modec_out_ << '\n';
			modec_out_ << '\n';

			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;
				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
					modec_out_.width(18);
					if (abs(n_vector_[size_tot][i]) < cutoff || n_vector_[size_tot][i] < 0.0)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
					}

					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						if (abs(n_vector_[j][i]) < cutoff || n_vector_[j][i] < 0.0)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[j][i]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
						}
					}
					modec_out_ << '\n';
				}
			}

			modec_out_.close();

		}

		if (if_print_toxicity_ == true)
		{
			vector<double > lamda;
			lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

			vector<double > TOXCI;
			TOXCI = ModecNuclideLibrary.nuclide_library_vector_[10];

			output_filename_ = input_filename_ + ".toxicity";

			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Radiotoxicity, Sv ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';

			modec_out_.setf(ios::left);
			int size = evolution_mode_.size();
			int size_tot = 0;
			for (int i = 0; i < size; ++i)
			{
				size_tot += substep_[i];
			}

			vector<double> Cum_time;
			Cum_time.resize(size_tot + 1);

			int count = 0;
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < substep_[i]; ++j)
				{
					Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
					count++;
				}
			}


			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			if (print_mode_ == 0)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << Cum_time[size_tot];
				modec_out_.width(8);
				modec_out_ << "s";
			}
			else
			{
				for (int j = 1; j <= size_tot; ++j)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[j];
					modec_out_.width(8);
					modec_out_ << "s";
				}
			}

			modec_out_ << '\n';
			modec_out_ << '\n';

			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;
				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] * TOXCI[i];
					modec_out_.width(18);
					if (abs(n_vector_[size_tot][i]) < cutoff || n_vector_[size_tot][i] < 0.0)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant*lamda[i] * TOXCI[i];
					}

					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						if (abs(n_vector_[j][i]) < cutoff || n_vector_[j][i] < 0.0)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[j][i]*Avogadro_Constant*lamda[i] * TOXCI[i];
						}
					}
					modec_out_ << '\n';
				}
			}

			modec_out_.close();

		}

		if (if_print_fission_rate_ == 2)
		{
			output_filename_ = input_filename_ + ".NeuProdRate";

			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Neutron Production Rate, s-1***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';

			modec_out_.setf(ios::left);
			int size = evolution_mode_.size();
			int size_tot = 0;
			for (int i = 0; i < size; ++i)
			{
				size_tot += substep_[i];
			}

			vector<double> Cum_time;
			Cum_time.resize(size_tot + 1);

			int count = 0;
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < substep_[i]; ++j)
				{
					Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
					count++;
				}
			}


			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			if (print_mode_ == 0)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << Cum_time[size_tot];
				modec_out_.width(8);
				modec_out_ << "s";
			}
			else
			{
				for (int j = 1; j <= size_tot; ++j)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[j];
					modec_out_.width(8);
					modec_out_ << "s";
				}
			}

			modec_out_ << '\n';
			modec_out_ << '\n';

			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;
				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24 * flux_vector_[0];
					modec_out_.width(18);

					if (abs(n_vector_[size_tot][i]) < cutoff || n_vector_[size_tot][i] < 0.0)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24* flux_vector_[size_tot];
					}
					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						if (n_vector_[j][i] < cutoff || n_vector_[j][i] < 0.0)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[j][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[9][i] * 1.0e-24* flux_vector_[j]; // 单位为Ci
						}
					}
					modec_out_ << '\n';
				}
			}
			modec_out_ << '\n';
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "Total: ";
			for (int j = 0; j <= size_tot; ++j)
			{
				modec_out_.width(18);

				modec_out_ << prod_neu[j]; // 单位为Ci

			}
			modec_out_ << '\n';

			modec_out_.close();

		}

		if (if_print_absorption_rate_ == 2)
		{
			output_filename_ = input_filename_ + ".NeuAbsRate";

			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Neutron Absorption Rate, s-1 ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';

			modec_out_.setf(ios::left);
			int size = evolution_mode_.size();
			int size_tot = 0;
			for (int i = 0; i < size; ++i)
			{
				size_tot += substep_[i];
			}

			vector<double> Cum_time;
			Cum_time.resize(size_tot + 1);

			int count = 0;
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < substep_[i]; ++j)
				{
					Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
					count++;
				}
			}


			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			if (print_mode_ == 0)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << Cum_time[size_tot];
				modec_out_.width(8);
				modec_out_ << "s";
			}
			else
			{
				for (int j = 1; j <= size_tot; ++j)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[j];
					modec_out_.width(8);
					modec_out_ << "s";
				}
			}

			modec_out_ << '\n';
			modec_out_ << '\n';

			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;
				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24* flux_vector_[0];
					modec_out_.width(18);
					if (n_vector_[size_tot][i] < cutoff)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[size_tot][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24* flux_vector_[size_tot];
					}
					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						if (n_vector_[j][i] < cutoff)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[j][i]*Avogadro_Constant * ModecNuclideLibrary.nuclide_library_vector_[6][i] * 1.0e-24* flux_vector_[j];
						}
					}
					modec_out_ << '\n';
				}
			}

			modec_out_ << '\n';
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "Total: ";
			for (int j = 0; j <= size_tot; ++j)
			{
				modec_out_.width(18);

				modec_out_ << absorption_neu[j]; // 单位为Ci

			}
			modec_out_ << '\n';

			modec_out_.close();

		}
		if (if_tracking_stockage == true)
		{
			output_filename_ = input_filename_ + ".stockage.concentration";

			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Concentrations, mol ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';
			modec_out_.setf(ios::left);
			int size = evolution_mode_.size();
			int size_tot = 0;
			vector<string> time_unit;
			for (int i = 0; i < size; ++i)
			{
				size_tot += substep_[i];
			}

			vector<double> Cum_time;
			Cum_time.resize(size_tot + 1);

			int count = 0;
			for (int i = 0; i < size; ++i)
			{
				for (int j = 0; j < substep_[i]; ++j)
				{
					Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
					time_unit.push_back(time_unit_[i]);
					count++;
				}
			}


			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "Time: ";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			if (print_mode_ == 0)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::right);
				if (time_unit[size_tot - 1] == "d" || time_unit[size_tot - 1] == "day" || time_unit[size_tot - 1] == "days")
				{
					Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0;
				}
				if (time_unit[size_tot - 1] == "m" || time_unit[size_tot - 1] == "month" || time_unit[size_tot - 1] == "months")
				{
					Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0 / 30.0;
				}
				else if (time_unit[size_tot - 1] == "y" || time_unit[size_tot - 1] == "year" || time_unit[size_tot - 1] == "years")
				{
					Cum_time[size_tot] = Cum_time[size_tot] / 24.0 / 3600.0 / 365.25;
					time_unit[size_tot - 1] = "yr";
				}

				if (Cum_time[size_tot] >= 1.0e+3)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
				}
				else
				{
					modec_out_.unsetf(ios::scientific | ios::uppercase);
					modec_out_.precision(5);
				}
				
				modec_out_ << Cum_time[size_tot];
				modec_out_.unsetf(ios::right);
				modec_out_ << " ";
				modec_out_.width(7);
				modec_out_ << time_unit[size_tot - 1];
			}
			else
			{
				for (int j = 1; j <= size_tot; ++j)
				{
					modec_out_.setf(ios::right);
					if (time_unit[j - 1] == "d" || time_unit[j - 1] == "day" || time_unit[j - 1] == "days")
					{
						Cum_time[j] = Cum_time[j] / 24.0 / 3600.0;
					}
					if (time_unit[j - 1] == "m" || time_unit[j - 1] == "month" || time_unit[j - 1] == "months")
					{
						Cum_time[j] = Cum_time[j] / 24.0 / 3600.0 / 30.0;
					}
					else if (time_unit[j - 1] == "y" || time_unit[j - 1] == "year" || time_unit[j - 1] == "years")
					{
						Cum_time[j] = Cum_time[j] / 24.0 / 3600.0 / 365.25;
						time_unit[j - 1] = "yr";
					}
					modec_out_.width(10);
					if (Cum_time[j] >= 1.0e+3)
					{
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(4);
					}
					else
					{
						modec_out_.unsetf(ios::scientific | ios::uppercase);
						modec_out_.precision(5);
					}
					
					modec_out_ << Cum_time[j];
					modec_out_.unsetf(ios::right);
					modec_out_ << " ";
					modec_out_.width(7);
					modec_out_ << time_unit[j - 1];
				}
			}

			modec_out_ << '\n';

			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "Power: ";

			if (print_mode_ == 0)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << 0.0;
				modec_out_.width(8);
				modec_out_ << "MW";

				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << 0.0;
				modec_out_.width(8);
				modec_out_ << "MW";
			}
			else
			{
				/*modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << 0.0;
				modec_out_.width(8);
				modec_out_ << "MW";*/

				for (int j = 0; j <= size_tot; ++j)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << 0.0;
					modec_out_.width(8);
					modec_out_ << " MW";
				}
			}

			modec_out_ << '\n';

			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "flux_: ";

			if (print_mode_ == 0)
			{
				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << 0.0;
				modec_out_.width(8);
				modec_out_ << "";

				modec_out_.width(10);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(4);
				modec_out_ << 0.0;
				modec_out_.width(8);
				modec_out_ << "";
			}
			else
			{
				//modec_out_.width(10);
				//modec_out_.setf(ios::scientific | ios::uppercase);
				//modec_out_.precision(4);
				//modec_out_ << 0.0;
				//modec_out_.width(8);
				//modec_out_ << "";
				for (int j = 0; j <= size_tot; ++j)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << 0.0;
					modec_out_.width(8);
					modec_out_ << "";
				}
			}
			modec_out_ << '\n';
			modec_out_ << '\n';

			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;
				if (print_mode_ == 0)
				{
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i + Nucl_size];
					modec_out_.width(18);
					if (n_vector_[size_tot][i + Nucl_size] < cutoff)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[size_tot][i + Nucl_size];
					}

					modec_out_ << '\n';
				}
				if (print_mode_ == 1)
				{
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					for (int j = 0; j <= size_tot; ++j)
					{
						modec_out_.width(18);
						if (n_vector_[j][i + Nucl_size] < cutoff)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[j][i + Nucl_size];
						}
					}
					modec_out_ << '\n';
				}
			}

			modec_out_.close();

			if (if_print_stockage_activity_ == true)
			{

				vector<double > lamda;
				lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

				output_filename_ = input_filename_ + ".stockage.radioactivity";

				output_file_ = work_direc_ + output_filename_;
				modec_out_.open(output_file_);
				modec_out_ << "*********************************** Nuclide Radioactivity, Ci ***********************************";
				modec_out_ << '\n';
				modec_out_ << '\n';

				modec_out_.setf(ios::left);
				int size = evolution_mode_.size();
				int size_tot = 0;
				for (int i = 0; i < size; ++i)
				{
					size_tot += substep_[i];
				}

				vector<double> Cum_time;
				Cum_time.resize(size_tot + 1);

				int count = 0;
				for (int i = 0; i < size; ++i)
				{
					for (int j = 0; j < substep_[i]; ++j)
					{
						Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
						count++;
					}
				}


				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(18);
				modec_out_ << "ModecInitial";

				if (print_mode_ == 0)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[size_tot];
					modec_out_.width(8);
					modec_out_ << "s";
				}
				else
				{
					for (int j = 1; j <= size_tot; ++j)
					{
						modec_out_.width(10);
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(4);
						modec_out_ << Cum_time[j];
						modec_out_.width(8);
						modec_out_ << "s";
					}
				}

				modec_out_ << '\n';
				modec_out_ << '\n';

				for (unsigned int i = 0; i < Nucl_size; ++i)
				{
					modec_out_.width(8);
					int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
					modec_out_ << nucl_id;
					int NZ = nucl_id / 10000;
					int NA = (nucl_id - NZ * 10000) / 10;
					stringstream ss;
					ss << NA;
					int NG = nucl_id - nucl_id / 10 * 10;
					string Nucl_Name;
					Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
					if (NG == 1) { Nucl_Name += "m"; }
					modec_out_.width(8);
					modec_out_ << Nucl_Name;
					if (print_mode_ == 0)
					{
						modec_out_.width(18);
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(9);
						modec_out_ << n_vector_[0][i + Nucl_size]*Avogadro_Constant*lamda[i] / (3.7e10);
						modec_out_.width(18);
						if (n_vector_[size_tot][i + Nucl_size] < cutoff)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[size_tot][i + Nucl_size]*Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
						}

						modec_out_ << '\n';
					}
					if (print_mode_ == 1)
					{
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(9);
						for (int j = 0; j <= size_tot; ++j)
						{
							modec_out_.width(18);
							if (n_vector_[j][i + Nucl_size] < cutoff)
							{
								modec_out_ << 0.0;
							}
							else
							{
								modec_out_ << n_vector_[j][i + Nucl_size]*Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
							}
						}
						modec_out_ << '\n';
					}
				}

				modec_out_.close();

			}

			if (if_print_stockage_decayenergy_ == true)
			{
				vector<double > lamda;
				lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

				vector<double > Q;
				Q = ModecNuclideLibrary.nuclide_library_vector_[3];

				output_filename_ = input_filename_ + ".stockage.decayheat";
				output_file_ = work_direc_ + output_filename_;
				modec_out_.open(output_file_);
				modec_out_ << "*********************************** Nuclide Decay Heat, Watts ***********************************";
				modec_out_ << '\n';
				modec_out_ << '\n';

				modec_out_.setf(ios::left);
				int size = evolution_mode_.size();
				int size_tot = 0;
				for (int i = 0; i < size; ++i)
				{
					size_tot += substep_[i];
				}

				vector<double> Cum_time;
				Cum_time.resize(size_tot + 1);

				int count = 0;
				for (int i = 0; i < size; ++i)
				{
					for (int j = 0; j < substep_[i]; ++j)
					{
						Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
						count++;
					}
				}


				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(18);
				modec_out_ << "ModecInitial";

				if (print_mode_ == 0)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[size_tot];
					modec_out_.width(8);
					modec_out_ << "s";
				}
				else
				{
					for (int j = 1; j <= size_tot; ++j)
					{
						modec_out_.width(10);
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(4);
						modec_out_ << Cum_time[j];
						modec_out_.width(8);
						modec_out_ << "s";
					}
				}

				modec_out_ << '\n';
				modec_out_ << '\n';

				for (unsigned int i = 0; i < Nucl_size; ++i)
				{
					modec_out_.width(8);
					int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
					modec_out_ << nucl_id;
					int NZ = nucl_id / 10000;
					int NA = (nucl_id - NZ * 10000) / 10;
					stringstream ss;
					ss << NA;
					int NG = nucl_id - nucl_id / 10 * 10;
					string Nucl_Name;
					Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
					if (NG == 1) { Nucl_Name += "m"; }
					modec_out_.width(8);
					modec_out_ << Nucl_Name;
					if (print_mode_ == 0)
					{
						modec_out_.width(18);
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(9);
						modec_out_ << n_vector_[0][i + Nucl_size]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6;
						modec_out_.width(18);
						if (n_vector_[size_tot][i + Nucl_size] < cutoff)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[size_tot][i + Nucl_size]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6; // 转换成单位 W
						}

						modec_out_ << '\n';
					}
					if (print_mode_ == 1)
					{
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(9);
						for (int j = 0; j <= size_tot; ++j)
						{
							modec_out_.width(18);
							if (n_vector_[j][i + Nucl_size] < cutoff)
							{
								modec_out_ << 0.0;
							}
							else
							{
								modec_out_ << n_vector_[j][i + Nucl_size]*Avogadro_Constant*lamda[i]*Q[i]* Electron_Coulomb*1.0e6;
							}
						}
						modec_out_ << '\n';
					}
				}

				modec_out_.close();

			}

			if (if_print_stockage_ampc_ == true)
			{
				vector<double > lamda;
				lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

				vector<double > AMPC;
				AMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

				output_filename_ = input_filename_ + ".stockage.ampc";
				output_file_ = work_direc_ + output_filename_;
				modec_out_.open(output_file_);
				modec_out_ << "*********************************** Nuclide Radioactivity, m3-air/Bq ***********************************";
				modec_out_ << '\n';
				modec_out_ << '\n';

				modec_out_.setf(ios::left);
				int size = evolution_mode_.size();
				int size_tot = 0;
				for (int i = 0; i < size; ++i)
				{
					size_tot += substep_[i];
				}

				vector<double> Cum_time;
				Cum_time.resize(size_tot + 1);

				int count = 0;
				for (int i = 0; i < size; ++i)
				{
					for (int j = 0; j < substep_[i]; ++j)
					{
						Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
						count++;
					}
				}


				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(18);
				modec_out_ << "ModecInitial";

				if (print_mode_ == 0)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[size_tot];
					modec_out_.width(8);
					modec_out_ << "s";
				}
				else
				{
					for (int j = 1; j <= size_tot; ++j)
					{
						modec_out_.width(10);
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(4);
						modec_out_ << Cum_time[j];
						modec_out_.width(8);
						modec_out_ << "s";
					}
				}

				modec_out_ << '\n';
				modec_out_ << '\n';

				for (unsigned int i = 0; i < Nucl_size; ++i)
				{
					modec_out_.width(8);
					int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
					modec_out_ << nucl_id;
					int NZ = nucl_id / 10000;
					int NA = (nucl_id - NZ * 10000) / 10;
					stringstream ss;
					ss << NA;
					int NG = nucl_id - nucl_id / 10 * 10;
					string Nucl_Name;
					Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
					if (NG == 1) { Nucl_Name += "m"; }
					modec_out_.width(8);
					modec_out_ << Nucl_Name;
					if (print_mode_ == 0)
					{
						modec_out_.width(18);
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(9);
						modec_out_ << n_vector_[0][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
						modec_out_.width(18);
						if (n_vector_[size_tot][i + Nucl_size] < cutoff)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[size_tot][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
						}

						modec_out_ << '\n';
					}
					if (print_mode_ == 1)
					{
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(9);
						for (int j = 0; j <= size_tot; ++j)
						{
							modec_out_.width(18);
							if (n_vector_[j][i + Nucl_size] < cutoff)
							{
								modec_out_ << 0.0;
							}
							else
							{
								modec_out_ << n_vector_[j][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
							}
						}
						modec_out_ << '\n';
					}
				}

				modec_out_.close();

			}

			if (if_print_stockage_wmpc_ == true)
			{
				vector<double > lamda;
				lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

				vector<double > WMPC;
				WMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

				output_filename_ = input_filename_ + ".stockage.wmpc";
				output_file_ = work_direc_ + output_filename_;
				modec_out_.open(output_file_);
				modec_out_ << "*********************************** Nuclide Radioactivity, m3-air/Bq ***********************************";
				modec_out_ << '\n';
				modec_out_ << '\n';

				modec_out_.setf(ios::left);
				int size = evolution_mode_.size();
				int size_tot = 0;
				for (int i = 0; i < size; ++i)
				{
					size_tot += substep_[i];
				}

				vector<double> Cum_time;
				Cum_time.resize(size_tot + 1);

				int count = 0;
				for (int i = 0; i < size; ++i)
				{
					for (int j = 0; j < substep_[i]; ++j)
					{
						Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
						count++;
					}
				}


				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(18);
				modec_out_ << "ModecInitial";

				if (print_mode_ == 0)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[size_tot];
					modec_out_.width(8);
					modec_out_ << "s";
				}
				else
				{
					for (int j = 1; j <= size_tot; ++j)
					{
						modec_out_.width(10);
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(4);
						modec_out_ << Cum_time[j];
						modec_out_.width(8);
						modec_out_ << "s";
					}
				}

				modec_out_ << '\n';
				modec_out_ << '\n';

				for (unsigned int i = 0; i < Nucl_size; ++i)
				{
					modec_out_.width(8);
					int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
					modec_out_ << nucl_id;
					int NZ = nucl_id / 10000;
					int NA = (nucl_id - NZ * 10000) / 10;
					stringstream ss;
					ss << NA;
					int NG = nucl_id - nucl_id / 10 * 10;
					string Nucl_Name;
					Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
					if (NG == 1) { Nucl_Name += "m"; }
					modec_out_.width(8);
					modec_out_ << Nucl_Name;
					if (print_mode_ == 0)
					{
						modec_out_.width(18);
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(9);
						modec_out_ << n_vector_[0][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
						modec_out_.width(18);
						if (n_vector_[size_tot][i + Nucl_size] < cutoff)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[size_tot][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
						}

						modec_out_ << '\n';
					}
					if (print_mode_ == 1)
					{
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(9);
						for (int j = 0; j <= size_tot; ++j)
						{
							modec_out_.width(18);
							if (n_vector_[j][i + Nucl_size] < cutoff)
							{
								modec_out_ << 0.0;
							}
							else
							{
								modec_out_ << n_vector_[j][i + Nucl_size]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
							}
						}
						modec_out_ << '\n';
					}
				}

				modec_out_.close();

			}

			if (if_print_stockage_toxicity_ == true)
			{
				vector<double > lamda;
				lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

				vector<double > TOXCI;
				TOXCI = ModecNuclideLibrary.nuclide_library_vector_[10];

				output_filename_ = input_filename_ + ".stockage.toxicity";
				output_file_ = work_direc_ + output_filename_;
				modec_out_.open(output_file_);
				modec_out_ << "*********************************** Nuclide Radiotoxicity, Sv ***********************************";
				modec_out_ << '\n';
				modec_out_ << '\n';

				modec_out_.setf(ios::left);
				int size = evolution_mode_.size();
				int size_tot = 0;
				for (int i = 0; i < size; ++i)
				{
					size_tot += substep_[i];
				}

				vector<double> Cum_time;
				Cum_time.resize(size_tot + 1);

				int count = 0;
				for (int i = 0; i < size; ++i)
				{
					for (int j = 0; j < substep_[i]; ++j)
					{
						Cum_time[count + 1] = Cum_time[count] + burnup_time_[i];
						count++;
					}
				}


				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(8);
				modec_out_ << "";
				modec_out_.width(18);
				modec_out_ << "ModecInitial";

				if (print_mode_ == 0)
				{
					modec_out_.width(10);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(4);
					modec_out_ << Cum_time[size_tot];
					modec_out_.width(8);
					modec_out_ << "s";
				}
				else
				{
					for (int j = 1; j <= size_tot; ++j)
					{
						modec_out_.width(10);
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(4);
						modec_out_ << Cum_time[j];
						modec_out_.width(8);
						modec_out_ << "s";
					}
				}

				modec_out_ << '\n';
				modec_out_ << '\n';

				for (unsigned int i = 0; i < Nucl_size; ++i)
				{
					modec_out_.width(8);
					int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
					modec_out_ << nucl_id;
					int NZ = nucl_id / 10000;
					int NA = (nucl_id - NZ * 10000) / 10;
					stringstream ss;
					ss << NA;
					int NG = nucl_id - nucl_id / 10 * 10;
					string Nucl_Name;
					Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
					if (NG == 1) { Nucl_Name += "m"; }
					modec_out_.width(8);
					modec_out_ << Nucl_Name;
					if (print_mode_ == 0)
					{
						modec_out_.width(18);
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(9);
						modec_out_ << n_vector_[0][i + Nucl_size]*Avogadro_Constant*lamda[i] * TOXCI[i];
						modec_out_.width(18);
						if (n_vector_[size_tot][i + Nucl_size] < cutoff)
						{
							modec_out_ << 0.0;
						}
						else
						{
							modec_out_ << n_vector_[size_tot][i + Nucl_size]*Avogadro_Constant*lamda[i] * TOXCI[i];
						}

						modec_out_ << '\n';
					}
					if (print_mode_ == 1)
					{
						modec_out_.setf(ios::scientific | ios::uppercase);
						modec_out_.precision(9);
						for (int j = 0; j <= size_tot; ++j)
						{
							modec_out_.width(18);
							if (n_vector_[j][i + Nucl_size] < cutoff)
							{
								modec_out_ << 0.0;
							}
							else
							{
								modec_out_ << n_vector_[j][i + Nucl_size]*Avogadro_Constant*lamda[i] * TOXCI[i];
							}
						}
						modec_out_ << '\n';
					}
				}

				modec_out_.close();

			}
		}
	}
	else
	{
		if (input_filename_ != "BURNInp")
		{
			output_filename_ = input_filename_ + ".eql.concentration";
		}
		output_file_ = work_direc_ + output_filename_;
		modec_out_.open(output_file_);
		modec_out_ << "*********************************** Nuclide Equilibrium Concentrations, " << dens_unit_ << " ***********************************";
		modec_out_ << '\n';
		modec_out_ << '\n';
		modec_out_.setf(ios::left);

		modec_out_.width(8);
		modec_out_ << "";
		modec_out_.width(8);
		modec_out_ << "";
		modec_out_.width(18);
		modec_out_ << "ModecInitial";

		modec_out_.width(18);
		modec_out_ << "Equilibrium";

		modec_out_ << '\n';
///////////////////////  flux_ ///////////////////
		modec_out_.width(8);
		modec_out_ << "";
		modec_out_.width(8);
		modec_out_ << "flux_: ";

		modec_out_.width(10);
		modec_out_.setf(ios::scientific | ios::uppercase);
		modec_out_.precision(4);
		modec_out_ << flux_vector_[0];
		modec_out_.width(8);
		modec_out_ << "";

		modec_out_.width(10);
		modec_out_.setf(ios::scientific | ios::uppercase);
		modec_out_.precision(4);
		modec_out_ << flux_vector_[1];
		modec_out_.width(8);
		modec_out_ << "";

		modec_out_ << '\n';

///////////////////////  power ///////////////////
		modec_out_.width(8);
		modec_out_ << "";
		modec_out_.width(8);
		modec_out_ << "Power: ";

		modec_out_.width(10);
		modec_out_.setf(ios::scientific | ios::uppercase);
		modec_out_.precision(4);
		modec_out_ << power_vector_[0];
		modec_out_.width(8);
		modec_out_ << "";

		modec_out_.width(10);
		modec_out_.setf(ios::scientific | ios::uppercase);
		modec_out_.precision(4);
		modec_out_ << power_vector_[1];
		modec_out_.width(8);
		modec_out_ << "";
		modec_out_ << '\n';

		modec_out_ << '\n';
////////////////////////////////////////////////////////

		for (unsigned int i = 0; i < Nucl_size; ++i)
		{
			modec_out_.width(8);
			int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
			modec_out_ << nucl_id;
			int NZ = nucl_id / 10000;
			int NA = (nucl_id - NZ * 10000) / 10;
			stringstream ss;
			ss << NA;
			int NG = nucl_id - nucl_id / 10 * 10;
			string Nucl_Name;
			Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
			if (NG == 1) { Nucl_Name += "m"; }
			modec_out_.width(8);
			modec_out_ << Nucl_Name;

			double convert_coeff;

			if (dens_unit_ == "mol") { convert_coeff = 1.0; }
			else if (dens_unit_ == "g") { convert_coeff = double(NA); }
			else if (dens_unit_ == "kg") { convert_coeff = double(NA) / 1000.0; }
			else if (dens_unit_ == "atom") { convert_coeff = Avogadro_Constant; }
			else if (dens_unit_ == "atom/(barn-cm)") { convert_coeff = Avogadro_Constant * (1.0e-24); }

			modec_out_.width(18);
			modec_out_.setf(ios::scientific | ios::uppercase);
			modec_out_.precision(9);
			modec_out_ << n_vector_[0][i] * convert_coeff;
			modec_out_.width(18);
			if (n_vector_[1][i] < cutoff)
			{
				modec_out_ << 0.0;
			}
			else
			{
				modec_out_ << n_vector_[1][i] * convert_coeff;
			}

			modec_out_ << '\n';

		}

		modec_out_.close();

		if (if_print_activity_ == true)
		{

			vector<double > lamda;
			lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

			output_filename_ = input_filename_ + ".eql.radioactivity";
			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Equilibrium Radioactivity, Ci ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';
			modec_out_.setf(ios::left);

			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			modec_out_.width(18);
			modec_out_ << "Equilibrium";

			modec_out_ << '\n';

			modec_out_ << '\n';
			modec_out_ << '\n';

			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;
				
					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / (3.7e10);
					modec_out_.width(18);
					if (n_vector_[1][i] < cutoff)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[1][i]*Avogadro_Constant*lamda[i] / (3.7e10); // 单位为Ci
					}

					modec_out_ << '\n';

			}

			modec_out_.close();

		}

		if (if_print_decayenergy_ == true)
		{
			vector<double > lamda;
			lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

			vector<double > Q;
			Q = ModecNuclideLibrary.nuclide_library_vector_[3];

			output_filename_ = input_filename_ + ".eql.decayheat";
			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Equilibrium Decay Heat, Watts ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';
			modec_out_.setf(ios::left);

			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			modec_out_.width(18);
			modec_out_ << "Equilibrium";

			modec_out_ << '\n';

			modec_out_ << '\n';
			modec_out_ << '\n';


			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;

					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6;
					modec_out_.width(18);
					if (n_vector_[1][i] < cutoff)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[1][i]*Avogadro_Constant*lamda[i] * Q[i] * Electron_Coulomb*1.0e6; // 转换成单位 W
					}

					modec_out_ << '\n';
			}

			modec_out_.close();

		}

		if (if_print_ampc_ == true)
		{
			vector<double > lamda;
			lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

			vector<double > AMPC;
			AMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

			output_filename_ = input_filename_ + ".eql.ampc";
			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Equilibrium Radioactivity, m3-air/Bq ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';
			modec_out_.setf(ios::left);

			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			modec_out_.width(18);
			modec_out_ << "Equilibrium";

			modec_out_ << '\n';

			modec_out_ << '\n';
			modec_out_ << '\n';


			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;

					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
					modec_out_.width(18);
					if (n_vector_[1][i] < cutoff)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[1][i]*Avogadro_Constant*lamda[i] / 3.7e10 / AMPC[i];
					}

					modec_out_ << '\n';

			}

			modec_out_.close();

		}

		if (if_print_wmpc_ == true)
		{
			vector<double > lamda;
			lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

			vector<double > WMPC;
			WMPC = ModecNuclideLibrary.nuclide_library_vector_[4];

			output_filename_ = input_filename_ + ".eql.wmpc";
			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Equilibrium Radioactivity, m3-air/Bq ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';
			modec_out_.setf(ios::left);

			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			modec_out_.width(18);
			modec_out_ << "Equilibrium";

			modec_out_ << '\n';

			modec_out_ << '\n';
			modec_out_ << '\n';


			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;

					modec_out_.width(18);
					modec_out_.setf(ios::scientific | ios::uppercase);
					modec_out_.precision(9);
					modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
					modec_out_.width(18);
					if (n_vector_[1][i] < cutoff)
					{
						modec_out_ << 0.0;
					}
					else
					{
						modec_out_ << n_vector_[1][i]*Avogadro_Constant*lamda[i] / 3.7e10 / WMPC[i];
					}

					modec_out_ << '\n';
			
			}

			modec_out_.close();

		}

		if (if_print_toxicity_ == true)
		{
			vector<double > lamda;
			lamda = ModecNuclideLibrary.nuclide_library_vector_[1];

			vector<double > TOXIC;
			TOXIC = ModecNuclideLibrary.nuclide_library_vector_[10];

			output_filename_ = input_filename_ + ".eql.toxicity";
			output_file_ = work_direc_ + output_filename_;
			modec_out_.open(output_file_);
			modec_out_ << "*********************************** Nuclide Equilibrium Radiotoxicity, Sv ***********************************";
			modec_out_ << '\n';
			modec_out_ << '\n';
			modec_out_.setf(ios::left);

			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(8);
			modec_out_ << "";
			modec_out_.width(18);
			modec_out_ << "ModecInitial";

			modec_out_.width(18);
			modec_out_ << "Equilibrium";

			modec_out_ << '\n';

			modec_out_ << '\n';
			modec_out_ << '\n';


			for (unsigned int i = 0; i < Nucl_size; ++i)
			{
				modec_out_.width(8);
				int nucl_id = ModecNuclideLibrary.nuclide_list_[i];
				modec_out_ << nucl_id;
				int NZ = nucl_id / 10000;
				int NA = (nucl_id - NZ * 10000) / 10;
				stringstream ss;
				ss << NA;
				int NG = nucl_id - nucl_id / 10 * 10;
				string Nucl_Name;
				Nucl_Name = ModecNuclideLibrary.element_name_list_[NZ - 1] + ss.str();
				if (NG == 1) { Nucl_Name += "m"; }
				modec_out_.width(8);
				modec_out_ << Nucl_Name;

				modec_out_.width(18);
				modec_out_.setf(ios::scientific | ios::uppercase);
				modec_out_.precision(9);
				modec_out_ << n_vector_[0][i]*Avogadro_Constant*lamda[i] * TOXIC[i];
				modec_out_.width(18);
				if (n_vector_[1][i] < cutoff)
				{
					modec_out_ << 0.0;
				}
				else
				{
					modec_out_ << n_vector_[1][i]*Avogadro_Constant*lamda[i] * TOXIC[i];
				}

				modec_out_ << '\n';

			}

			modec_out_.close();

		}
	}
}