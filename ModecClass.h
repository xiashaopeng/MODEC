#ifndef _RUN_MODEC_H
#define _RUN_MODEC_H
//#include <iostream>
//#include <fstream>
#include "SolveTrans.h"
//#include "SparseMatrix.h"
//#include "NuclList.h"
#include "IntegralMethods.h"
class ModecClass
{
public: 
	string info_message_; // 存储警告或者错误信息

	ifstream modec_inp_;							// 读取BURN输入卡
	ofstream modec_out_;							// BURN输出流

	string work_direc_;//工作文件夹
	string input_filename_ = "modec.input";
	string input_file_;			// 缺省的BURN输入卡名称

	string output_filename_;  // 输出文件名，为输入卡文件名加上相应后缀
	string output_file_;		// 缺省的BURN输入卡名称
	string dens_unit_ = "mol";					// 核素浓度单位 mol 为缺省值
	vector<string> time_unit_; 					// 时间单位 d 为缺省值


	vector<double>  burnup_time_;				// 每个子步的燃耗时间
	vector<int> substep_;						// 子步个数
	
	vector<int> evolution_mode_;					  // 燃耗模式， 0：衰变演化    1：定通量演化     2：定功率演化     3：自定义流程
	vector<double> evolution_value_;                  // 对应燃耗模式，0.0代表纯衰变，其余对应通量(10)或者功率(2)

	int solver_selection_ = 1;						// 0: TTA方法   1：CRAM方法（缺省值）
	int print_mode_ = 0;							// 0：只输出最后一步的结果（缺省值）   1：输出每个子步的结果
	string decay_library_name_, fission_yields_library_name_;
	string depth_library_name_, couple_library_name_;					// 衰变数据库和裂变产物数据库名称
	int lib_tag_; // 数据库标志 lib_tag_ = 0纯衰变且采用decay_library_name_
				  //           lib_tag_ = 1采用depth_library_name_
				  //           lib_tag_ = 2采用couple_library_name_	

	vector<vector<double > >   n_vector_;  // 保存每个燃耗子步的计算结果，用于输出；
	vector<double >    power_vector_;  // 保存每个燃耗子步的功率结果，用于输出；
	vector<double >    flux_vector_;  // 保存每个燃耗子步的通量结果，用于输出；


	//  三个bool变量，分别用于判断是否计算kinf,裂变反应率和中子吸收率
	int if_print_kinf_ = 0;
	int if_print_fission_rate_ = 0;
	int if_print_absorption_rate_ = 0;

	vector<double> kinf_vector_;  // 保存每个子步用反应率计算的kinf
	vector<double> fission_rate_vector_;  // 保存每个子步用反应率计算的裂变反应率
	vector<double> absorption_rate_vector_; // 保存每步的总的中子吸收率

	/// 判断input是否读入成功 ///
	int if_read_density_tag_ = 0;
	int if_read_time_tag_ = 0;
	int if_read_mode_tag_ = 0;

	bool if_calculate_equilibrium_ = false;

	bool if_print_activity_ = false;
	bool if_print_decayenergy_ = false;
	bool if_print_ampc_ = false;
	bool if_print_wmpc_ = false;
	bool if_print_toxicity_ = false; // 增加Sv毒性输出信息

	bool if_print_stockage_activity_ = false;
	bool if_print_stockage_decayenergy_ = false;
	bool if_print_stockage_ampc_ = false;
	bool if_print_stockage_wmpc_ = false;
	bool if_print_stockage_toxicity_ = false; // 增加Sv毒性输出信息

	//////////////////////////////   后处理相关变量参数  /////////////////////////////////

	bool if_continously_remove_ = false;				// 是否有连续后处理，如果有则为真，缺省值为false

	bool if_variable_feeding_ = false;					// 是否有非定常添料率（添料率随燃料性质而变）在线添料，如果有则为真，缺省值为false

	bool if_constant_online_feeding_ = false;			// 是否有固定添料率的在线添料，如果有则为真，缺省值为false

	bool if_keeping_eutectic_stable_ = false;		// 是否保持熔盐共熔点稳定，如果保证则为真，缺省值为false

	bool if_tracking_stockage = false;			// 是否追踪堆外提取裂变产物和Pa233的演化，如果追踪则为真，缺省值为false

	int remove_group_number_;						// 后处理元素分组
	vector<int> remove_group_vector_;					// 各组的后处理元素个数
	vector<double> remove_rate_vector_;					// 各组后处理元素的后处理速率
	vector<vector<int> > remove_element_vector_;		// 各组后处理元素的原子序数 Z

	//         非定常添料率相关变量 -- 该模块已弃用       //
	int variable_feeding_group_num_;						// 连续添料核素分组
	vector<int> variable_feeding_group_vector_;					// 各组的添料核素个数
	double variable_feeding_ratio_;						// 各组添料核素的比例份额
	vector<vector<int> > variable_feeding_nuclide_id_vector_;		// 各组的添料核素的ID
	vector<vector<double> > variable_feeding_nuclide_ratio_vector_;		// 各组的添料核素的浓度份额
	

//////////////////// 固定添料率相关变量 ////////////////////////////
	int constant_feeding_calculation_methods_ = 2; // 添料率常数的计算方法，=1表示采用数值积分方法；=2表示采用增广矩阵方法求解
	int constant_feeding_nuclide_num_;
	vector<int> constant_feeding_nuclide_id_vector_;
	vector<double> constant_feeding_rate_; // 单位为mol/s

	vector<double > constant_feeding_vector_;
	vector<double> gauss_legendre_weight_;		//高斯-勒让德求积权重
	vector<double> gauss_legendre_abscissa_;		//高斯-勒让德求积点
	//vector<double> GL_weight_eql;		//高斯-勒让德求积权重
	//vector<double> GL_abscissa_eql;	//高斯-勒让德求积点


	// 重金属核素的总的消失率，包括裂变消失和后处理移除消失，用于在线添料率的计算
	double heavy_nuclide_loss_rate_ = 0; 

	// 裂变产物的总的产生率，用裂变产生率减去后处理移除速率，用于保证共熔点的计算，在Li盐体系下，保证Li加上裂变产物mol浓度不变
	double fission_products_production_rate_ = 0; 

	////////////  修正裂变份额相关参数  ////////////
	int nearest_neighbor_;
	double yield_factor_;

public:
	void RunModec(int argc, char *argv[])
	{
		ModecInitial(argc,argv);		
		BuildSpMat();

		////////// 先给出初始时刻的通量和功率 //////////
		if (evolution_mode_[0] == 1)
		{
			ModecNuclideLibrary.flux_ = evolution_value_[0];
		}
		if (evolution_mode_[0] == 2)
		{
			ModecNuclideLibrary.specified_power_ = evolution_value_[0];
		}
		ModecNuclideLibrary.CalculateFlux(evolution_mode_[0]);
		flux_vector_.push_back(ModecNuclideLibrary.flux_);
		power_vector_.push_back(ModecNuclideLibrary.specified_power_);
		///////////////////////////////////////////////
		InfoMessage::ends.push_back(clock());

		if (if_calculate_equilibrium_ == false)
		{
			ModecProcedure();
		}
		else
		{
			CalEquilibrium(evolution_mode_[0]);
		}
		InfoMessage::ends.push_back(clock());
		ModecOutput();
		InfoMessage::ends.push_back(clock());
		return;
	};


//////////////////////////////////////////////////////   ModecClass 私有成员函数及私有成员变量   ////////////////////////////////////////////////////
private:
	SpMat TransMatrixDecay, TransMatrixCrossSection, TransMatrixFissionYields; // 核素演化计算中必须包含衰变，核素截面为可选项，如果是纯衰变则不需要

	SpMat TransMatrixReprocess, TransMatrixStockage; // 当追踪堆外核素演化时，需要用到的两个矩阵，分别对应大矩阵的左下和右下的分块矩阵

	SparseMatrixMCS TtaMatrixDecay, TtaMatrixCrossSection, TtaMatrixFissionYields; // 用于TTA方法的稀疏矩阵存储格式

	SolveTrans Solver; // 选择求解器：CRAM

	NuclLibrary ModecNuclideLibrary; //数据库建立 
private:
	void ModecInitial(int argc, char *argv[]);
	void ModecProcedure();
	void ModecOutput();

	void BuildSpMat();

	void CalEquilibrium(int mode);

	void Evolution(int mode, double time, int subtime);


	void TransitionMatrixOutput(SpMat matrix) // 指定矩阵元素输出到文件中
	{
		int matrix_dimension = ModecNuclideLibrary.nuclide_number_;

		ofstream matrix_output;
		matrix_output.open("transition_matrix.out");
		for (int _row = 0; _row < matrix_dimension; ++_row)
		{
			for (int _col = 0; _col < matrix_dimension; ++_col)
			{
				matrix_output<< matrix.Element(_row,_col)<<" ";
			}
			matrix_output << '\n';
		}
		exit(0);
	}

private: // 矩阵构建函数

	/////////////////////////////// 读取MODEC本身数据库 //////////////////////////////////////////////
	void DecayToSpMat();
	void DecayToSpMatForTta();
	void XSfromTriton();
	void CalculateEffectiveFissionYields();
	
	///////////////////////////////  读取Depth程序数据库：ORIGEN2/ORIGENS  /////////////////////////////////////
	void ReadFromDepthLib();
	void ReadFromDepthLibForTta(); // TTA方法的矩阵读取
	void ConstructFissionYieldsSpMat();
	void ConstructFissionYieldsSpMatForTta();

	///////////////////////////////  读取Couple数据库：ORIGENS  ///////////////////////////////////////////////
	void ReadFromCouple();
	void ReadFromCoupleForTta();

	//////////////////////////////////在线后处理和连续添料的处理模块////////////////////////////////////////////////
	void AddOnlineReprocessingCoeffi();
	void ContinuouslyFeeding();

};
#endif
