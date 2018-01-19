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
	string info_message_; // �洢������ߴ�����Ϣ

	ifstream modec_inp_;							// ��ȡBURN���뿨
	ofstream modec_out_;							// BURN�����

	string work_direc_;//�����ļ���
	string input_filename_ = "modec.input";
	string input_file_;			// ȱʡ��BURN���뿨����

	string output_filename_;  // ����ļ�����Ϊ���뿨�ļ���������Ӧ��׺
	string output_file_;		// ȱʡ��BURN���뿨����
	string dens_unit_ = "mol";					// ����Ũ�ȵ�λ mol Ϊȱʡֵ
	vector<string> time_unit_; 					// ʱ�䵥λ d Ϊȱʡֵ


	vector<double>  burnup_time_;				// ÿ���Ӳ���ȼ��ʱ��
	vector<int> substep_;						// �Ӳ�����
	
	vector<int> evolution_mode_;					  // ȼ��ģʽ�� 0��˥���ݻ�    1����ͨ���ݻ�     2���������ݻ�     3���Զ�������
	vector<double> evolution_value_;                  // ��Ӧȼ��ģʽ��0.0����˥�䣬�����Ӧͨ��(10)���߹���(2)

	int solver_selection_ = 1;						// 0: TTA����   1��CRAM������ȱʡֵ��
	int print_mode_ = 0;							// 0��ֻ������һ���Ľ����ȱʡֵ��   1�����ÿ���Ӳ��Ľ��
	string decay_library_name_, fission_yields_library_name_;
	string depth_library_name_, couple_library_name_;					// ˥�����ݿ���ѱ�������ݿ�����
	int lib_tag_; // ���ݿ��־ lib_tag_ = 0��˥���Ҳ���decay_library_name_
				  //           lib_tag_ = 1����depth_library_name_
				  //           lib_tag_ = 2����couple_library_name_	

	vector<vector<double > >   n_vector_;  // ����ÿ��ȼ���Ӳ��ļ����������������
	vector<double >    power_vector_;  // ����ÿ��ȼ���Ӳ��Ĺ��ʽ�������������
	vector<double >    flux_vector_;  // ����ÿ��ȼ���Ӳ���ͨ����������������


	//  ����bool�������ֱ������ж��Ƿ����kinf,�ѱ䷴Ӧ�ʺ�����������
	int if_print_kinf_ = 0;
	int if_print_fission_rate_ = 0;
	int if_print_absorption_rate_ = 0;

	vector<double> kinf_vector_;  // ����ÿ���Ӳ��÷�Ӧ�ʼ����kinf
	vector<double> fission_rate_vector_;  // ����ÿ���Ӳ��÷�Ӧ�ʼ�����ѱ䷴Ӧ��
	vector<double> absorption_rate_vector_; // ����ÿ�����ܵ�����������

	/// �ж�input�Ƿ����ɹ� ///
	int if_read_density_tag_ = 0;
	int if_read_time_tag_ = 0;
	int if_read_mode_tag_ = 0;

	bool if_calculate_equilibrium_ = false;

	bool if_print_activity_ = false;
	bool if_print_decayenergy_ = false;
	bool if_print_ampc_ = false;
	bool if_print_wmpc_ = false;
	bool if_print_toxicity_ = false; // ����Sv���������Ϣ

	bool if_print_stockage_activity_ = false;
	bool if_print_stockage_decayenergy_ = false;
	bool if_print_stockage_ampc_ = false;
	bool if_print_stockage_wmpc_ = false;
	bool if_print_stockage_toxicity_ = false; // ����Sv���������Ϣ

	//////////////////////////////   ������ر�������  /////////////////////////////////

	bool if_continously_remove_ = false;				// �Ƿ������������������Ϊ�棬ȱʡֵΪfalse

	bool if_variable_feeding_ = false;					// �Ƿ��зǶ��������ʣ���������ȼ�����ʶ��䣩�������ϣ��������Ϊ�棬ȱʡֵΪfalse

	bool if_constant_online_feeding_ = false;			// �Ƿ��й̶������ʵ��������ϣ��������Ϊ�棬ȱʡֵΪfalse

	bool if_keeping_eutectic_stable_ = false;		// �Ƿ񱣳����ι��۵��ȶ��������֤��Ϊ�棬ȱʡֵΪfalse

	bool if_tracking_stockage = false;			// �Ƿ�׷�ٶ�����ȡ�ѱ�����Pa233���ݻ������׷����Ϊ�棬ȱʡֵΪfalse

	int remove_group_number_;						// ����Ԫ�ط���
	vector<int> remove_group_vector_;					// ����ĺ���Ԫ�ظ���
	vector<double> remove_rate_vector_;					// �������Ԫ�صĺ�������
	vector<vector<int> > remove_element_vector_;		// �������Ԫ�ص�ԭ������ Z

	//         �Ƕ�����������ر��� -- ��ģ��������       //
	int variable_feeding_group_num_;						// �������Ϻ��ط���
	vector<int> variable_feeding_group_vector_;					// ��������Ϻ��ظ���
	double variable_feeding_ratio_;						// �������Ϻ��صı����ݶ�
	vector<vector<int> > variable_feeding_nuclide_id_vector_;		// ��������Ϻ��ص�ID
	vector<vector<double> > variable_feeding_nuclide_ratio_vector_;		// ��������Ϻ��ص�Ũ�ȷݶ�
	

//////////////////// �̶���������ر��� ////////////////////////////
	int constant_feeding_calculation_methods_ = 2; // �����ʳ����ļ��㷽����=1��ʾ������ֵ���ַ�����=2��ʾ����������󷽷����
	int constant_feeding_nuclide_num_;
	vector<int> constant_feeding_nuclide_id_vector_;
	vector<double> constant_feeding_rate_; // ��λΪmol/s

	vector<double > constant_feeding_vector_;
	vector<double> gauss_legendre_weight_;		//��˹-���õ����Ȩ��
	vector<double> gauss_legendre_abscissa_;		//��˹-���õ������
	//vector<double> GL_weight_eql;		//��˹-���õ����Ȩ��
	//vector<double> GL_abscissa_eql;	//��˹-���õ������


	// �ؽ������ص��ܵ���ʧ�ʣ������ѱ���ʧ�ͺ����Ƴ���ʧ���������������ʵļ���
	double heavy_nuclide_loss_rate_ = 0; 

	// �ѱ������ܵĲ����ʣ����ѱ�����ʼ�ȥ�����Ƴ����ʣ����ڱ�֤���۵�ļ��㣬��Li����ϵ�£���֤Li�����ѱ����molŨ�Ȳ���
	double fission_products_production_rate_ = 0; 

	////////////  �����ѱ�ݶ���ز���  ////////////
	int nearest_neighbor_;
	double yield_factor_;

public:
	void RunModec(int argc, char *argv[])
	{
		ModecInitial(argc,argv);		
		BuildSpMat();

		////////// �ȸ�����ʼʱ�̵�ͨ���͹��� //////////
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


//////////////////////////////////////////////////////   ModecClass ˽�г�Ա������˽�г�Ա����   ////////////////////////////////////////////////////
private:
	SpMat TransMatrixDecay, TransMatrixCrossSection, TransMatrixFissionYields; // �����ݻ������б������˥�䣬���ؽ���Ϊ��ѡ�����Ǵ�˥������Ҫ

	SpMat TransMatrixReprocess, TransMatrixStockage; // ��׷�ٶ�������ݻ�ʱ����Ҫ�õ����������󣬷ֱ��Ӧ���������º����µķֿ����

	SparseMatrixMCS TtaMatrixDecay, TtaMatrixCrossSection, TtaMatrixFissionYields; // ����TTA������ϡ�����洢��ʽ

	SolveTrans Solver; // ѡ���������CRAM

	NuclLibrary ModecNuclideLibrary; //���ݿ⽨�� 
private:
	void ModecInitial(int argc, char *argv[]);
	void ModecProcedure();
	void ModecOutput();

	void BuildSpMat();

	void CalEquilibrium(int mode);

	void Evolution(int mode, double time, int subtime);


	void TransitionMatrixOutput(SpMat matrix) // ָ������Ԫ��������ļ���
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

private: // ���󹹽�����

	/////////////////////////////// ��ȡMODEC�������ݿ� //////////////////////////////////////////////
	void DecayToSpMat();
	void DecayToSpMatForTta();
	void XSfromTriton();
	void CalculateEffectiveFissionYields();
	
	///////////////////////////////  ��ȡDepth�������ݿ⣺ORIGEN2/ORIGENS  /////////////////////////////////////
	void ReadFromDepthLib();
	void ReadFromDepthLibForTta(); // TTA�����ľ����ȡ
	void ConstructFissionYieldsSpMat();
	void ConstructFissionYieldsSpMatForTta();

	///////////////////////////////  ��ȡCouple���ݿ⣺ORIGENS  ///////////////////////////////////////////////
	void ReadFromCouple();
	void ReadFromCoupleForTta();

	//////////////////////////////////���ߺ�����������ϵĴ���ģ��////////////////////////////////////////////////
	void AddOnlineReprocessingCoeffi();
	void ContinuouslyFeeding();

};
#endif
