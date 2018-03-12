/**
*  @file     ModecClass.h
*  @brief    MODEC程序流程执行控制模块
*  @author   夏少鹏
*  @email    xiashaopeng@sinap.ac.cn
*  @version  1.0.0
*  @date     01/30/2018
*/

#ifndef _RUN_MODEC_H
#define _RUN_MODEC_H
#include <sstream>
#include <iomanip>
#include "SolveTrans.h"
#include "IntegralMethods.h"
#include "tinyxml2.h"

/**
 * @brief MODEC的基础类
 * 包含了程序执行需要的所有参数的申明和定义
*/
class ModecClass {
  public:
    /**
    * @name MODEC计算基本参数
    * @{
    */
    string info_message_; 							///< 存储警告信息或者错误信息

    ifstream modec_inp_;							///< MODEC输入文件流
    ofstream modec_out_;							///< MODEC输出文件流

    string work_direc_;							///< 工作文件夹
    string input_filename_ = "modec.xml";					///< 缺省的输入卡文件名
    string input_file_;							///< 缺省的包含工作文件夹地址的输入卡文件名

    string output_filename_;  						///< 输出文件名。命令规则为输入卡文件名加上相应后缀
    string output_file_;							///< 包含工作文件夹地址的输出文件名
    string dens_unit_ = "mol";						///< 核素浓度单位: mol 为缺省值
    vector<string> time_unit_; 						///< 时间单位

    vector<double>  burnup_time_;						///< 每个子步的燃耗时间
    vector<int> substep_;							///< 子步总个数

    vector<int> evolution_mode_;					  	///< 燃耗模式: 0：衰变演化  1：定通量演化  2：定功率演化  3：自定义流程
    vector<double> evolution_value_;                  			///< 对应燃耗模式，存储对应通量或者功率数值

    int if_flow_mode_ = 0;							///< 判断是否求解流动燃耗模型: 0: 不流动模型; 1: 流动模型
    vector<double> residue_time_;						///< 流动燃耗模型下，每个燃耗区的滞留时间

    int solver_selection_ = 1;						///< 燃耗求解器。0: TTA方法   1：CRAM方法（缺省值）
    /** @} */

    /**
    * @name 燃耗数据库参数设置
    * @{
    */
    int lib_tag_; 								///< 数据库标志： 0:纯衰变且采用decay_library_name_; 1:采用depth_library_name_; 2:采用couple_library_name_
    string decay_library_name_, fission_yields_library_name_;		///< 衰变数据库名称和裂变产物库名称
    string depth_library_name_, couple_library_name_;			///< DEPTH数据库名称和ORIGEN-S数据库名称
    /** @} */

    /**
    * @name 燃耗计算存储数组
    * @{
    */
    vector<vector<double > >   n_vector_;  					///< 保存每个燃耗子步的计算结果，用于输出；
    vector<double >    power_vector_; 					///< 保存每个燃耗子步的功率结果，用于输出；
    vector<double >    flux_vector_;  					///< 保存每个燃耗子步的通量结果，用于输出；

    vector<double> kinf_vector_; 	 					///< 保存每个子步计算的kinf
    vector<double> fission_rate_vector_;  					///< 保存每个子步计算的中子产生率
    vector<double> absorption_rate_vector_; 				///< 保存每个子步的总的中子吸收率
    /** @} */

    /**
    * @name 输出格式参数
    * @{
    */
    int print_mode_ = 0;							///< 输出格式: 0：只输出最后一步的结果（缺省值）   1：输出每个子步的结果

    int if_print_kinf_ = 0;							///< 判断是否计算kinf。0: 不计算; 1: 计算
    int if_print_fission_rate_ = 0;						///< 判断是否计算中子产生率。0: 不计算; 1:计算总体的中子产生率并输出在浓度文件中; 2:计算每个核素的中子产生率并以单独文件输出
    int if_print_absorption_rate_ = 0;					///< 判断是否计算中子吸收率。0: 不计算; 1:计算总体的中子吸收率并输出在浓度文件中; 2:计算每个核素的中子吸收率并以单独文件输出
    int if_read_density_tag_ = 0;						///< 判断是否成功读入核素密度
    int if_read_time_tag_ = 0;						///< 判断是否成功读入燃耗时间
    int if_read_mode_tag_ = 0;						///< 判断是否读入燃耗模式

    bool if_calculate_equilibrium_ = false;					///< 判断是否计算平衡态

    bool if_print_activity_ = false;					///< 判断是否输出放射性活度结果
    bool if_print_decayenergy_ = false;					///< 判断是否输出衰变热结果
    bool if_print_ampc_ = false;						///< 判断是否输出AMPC毒性结果
    bool if_print_wmpc_ = false;						///< 判断是否输出WMPC毒性结果
    bool if_print_toxicity_ = false;					///< 判断是否输出Sv毒性结果

    bool if_print_stockage_activity_ = false;				///< 判断是否输出堆外存储区核素放射性活度
    bool if_print_stockage_decayenergy_ = false;				///< 判断是否输出堆外存储区核素衰变热
    bool if_print_stockage_ampc_ = false;					///< 判断是否输出堆外存储区核素AMPC毒性
    bool if_print_stockage_wmpc_ = false;					///< 判断是否输出堆外存储区核素WMPC毒性
    bool if_print_stockage_toxicity_ = false;				///< 判断是否输出堆外存储区核素Sv毒性
    /** @} */

    /**
    * @name 在线后处理参数设置
    * @{
    */
    bool if_continously_remove_ = false;					///< 判断是否有连续后处理，如果有则为真，缺省值为false
    bool if_constant_online_feeding_ = false;				///< 判断是否有固定添料率的在线添料，如果有则为真，缺省值为false

    bool if_tracking_stockage = false;					///< 判断是否追踪堆外存储区核素的演化，如果追踪则为真，缺省值为false

    int remove_group_number_;						///< 后处理元素组数
    vector<int> remove_group_vector_;					///< 各组的后处理元素个数
    vector<double> remove_rate_vector_;					///< 各组后处理元素的后处理速率
    vector<vector<int> > remove_element_vector_;				///< 各组后处理元素的原子序数Z
    /** @} */

    /**
    * @name 非定常在线后处理参数设置(已弃用)
    * @{
    */
    bool if_variable_feeding_ = false;					///< (弃用)判断是否有非定常添料率（添料率随燃料性质而变）在线添料，如果有则为真，缺省值为false
    bool if_keeping_eutectic_stable_ = false;				///< (弃用)判断是否保持熔盐共熔点稳定，如果保证则为真，缺省值为false
    int variable_feeding_group_num_;					///< (弃用)连续添料核素分组
    vector<int> variable_feeding_group_vector_;				///< (弃用)各组的添料核素个数
    double variable_feeding_ratio_;						///< (弃用)各组添料核素的比例份额
    vector<vector<int> > variable_feeding_nuclide_id_vector_;		///< (弃用)各组的添料核素的ID
    vector<vector<double> > variable_feeding_nuclide_ratio_vector_;		///< (弃用)各组的添料核素的浓度份额
    /** @} */

    /**
    * @name 连续添料参数设置
    * @{
    */
    int constant_feeding_calculation_methods_ = 2; 				///< 添料率常数的计算方法: =1表示采用数值积分方法；=2表示采用增广矩阵方法求解
    int constant_feeding_nuclide_num_;					///< 添料核素数量
    vector<int> constant_feeding_nuclide_id_vector_;			///< 添料核素ID向量
    vector<double> constant_feeding_rate_; 					///< 添料率向量，单位为mol/s

    vector<double > constant_feeding_vector_;				///< 用于数值积分法的添料率向量
    vector<double> gauss_legendre_weight_;					///< 高斯-勒让德求积权重
    vector<double> gauss_legendre_abscissa_;				///< 高斯-勒让德求积点
    /** @} */

    /**
    * @name 修正裂变产额参数
    * @{
    */
    int nearest_neighbor_;
    double yield_factor_;
    /** @} */

    /// 重金属核素的总的消失率，包括裂变消失和后处理移除消失，用于在线添料率的计算
    double heavy_nuclide_loss_rate_ = 0;

    /// 裂变产物的总的产生率，用裂变产生率减去后处理移除速率，用于保证共熔点的计算，在Li盐体系下，保证Li加上裂变产物mol浓度不变
    double fission_products_production_rate_ = 0;

  public:
    /**
    * @brief MODEC外部接口函数
    */
    void RunModec(int argc, char *argv[]) {
        ModecInitial(argc,argv);
        BuildSpMat();

        ////////// 先给出初始时刻的通量和功率 //////////
        if (evolution_mode_[0] == 1) {
            ModecNuclideLibrary.flux_ = evolution_value_[0];
        }
        if (evolution_mode_[0] == 2) {
            ModecNuclideLibrary.specified_power_ = evolution_value_[0];
        }
        ModecNuclideLibrary.CalculateFlux(evolution_mode_[0]);
        flux_vector_.push_back(ModecNuclideLibrary.flux_);
        power_vector_.push_back(ModecNuclideLibrary.specified_power_);
        ///////////////////////////////////////////////
        InfoMessage::ends.push_back(clock());

        if (if_calculate_equilibrium_ == false) {
            ModecProcedure();
        } else {
            CalEquilibrium(evolution_mode_[0]);
        }
        InfoMessage::ends.push_back(clock());
        ModecOutputXml(); // 测试XML输出
        InfoMessage::ends.push_back(clock());
        return;
    };


  private:
    /**
    * @name 基于CRAM方法的燃耗方程系数矩阵
    * @{
    */
    SpMat TransMatrixDecay;								///< 用于CRAM方法的衰变系数存储矩阵
    SpMat TransMatrixPureDecay;							///< 用于流动计算的纯衰变矩阵
    SpMat TransMatrixCrossSection;							///< 用于CRAM方法的截面数据存储矩阵
    SpMat TransMatrixFissionYields; 						///< 用于CRAM方法的裂变产物存储矩阵（只用于读取DepthLib）

    SpMat TransMatrixReprocess;							///< 用于追踪后处理的堆外核素演化时，存储堆内到堆外的迁移矩阵
    SpMat TransMatrixStockage;							///< 用于追踪后处理的堆外核素演化时，存储堆外衰变的衰变矩阵
    /** @} */

    /**
    * @name 基于TTA方法的燃耗方程系数矩阵
    * @{
    */
    SparseMatrixMCS TtaMatrixDecay;							///< 用于TTA方法的衰变系数邻接稀疏矩阵
    SparseMatrixMCS TtaMatrixCrossSection;						///< 用于TTA方法的截面数据邻接稀疏矩阵
    SparseMatrixMCS TtaMatrixFissionYields; 					///< 用于TTA方法的裂变产额数据邻接稀疏矩阵
    /** @} */

    SolveTrans Solver; 								///< 选择求解器：0: TTA; 1: CRAM

    NuclLibrary ModecNuclideLibrary; 						///< MODEC燃耗数据库

  private:
    /**
    * @brief 输入卡读入模块
    */
    void ModecInitial(int argc, char **argv);

    /**
    * @brief 燃耗计算控制模块
    */
    void ModecProcedure();

    /**
    * @brief 结果输出模块
    */
    void ModecOutput();

    /**
    * @brief 结果输出模块xml
    */
    void ModecOutputXml();

    /**
    * @brief 燃耗矩阵构建控制模块
    */
    void BuildSpMat();

    /**
    * @brief 平衡态计算控制模块
    */
    void CalEquilibrium(int mode);

    /**
    * @brief 整体流程控制模块
    */
    void Evolution(int mode, double time, int subtime);

    /**
    * @brief 矩阵输出函数
    * 输出燃耗矩阵到文件中，用于研究之用
    */
    void TransitionMatrixOutput(SpMat matrix) { // 指定矩阵元素输出到文件中
        int matrix_dimension = ModecNuclideLibrary.nuclide_number_;

        ofstream matrix_output;
        matrix_output.open("transition_matrix.out");
        for (int _row = 0; _row < matrix_dimension; ++_row) {
            for (int _col = 0; _col < matrix_dimension; ++_col) {
                matrix_output<< matrix.Element(_row,_col)<<" ";
            }
            matrix_output << '\n';
        }
        exit(0);
    }

  private: // 矩阵构建函数

    /**
    * @name 读取MODEC自带数据库构建矩阵
    * @{
    */

    /**
    * @brief 读取自带数据库构建CRAM用衰变矩阵
    */
    void DecayToSpMat();
    /**
    * @brief 读取自带数据库构建TTA用衰变矩阵
    */
    void DecayToSpMatForTta();
    /**
    * @brief 读取TRITON输出文件构建CRAM用燃耗矩阵
    */
    void XSfromTriton();
    /**
    * @brief 读取自带裂变份额数据库
    */
    void CalculateEffectiveFissionYields();

    /** @} */

    /**
    * @name 读取DEPTH数据库构建矩阵
    * @{
    */
    /**
    * @brief 读取DEPTH数据库构建CRAM用燃耗矩阵
    */
    void ReadFromDepthLib();
    /**
    * @brief 读取DEPTH数据库构建TTA用燃耗矩阵
    */
    void ReadFromDepthLibForTta(); ///< TTA方法的矩阵读取
    /**
    * @brief 读取DEPTH数据库构建CRAM用裂变份额矩阵
    */
    void ConstructFissionYieldsSpMat();
    /**
    * @brief 读取自带数据库构建TTA用裂变份额矩阵
    */
    void ConstructFissionYieldsSpMatForTta();
    /** @} */

    /**
    * @name 读取COUPLE加工的ORIGENS格式数据库构建矩阵
    * @{
    */
    /**
    * @brief 读取COUPLE数据库构建CRAM用燃耗矩阵
    */
    void ReadFromCouple();
    /**
    * @brief 读取COUPLE数据库构建TTA用燃耗矩阵
    */
    void ReadFromCoupleForTta();
    /** @} */

    /**
    * @name 在线后处理以及连续添料处理模块
    * @{
    */
    /**
    * @brief 处理在线后处理系数
    */
    void AddOnlineReprocessingCoeffi();
    /**
    * @brief 处理连续添料模块(弃用)
    */
    void ContinuouslyFeeding();
    /** @} */
};
#endif
