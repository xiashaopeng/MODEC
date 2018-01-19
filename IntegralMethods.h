#ifndef INTEGRAL_METHODS_H
#define INTEGRAL_METHODS_H

#include <iostream>
#include <vector>

using namespace std;
class GaussLegendreQuadrature    //  (-1,1)区间数值积分，用于添料率常数的处理
{
public:
	static int GL_order;
	static vector<vector<double> > gauss_legendre_weight_;	// 高斯权重
	static vector<vector<double> > gauss_legendre_abscissa_; // 高斯点
	static vector<vector<double> > Value_GL_weight();
	static vector<vector<double> > Value_GL_abscissa();
};

class GaussLaguerreQuadrature  // (0,infinite)区间数值积分，用于平衡态核素浓度的计算
{
public:
	static int GLa_order;
	static vector<vector<double> > GLa_weight;	// 高斯权重
	static vector<vector<double> > GLa_abscissa; // 高斯点
	static vector<vector<double> > Value_GLa_weight();
	static vector<vector<double> > Value_GLa_abscissa();
};

#endif
