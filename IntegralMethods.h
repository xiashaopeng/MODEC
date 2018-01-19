#ifndef INTEGRAL_METHODS_H
#define INTEGRAL_METHODS_H

#include <iostream>
#include <vector>

using namespace std;
class GaussLegendreQuadrature    //  (-1,1)������ֵ���֣����������ʳ����Ĵ���
{
public:
	static int GL_order;
	static vector<vector<double> > gauss_legendre_weight_;	// ��˹Ȩ��
	static vector<vector<double> > gauss_legendre_abscissa_; // ��˹��
	static vector<vector<double> > Value_GL_weight();
	static vector<vector<double> > Value_GL_abscissa();
};

class GaussLaguerreQuadrature  // (0,infinite)������ֵ���֣�����ƽ��̬����Ũ�ȵļ���
{
public:
	static int GLa_order;
	static vector<vector<double> > GLa_weight;	// ��˹Ȩ��
	static vector<vector<double> > GLa_abscissa; // ��˹��
	static vector<vector<double> > Value_GLa_weight();
	static vector<vector<double> > Value_GLa_abscissa();
};

#endif
