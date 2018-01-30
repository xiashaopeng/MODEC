#include <iostream>
#include "NuclList.h"
using namespace std;

const double Electron_Coulomb = 1.6021766208E-19;//1.6021766208E-19; // 电子电量        Source: 2014 CODATA  
const double Avogadro_Constant = 6.022140857E+23; // 阿伏加德罗常数   Source: 2014 CODATA  

void NuclLibrary::CalculateFlux(int mode)
{
	//int mode;
	if (mode == 0)
	{
		specified_power_ = 0.0;
		flux_ = 0.0;
	}

	/////////// Constant Flux ////////////////////////
	if (mode == 1)
	{
		double flux = CalculatePowerKernel();
		specified_power_ = flux_ * Electron_Coulomb * flux / 1.0e24 ; // 计算每一步的功率
	}
	////////////////////////////////////////////////////

	/////////// Constant Power ////////////////////////
	if (mode == 2)
	{
		double flux = CalculatePowerKernel();
		flux_ = 1.0 / Electron_Coulomb * specified_power_/ flux*1.0e24; // 此处最后*1.0e24的原因是TRITON给出的截面单位为barn
	}
	////////////////////////////////////////////////////
}

double NuclLibrary::CalculatePowerKernel()
{
	double flux = 0.0;
	for (unsigned int i = 0; i < nuclide_number_; ++i)// 此处默认为采用mol作为核素浓度单位
	{
		if (nuclide_list_[i] == 902300)
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (190.00*nuclide_library_vector_[7][i] + 5.010*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 902320)
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (189.21*nuclide_library_vector_[7][i] + 4.786*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 902330)//Th233
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (190.00*nuclide_library_vector_[7][i] + 6.080*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 912310)//pa231
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (190.00*nuclide_library_vector_[7][i] + 5.660*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 912330)//Pa233
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (189.10*nuclide_library_vector_[7][i] + 5.197*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 922320)//U232
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (200.00*nuclide_library_vector_[7][i] + 5.930*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 922330)//U233
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (191.29*nuclide_library_vector_[7][i] + 6.841*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 922340)//U234
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (190.30*nuclide_library_vector_[7][i] + 5.297*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 922350)//U235
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (194.02*nuclide_library_vector_[7][i] + 6.545*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 922360)//U236
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (192.80*nuclide_library_vector_[7][i] + 5.124*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 922380)//U238
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (198.12*nuclide_library_vector_[7][i] + 4.804*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 932370)//Np237
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (195.10*nuclide_library_vector_[7][i] + 5.490*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 932390)//Np239
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (200.00*nuclide_library_vector_[7][i] + 4.970*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 942380)//Pu238
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (197.80*nuclide_library_vector_[7][i] + 5.550*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 942390)//Pu239
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (200.05*nuclide_library_vector_[7][i] + 6.533*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 942400)//Pu240
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (199.79*nuclide_library_vector_[7][i] + 5.241*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 942410)//Pu241
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (202.22*nuclide_library_vector_[7][i] + 6.301*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 942420)//Pu242
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (200.62*nuclide_library_vector_[7][i] + 5.071*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 942430)//Pu243
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (200.00*nuclide_library_vector_[7][i] + 6.020*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 952410)//Am241
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (202.30*nuclide_library_vector_[7][i] + 5.529*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 952421)//Am242m
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (202.29*nuclide_library_vector_[7][i] + 6.426*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 952430)//Am243
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (202.10*nuclide_library_vector_[7][i] + 5.363*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 962440)//Cm244
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (200.00*nuclide_library_vector_[7][i] + 6.451*nuclide_library_vector_[8][i]);
		}
		else if (nuclide_list_[i] == 962450)//Cm245
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * (200.00*nuclide_library_vector_[7][i] + 6.110*nuclide_library_vector_[8][i]);
		}

		else if (nuclide_list_[i] == 10010)//H1
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 2.225 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 50100)//B10
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 2.790 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 80160)//O16
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 4.143 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 260560)//Fe56
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 7.600 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 280580)//Ni58
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 9.020 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 400900)//Zr90
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 7.203 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 400910)//Zr91
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 8.635 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 400920)//Zr92
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 6.758 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 400960)//Zr96
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 5.571 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 420950)//Mo95
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 9.154 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 430950)//Tc95
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 7.710 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 441010)//Ru101
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 9.216 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 441030)//Ru103
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 6.999 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 441050)//Ru105
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 7.094 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 471090)//Ag109
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 6.825 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 541310)//Xe131
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 8.936 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 541350)//Xe135
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 7.880 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 551330)//Cs133
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 6.704 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 551340)//Cs134
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 6.550 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 601430)//Nd143
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 7.817 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 601450)//Nd145
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 7.565 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 611470)//Pm147
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 5.900 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 611480)//Pm148
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 7.266 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 611481)//Pm148m
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 7.266 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 621470)//Sm147
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 8.140 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 621490)//Sm149
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 7.982 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 621500)//Sm150
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 5.596 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 621510)//Sm151
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 8.258 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 621520)//Sm152
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 5.867 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 631530)//Eu153
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 6.444 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 631540)//Eu154
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 8.167 * nuclide_library_vector_[8][i];
		}
		else if (nuclide_list_[i] == 631550)//Eu155
		{
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 6.490 * nuclide_library_vector_[8][i];
		}
		else
		{
			if (nuclide_list_[i] > 890000)
			{
				/*double RE;
				int Z = nuclide_list_[i] / 10000;
				int A = (nuclide_list_[i] - Z * 10000) / 10;
				RE = 1.29927E-3*(pow(Z, 2.0)*pow(A, 0.5)) + 33.12;
				flux += nuclide_library_vector_[0][i] * Avogadro_Constant * RE *nuclide_library_vector_[7][i];*/
				flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 200.0 *nuclide_library_vector_[7][i];
			}
			flux += nuclide_library_vector_[0][i] * Avogadro_Constant * 5.0 * nuclide_library_vector_[8][i];
			//flux += nuclide_library_vector_[0][i] * Avogadro_Constant * nuclide_library_vector_[3][i] * nuclide_library_vector_[1][i]*1.0E+24;
		}
	}

	return flux;
}
