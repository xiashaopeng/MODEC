#include "ModecClass.h"

using namespace std;


////////////////////////////////////////////读取Couple的ft33f001库///////////////////////////////////////////////
void ModecClass::ReadFromCouple()
{
	vector<int> nuclide_number_;// 文件中所有核素ID
	vector<int> DecNum; // 核素对应的衰变来源个数
	vector<int> XsNum; // 核素对应的中子反应来源个数
	vector<vector<int> > DecProdID; // 衰变来源核素ID，对应于nuclide_number_的顺序
	vector<vector<int> > XsProdID; // 中子核反应来源核素ID，对应于nuclide_number_的顺序
	vector<vector<double> > DecReactRate; // 矩阵非对角非零元衰变反应率，对应DecProdNum
	vector<vector<double> > XsReactRate; // 矩阵非对角非零元中子反应反应率，对应XsProdNum
	vector<double> DecDiagonal; // 矩阵对角元的衰变反应率
	vector<double> XsDiagonal; // 矩阵对角元的中子反应反应率

	int size_NuclNum, size_Tot;
	int fiss_tag, fiss_num;
	int itemp;
	double dtemp;
	char ctemp[100];

	ifstream ReadCouple;
	ReadCouple.open(couple_library_name_);
	if (!ReadCouple)
	{
		InfoMessage::ErrorMessage("Position: void ModecClass::ReadFromCouple; \n Error: cannot open COUPLE library file! go and check the existence of the lib file.", 1);
	}
	ReadCouple >> size_NuclNum;
	nuclide_number_.resize(size_NuclNum);
	DecNum.resize(size_NuclNum);
	XsNum.resize(size_NuclNum);
	DecProdID.resize(size_NuclNum);
	XsProdID.resize(size_NuclNum);
	DecDiagonal.resize(size_NuclNum);
	XsDiagonal.resize(size_NuclNum);
	DecReactRate.resize(size_NuclNum);
	XsReactRate.resize(size_NuclNum);

	ReadCouple >> fiss_tag;
	ReadCouple >> fiss_num;
	ReadCouple >> itemp;
	ReadCouple >> size_Tot;
	while (!ReadCouple.eof())
	{
		ReadCouple >> dtemp;
		if (dtemp == 20)
		{
			ReadCouple.getline(ctemp, 100);
			ReadCouple.getline(ctemp, 100);
			break;
		}
	}
	while (!ReadCouple.eof()) // 读取所有核素的ID并存如nuclide_number_向量中
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> nuclide_number_[i];
			}
			break;
		}
	}
	while (!ReadCouple.eof()) // 读取衰变来源的总个数，并存入DecNum向量中
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> DecNum[i];
			}
			break;
		}
	}
	while (!ReadCouple.eof()) // 读取每个核素来源的总个数，并存入XsNum向量中
	{
		ReadCouple >> itemp;
		if (itemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> XsNum[i];
			}
			break;
		}
	}


	while (!ReadCouple.eof()) // Index in array NUCL of the parent of the transition stored in matrix A
	{
		ReadCouple >> itemp;
		if (itemp == size_Tot)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				for (unsigned int j = 0; j < XsNum[i]; ++j)
				{
					ReadCouple >> itemp;
					if (j < DecNum[i])
					{
						DecProdID[i].push_back(itemp);
					}
					else
					{
						XsProdID[i].push_back(itemp);
					}
				}
			}
			break;
		}
	}


	while (!ReadCouple.eof()) // 读取每个核素的衰变常数Radioactive decay constants
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> DecDiagonal[i];
				if (DecDiagonal[i] > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[1][nuclid] = DecDiagonal[i];
					}			
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // Compressed transition matrix
	{
		ReadCouple >> dtemp;
		if (dtemp == size_Tot)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				for (unsigned int j = 0; j < XsNum[i]; ++j)
				{
					ReadCouple >> dtemp;
					if (j < DecNum[i])
					{
						DecReactRate[i].push_back(dtemp);
					}
					else
					{
						XsReactRate[i].push_back(dtemp);
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取每个核素的中子反应消失率Capture cross sections (total removal)
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> XsDiagonal[i];
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取裂变截面 Fission cross sections
	{
		ReadCouple >> dtemp;
		if (dtemp == fiss_num)
		{
			for (int i = fiss_tag; i < fiss_tag + fiss_num; ++i)
			{
				ReadCouple >> dtemp;
				if (dtemp > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[7][nuclid] = dtemp;
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取中子产生截面 Neutron production cross sections 
	{
		ReadCouple >> dtemp;
		if (dtemp == fiss_tag + fiss_num)
		{
			for (int i = 0; i < fiss_tag + fiss_num; ++i)
			{
				ReadCouple >> dtemp;
				if (dtemp > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[9][nuclid] = dtemp;
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取衰变热系数 Recoverable decay energy values
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> dtemp;
				if (dtemp > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[3][nuclid] = dtemp;
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // Fractions of Q-values associated with photon emission
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取AMPC系数 Radioactive concentration guide values for air
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> dtemp;
				if (dtemp > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[4][nuclid] = dtemp;
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取WMPC系数 Radioactive concentration guide values for air
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> dtemp;
				if (dtemp > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[5][nuclid] = dtemp;
					}
				}
			}
			break;
		}
	}

	ReadCouple.close();
	///////////////////////////////////////////////////////////////////////////


	//////////////////////////////// 构造矩阵TransMatrixDecay和TransMatrixCrossSection ////////////////////
	int _row, _col;
	double _val;
	//int mat_dim = ModecNuclideLibrary.nuclide_library_vector_[0].size();

	for (int i = 0; i < size_NuclNum; ++i)
	{
		int ReactID = nuclide_number_[i];
		_row = ModecNuclideLibrary.GetNuclIndex(ReactID);
		_val = DecDiagonal[i];
		if ((!TransMatrixDecay.JudgeExist(_row,_row, _val)) && _val>0.0)
		{
			TransMatrixDecay.AddElement(_row, _row, -_val);
		}

		for (unsigned int j = 0; j < DecProdID[i].size(); ++j)
		{
			int ProdcID = nuclide_number_[DecProdID[i][j] - 1];
			_col = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[DecProdID[i][j] - 1]);
			_val = DecReactRate[i][j];
			/*if (ReactID == 541251 && ProdcID == 541260)
			{
				continue;
			}*/

			if ((!TransMatrixDecay.JudgeExist(_row, _col, _val)) && _val > 0.0)
			{
				TransMatrixDecay.AddElement(_row, _col, _val);
			}
		}
	}

	if (if_tracking_stockage == true)
	{
		TransMatrixStockage = TransMatrixDecay;
	}

	/* 在截面矩阵建立的过程中，需要将n,g截面保存下来，用于功率、通量转换的计算 */
	
	for (int i = 0; i < size_NuclNum; ++i)
	{
		int ReactID = nuclide_number_[i];
		_row = ModecNuclideLibrary.GetNuclIndex(ReactID);
		_val = XsDiagonal[i];
		if ((!TransMatrixCrossSection.JudgeExist(_row, _row, _val)) && _val>0.0)
		{
			TransMatrixCrossSection.AddElement(_row, _row, -_val);
			ModecNuclideLibrary.nuclide_library_vector_[6][_row] = _val; // 将核素的 n,t 截面保存下来
		}
		for (unsigned int j = 0; j < XsProdID[i].size(); ++j)
		{
			int ProdcID = nuclide_number_[XsProdID[i][j] - 1];
			_col = ModecNuclideLibrary.GetNuclIndex(ProdcID);
			_val = XsReactRate[i][j];
			/*if (ReactID == 541240 && (ProdcID==561280 || ProdcID==541230))
			{
				_val = 0.0;
			}*/
			/*if (ReactID == 541251 && ProdcID== 541260)
			{
				continue;
			}*/
			if ((!TransMatrixCrossSection.JudgeExist(_row, _col, _val)) && _val > 0.0)
			{
				TransMatrixCrossSection.AddElement(_row, _col, _val);
				if (((ReactID - ProdcID == 10) && (ProdcID % 10 ==0)) || ReactID - ProdcID == 9) //包括激发态核素和基态核素的n,g反应，注意，不包含n,gx反应
				{
					ModecNuclideLibrary.nuclide_library_vector_[8][_col] = _val; // 将核素的 n,g 截面保存下来
				}
			}
		}
	}
}

void ModecClass::ReadFromCoupleForTta()
{
	vector<int> nuclide_number_;// 文件中所有核素ID
	vector<int> DecNum; // 核素对应的衰变来源个数
	vector<int> XsNum; // 核素对应的中子反应来源个数
	vector<vector<int> > DecProdID; // 衰变来源核素ID，对应于nuclide_number_的顺序
	vector<vector<int> > XsProdID; // 中子核反应来源核素ID，对应于nuclide_number_的顺序
	vector<vector<double> > DecReactRate; // 矩阵非对角非零元衰变反应率，对应DecProdNum
	vector<vector<double> > XsReactRate; // 矩阵非对角非零元中子反应反应率，对应XsProdNum
	vector<double> DecDiagonal; // 矩阵对角元的衰变反应率
	vector<double> XsDiagonal; // 矩阵对角元的中子反应反应率

	int size_NuclNum, size_Tot;
	int fiss_tag, fiss_num;
	int itemp;
	double dtemp;
	char ctemp[100];

	ifstream ReadCouple;
	ReadCouple.open(couple_library_name_);
	if (!ReadCouple)
	{
		InfoMessage::ErrorMessage("Position: void ModecClass::ReadFromCouple; \n Error: cannot open COUPLE library file! go and check the existence of the lib file.", 1);
	}
	ReadCouple >> size_NuclNum;
	nuclide_number_.resize(size_NuclNum);
	DecNum.resize(size_NuclNum);
	XsNum.resize(size_NuclNum);
	DecProdID.resize(size_NuclNum);
	XsProdID.resize(size_NuclNum);
	DecDiagonal.resize(size_NuclNum);
	XsDiagonal.resize(size_NuclNum);
	DecReactRate.resize(size_NuclNum);
	XsReactRate.resize(size_NuclNum);

	ReadCouple >> fiss_tag;
	ReadCouple >> fiss_num;
	ReadCouple >> itemp;
	ReadCouple >> size_Tot;
	while (!ReadCouple.eof())
	{
		ReadCouple >> dtemp;
		if (dtemp == 20)
		{
			ReadCouple.getline(ctemp, 100);
			ReadCouple.getline(ctemp, 100);
			break;
		}
	}
	while (!ReadCouple.eof()) // 读取所有核素的ID并存如nuclide_number_向量中
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> nuclide_number_[i];
			}
			break;
		}
	}
	while (!ReadCouple.eof()) // 读取衰变来源的总个数，并存入DecNum向量中
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> DecNum[i];
			}
			break;
		}
	}
	while (!ReadCouple.eof()) // 读取每个核素来源的总个数，并存入XsNum向量中
	{
		ReadCouple >> itemp;
		if (itemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> XsNum[i];
			}
			break;
		}
	}


	while (!ReadCouple.eof()) // Index in array NUCL of the parent of the transition stored in matrix A
	{
		ReadCouple >> itemp;
		if (itemp == size_Tot)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				for (unsigned int j = 0; j < XsNum[i]; ++j)
				{
					ReadCouple >> itemp;
					if (j < DecNum[i])
					{
						DecProdID[i].push_back(itemp);
					}
					else
					{
						XsProdID[i].push_back(itemp);
					}
				}
			}
			break;
		}
	}


	while (!ReadCouple.eof()) // 读取每个核素的衰变常数Radioactive decay constants
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> DecDiagonal[i];
				if (DecDiagonal[i] > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[1][nuclid] = DecDiagonal[i];
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // Compressed transition matrix
	{
		ReadCouple >> dtemp;
		if (dtemp == size_Tot)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				for (unsigned int j = 0; j < XsNum[i]; ++j)
				{
					ReadCouple >> dtemp;
					if (j < DecNum[i])
					{
						DecReactRate[i].push_back(dtemp);
					}
					else
					{
						XsReactRate[i].push_back(dtemp);
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取每个核素的中子反应消失率Capture cross sections (total removal)
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> XsDiagonal[i];
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取裂变截面 Fission cross sections
	{
		ReadCouple >> dtemp;
		if (dtemp == fiss_num)
		{
			for (int i = fiss_tag; i < fiss_tag + fiss_num; ++i)
			{
				ReadCouple >> dtemp;
				if (dtemp > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[7][nuclid] = dtemp;
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取中子产生截面 Neutron production cross sections 
	{
		ReadCouple >> dtemp;
		if (dtemp == fiss_tag + fiss_num)
		{
			for (int i = 0; i < fiss_tag + fiss_num; ++i)
			{
				ReadCouple >> dtemp;
				if (dtemp > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[9][nuclid] = dtemp;
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取衰变热系数 Recoverable decay energy values
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> dtemp;
				if (dtemp > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[3][nuclid] = dtemp;
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // Fractions of Q-values associated with photon emission
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取AMPC系数 Radioactive concentration guide values for air
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> dtemp;
				if (dtemp > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[4][nuclid] = dtemp;
					}
				}
			}
			break;
		}
	}

	while (!ReadCouple.eof()) // 读取WMPC系数 Radioactive concentration guide values for air
	{
		ReadCouple >> dtemp;
		if (dtemp == size_NuclNum)
		{
			for (int i = 0; i < size_NuclNum; ++i)
			{
				ReadCouple >> dtemp;
				if (dtemp > 0)
				{
					int nuclid = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[i]);
					if (nuclid >= 0)
					{
						ModecNuclideLibrary.nuclide_library_vector_[5][nuclid] = dtemp;
					}
				}
			}
			break;
		}
	}

	ReadCouple.close();
	///////////////////////////////////////////////////////////////////////////


	//////////////////////////////// 构造矩阵TransMatrixDecay和TransMatrixCrossSection ////////////////////
	int _row, _col;
	double _val;
	//int mat_dim = ModecNuclideLibrary.nuclide_library_vector_[0].size();

	for (int i = 0; i < size_NuclNum; ++i)
	{
		int ReactID = nuclide_number_[i];
		_row = ModecNuclideLibrary.GetNuclIndex(ReactID);
		_val = DecDiagonal[i];
		if ((!TtaMatrixDecay.JudgeExist(_row, _row, _val)) && _val>0.0)
		{
			TtaMatrixDecay.AddElementCCS(_row, _row, -_val);
		}

		for (unsigned int j = 0; j < DecProdID[i].size(); ++j)
		{
			int ProdcID = nuclide_number_[DecProdID[i][j] - 1];
			_col = ModecNuclideLibrary.GetNuclIndex(nuclide_number_[DecProdID[i][j] - 1]);
			_val = DecReactRate[i][j];
			/*if (ReactID == 541251 && ProdcID == 541260)
			{
			continue;
			}*/

			if ((!TtaMatrixDecay.JudgeExist(_row, _col, _val)) && _val > 0.0)
			{
				TtaMatrixDecay.AddElementCCS(_row, _col, _val);
			}
		}
	}

	//if (if_tracking_stockage == true)
	//{
	//	TransMatrixStockage = TtaMatrixDecay;
	//}

	/* 在截面矩阵建立的过程中，需要将n,g截面保存下来，用于功率、通量转换的计算 */

	for (int i = 0; i < size_NuclNum; ++i)
	{
		int ReactID = nuclide_number_[i];
		_row = ModecNuclideLibrary.GetNuclIndex(ReactID);
		_val = XsDiagonal[i];
		if ((!TtaMatrixCrossSection.JudgeExist(_row, _row, _val)) && _val>0.0)
		{
			TtaMatrixCrossSection.AddElementCCS(_row, _row, -_val);
			ModecNuclideLibrary.nuclide_library_vector_[6][_row] = _val; // 将核素的 n,t 截面保存下来
		}
		for (unsigned int j = 0; j < XsProdID[i].size(); ++j)
		{
			int ProdcID = nuclide_number_[XsProdID[i][j] - 1];
			_col = ModecNuclideLibrary.GetNuclIndex(ProdcID);
			_val = XsReactRate[i][j];
			/*if (ReactID == 541240 && (ProdcID==561280 || ProdcID==541230))
			{
			_val = 0.0;
			}*/
			/*if (ReactID == 541251 && ProdcID== 541260)
			{
			continue;
			}*/
			if ((!TtaMatrixCrossSection.JudgeExist(_row, _col, _val)) && _val > 0.0)
			{
				TtaMatrixCrossSection.AddElementCCS(_row, _col, _val);
				if (((ReactID - ProdcID == 10) && (ProdcID % 10 == 0)) || ReactID - ProdcID == 9) //包括激发态核素和基态核素的n,g反应，注意，不包含n,gx反应
				{
					ModecNuclideLibrary.nuclide_library_vector_[8][_col] = _val; // 将核素的 n,g 截面保存下来
				}
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////