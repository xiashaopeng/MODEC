#include "ModecClass.h"

using namespace std;


////////////////////////////////////////////��ȡCouple��ft33f001��///////////////////////////////////////////////
void ModecClass::ReadFromCouple()
{
	vector<int> nuclide_number_;// �ļ������к���ID
	vector<int> DecNum; // ���ض�Ӧ��˥����Դ����
	vector<int> XsNum; // ���ض�Ӧ�����ӷ�Ӧ��Դ����
	vector<vector<int> > DecProdID; // ˥����Դ����ID����Ӧ��nuclide_number_��˳��
	vector<vector<int> > XsProdID; // ���Ӻ˷�Ӧ��Դ����ID����Ӧ��nuclide_number_��˳��
	vector<vector<double> > DecReactRate; // ����ǶԽǷ���Ԫ˥�䷴Ӧ�ʣ���ӦDecProdNum
	vector<vector<double> > XsReactRate; // ����ǶԽǷ���Ԫ���ӷ�Ӧ��Ӧ�ʣ���ӦXsProdNum
	vector<double> DecDiagonal; // ����Խ�Ԫ��˥�䷴Ӧ��
	vector<double> XsDiagonal; // ����Խ�Ԫ�����ӷ�Ӧ��Ӧ��

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
	while (!ReadCouple.eof()) // ��ȡ���к��ص�ID������nuclide_number_������
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
	while (!ReadCouple.eof()) // ��ȡ˥����Դ���ܸ�����������DecNum������
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
	while (!ReadCouple.eof()) // ��ȡÿ��������Դ���ܸ�����������XsNum������
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


	while (!ReadCouple.eof()) // ��ȡÿ�����ص�˥�䳣��Radioactive decay constants
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

	while (!ReadCouple.eof()) // ��ȡÿ�����ص����ӷ�Ӧ��ʧ��Capture cross sections (total removal)
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

	while (!ReadCouple.eof()) // ��ȡ�ѱ���� Fission cross sections
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

	while (!ReadCouple.eof()) // ��ȡ���Ӳ������� Neutron production cross sections 
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

	while (!ReadCouple.eof()) // ��ȡ˥����ϵ�� Recoverable decay energy values
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

	while (!ReadCouple.eof()) // ��ȡAMPCϵ�� Radioactive concentration guide values for air
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

	while (!ReadCouple.eof()) // ��ȡWMPCϵ�� Radioactive concentration guide values for air
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


	//////////////////////////////// �������TransMatrixDecay��TransMatrixCrossSection ////////////////////
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

	/* �ڽ���������Ĺ����У���Ҫ��n,g���汣�����������ڹ��ʡ�ͨ��ת���ļ��� */
	
	for (int i = 0; i < size_NuclNum; ++i)
	{
		int ReactID = nuclide_number_[i];
		_row = ModecNuclideLibrary.GetNuclIndex(ReactID);
		_val = XsDiagonal[i];
		if ((!TransMatrixCrossSection.JudgeExist(_row, _row, _val)) && _val>0.0)
		{
			TransMatrixCrossSection.AddElement(_row, _row, -_val);
			ModecNuclideLibrary.nuclide_library_vector_[6][_row] = _val; // �����ص� n,t ���汣������
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
				if (((ReactID - ProdcID == 10) && (ProdcID % 10 ==0)) || ReactID - ProdcID == 9) //��������̬���غͻ�̬���ص�n,g��Ӧ��ע�⣬������n,gx��Ӧ
				{
					ModecNuclideLibrary.nuclide_library_vector_[8][_col] = _val; // �����ص� n,g ���汣������
				}
			}
		}
	}
}

void ModecClass::ReadFromCoupleForTta()
{
	vector<int> nuclide_number_;// �ļ������к���ID
	vector<int> DecNum; // ���ض�Ӧ��˥����Դ����
	vector<int> XsNum; // ���ض�Ӧ�����ӷ�Ӧ��Դ����
	vector<vector<int> > DecProdID; // ˥����Դ����ID����Ӧ��nuclide_number_��˳��
	vector<vector<int> > XsProdID; // ���Ӻ˷�Ӧ��Դ����ID����Ӧ��nuclide_number_��˳��
	vector<vector<double> > DecReactRate; // ����ǶԽǷ���Ԫ˥�䷴Ӧ�ʣ���ӦDecProdNum
	vector<vector<double> > XsReactRate; // ����ǶԽǷ���Ԫ���ӷ�Ӧ��Ӧ�ʣ���ӦXsProdNum
	vector<double> DecDiagonal; // ����Խ�Ԫ��˥�䷴Ӧ��
	vector<double> XsDiagonal; // ����Խ�Ԫ�����ӷ�Ӧ��Ӧ��

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
	while (!ReadCouple.eof()) // ��ȡ���к��ص�ID������nuclide_number_������
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
	while (!ReadCouple.eof()) // ��ȡ˥����Դ���ܸ�����������DecNum������
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
	while (!ReadCouple.eof()) // ��ȡÿ��������Դ���ܸ�����������XsNum������
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


	while (!ReadCouple.eof()) // ��ȡÿ�����ص�˥�䳣��Radioactive decay constants
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

	while (!ReadCouple.eof()) // ��ȡÿ�����ص����ӷ�Ӧ��ʧ��Capture cross sections (total removal)
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

	while (!ReadCouple.eof()) // ��ȡ�ѱ���� Fission cross sections
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

	while (!ReadCouple.eof()) // ��ȡ���Ӳ������� Neutron production cross sections 
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

	while (!ReadCouple.eof()) // ��ȡ˥����ϵ�� Recoverable decay energy values
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

	while (!ReadCouple.eof()) // ��ȡAMPCϵ�� Radioactive concentration guide values for air
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

	while (!ReadCouple.eof()) // ��ȡWMPCϵ�� Radioactive concentration guide values for air
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


	//////////////////////////////// �������TransMatrixDecay��TransMatrixCrossSection ////////////////////
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

	/* �ڽ���������Ĺ����У���Ҫ��n,g���汣�����������ڹ��ʡ�ͨ��ת���ļ��� */

	for (int i = 0; i < size_NuclNum; ++i)
	{
		int ReactID = nuclide_number_[i];
		_row = ModecNuclideLibrary.GetNuclIndex(ReactID);
		_val = XsDiagonal[i];
		if ((!TtaMatrixCrossSection.JudgeExist(_row, _row, _val)) && _val>0.0)
		{
			TtaMatrixCrossSection.AddElementCCS(_row, _row, -_val);
			ModecNuclideLibrary.nuclide_library_vector_[6][_row] = _val; // �����ص� n,t ���汣������
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
				if (((ReactID - ProdcID == 10) && (ProdcID % 10 == 0)) || ReactID - ProdcID == 9) //��������̬���غͻ�̬���ص�n,g��Ӧ��ע�⣬������n,gx��Ӧ
				{
					ModecNuclideLibrary.nuclide_library_vector_[8][_col] = _val; // �����ص� n,g ���汣������
				}
			}
		}
	}
}
//////////////////////////////////////////////////////////////////////////////////////////////////////