#ifndef _SPARSE_MATRIX_H
#define _SPARSE_MATRIX_H

//#include <iostream>
//#include <vector>
//#include <cstdlib>
//#include <math.h>
#include "Complex.h"
#include "ErrorMessage.h"
using namespace std;
const int cNuclNum = 1693;
// �������л��ģʽ���洢����CRAM�㷨��ȼ��ϡ����󣬲�����ֱ�ӷ�LU�ֽ�������Է��������⡣

class SparseMatrixMCS // ˫����洢ϡ�����
{
public:
	vector< double > diagonal_val_; // �Խ�Ԫ��
	vector<vector< double> > col_val_; // ������ѹ����ÿ�е�Ԫ��ֵ
	vector<vector<int> > col_index_; // ������ѹ����ÿ�е�Ԫ�ص�λ��
	int element_amount_ = 0; // ����ǶԽǷ���Ԫ����
	int spmat_dimen_;
public:
	SparseMatrixMCS(int nuclide_number_) // ��ϡ�����Ľ�������ʼ��ϡ�����������
	{
		spmat_dimen_ = nuclide_number_;
		diagonal_val_.resize(nuclide_number_);
		col_val_.resize(nuclide_number_);
		col_index_.resize(nuclide_number_);
	}

	SparseMatrixMCS() // ��ϡ�����Ľ�������ʼ��ϡ�����������
	{
		spmat_dimen_ = cNuclNum;
		diagonal_val_.resize(cNuclNum);
		col_val_.resize(cNuclNum);
		col_index_.resize(cNuclNum);
	}

	void Reset()
	{
		element_amount_ = 0;
		for (int i = 0; i < spmat_dimen_; ++i)
		{
			diagonal_val_[i] = 0.0;
			col_val_[i].resize(0);
			col_index_[i].resize(0);
		}
	}

	void Resize(int nuclide_number_)
	{
		spmat_dimen_ = nuclide_number_;
		diagonal_val_.resize(nuclide_number_);
		col_val_.resize(nuclide_number_);
		col_index_.resize(nuclide_number_);
	}

	void AddElementCCS(const int &_row,const int &_col,const double &_val) // ϡ��������ѹ���洢
	{
		if (_row == -1 || _col == -1)
		{
			InfoMessage::ErrorMessage("Position: void SparseMatrixMCS::AddElementCCS; \n Warning: index of row or column equal -1.", 0);
			return;
		}
		if (_row == _col)
		{
			diagonal_val_[_row] = _val;
			return;
		}
		else
		{
			/* �����жϸ�λ���Ƿ��Ѿ�����Ԫ�أ������δ���ڣ�ֱ����ӣ�����Ѿ����ڣ����ӷ����� */
			int PL = -1; // ��ʶ���������ж��Ƿ����Ԫ��
			for (unsigned int i = 0; i < col_index_[_col].size(); ++i)
			{
				if (col_index_[_col][i] == _row)
				{
					PL = i; // ��ʾ�������к���
					break;
				}
			}
			if (PL == -1)
			{
				col_index_[_col].push_back(_row);
				col_val_[_col].push_back(_val);
				element_amount_++;
			}
			else
			{
				col_val_[_col][PL] += _val;
			}
		}
	}

	//void RemoveElementCCS(int _row, int _col) // ��ѹ���洢ϡ������ɾ��Ԫ����
	//{
	//	if (_row == -1 || _col == -1)
	//	{
	//		InfoMessage::ErrorMessage("����ľ������꣺-1", 0);
	//		return;
	//	}
	//	if (_row == _col)
	//	{
	//		diagonal_val_[_row] = 0.0;
	//		return;
	//	}
	//	else
	//	{
	//		/* �����жϸ�λ���Ƿ��Ѿ�����Ԫ�أ������δ���ڣ�ֱ����ӣ�����Ѿ����ڣ����ӷ����� */
	//		int PL = -1; // ��ʶ���������ж��Ƿ����Ԫ��
	//		for (unsigned int i = 0; i < col_index_[_col].size(); ++i)
	//		{
	//			if (col_index_[_col][i] == _row)
	//			{
	//				PL = i; // ��ʾ�������к���
	//				break;
	//			}
	//		}
	//		if (PL == -1)
	//		{
	//			return;
	//		}
	//		else
	//		{
	//			col_val_[_col][PL] += _val;
	//		}
	//	}
	//}

	double ElementCCS(const int &_row, const int & _col)
	{
		if (_row == _col)
		{
			return diagonal_val_[_row];
		}

		for (unsigned int i = 0; i < col_index_[_col].size(); ++i)
		{
			if (col_index_[_col][i] == _row)
			{
				return col_val_[_col][i];
			}
		}
		return 0.0;
	}

	vector<int> RowListIndex(const int & _col)
	{
		if (col_index_[_col].size() == 0)
		{
			return vector<int>{-1};
		}
		vector<int> row_list_index;
		for (unsigned int i = 0; i < col_index_[_col].size(); ++i)
		{
			row_list_index.push_back(col_index_[_col][i]);
		}
		return row_list_index;
	}

	vector<double> RowListValue(const int & _col)
	{
		if (col_index_[_col].size() == 0)
		{
			InfoMessage::ErrorMessage("Position: vector<double> SparseMatrixMCS::RowListValue; \n Error: no elements in the required column.",1);
		}
		vector<double> row_list_val;
		for (unsigned int i = 0; i < col_index_[_col].size(); ++i)
		{
			row_list_val.push_back(col_val_[_col][i]);
		}
		return row_list_val;
	}

	SparseMatrixMCS operator*(const double & time)
	{
		int size = diagonal_val_.size();
		SparseMatrixMCS matrix(*this); // ���ƹ���matrix
		for (int i = 0; i < size; ++i)
		{
			matrix.diagonal_val_[i] = diagonal_val_[i] * time;
		}
		for (int i = 0; i < size; ++i)
		{
			int size1 = col_val_[i].size();
			for (int j = 0; j < size1; ++j)
			{
				matrix.col_val_[i][j] = col_val_[i][j] * time;
			}
		}
		return matrix;
	}

	SparseMatrixMCS operator+(const SparseMatrixMCS& matrix)
	{
		SparseMatrixMCS summatrix(matrix); // ��matrix������this��ֵ��ԭ���ǿ��ǵ��ѱ����Ҫ����˥�����
		int size = summatrix.diagonal_val_.size();

		for (int i = 0; i < size; ++i)
		{
			summatrix.diagonal_val_[i] += diagonal_val_[i];

			int size1 = summatrix.col_val_[i].size();
			int size2 = col_val_[i].size();
			for (int j = 0; j < size2; ++j)
			{
				int jj = 0;
				for (jj = 0; jj < size1; ++jj)
				{
					if (summatrix.col_index_[i][jj] == col_index_[i][j])
					{
						summatrix.col_val_[i][jj] += col_val_[i][j];
						break;
					}
				}
				if (jj >= size1)
				{
					summatrix.col_index_[i].push_back(col_index_[i][j]);
					summatrix.col_val_[i].push_back(col_val_[i][j]);
					summatrix.element_amount_++;
				}
			}
		}
		return summatrix;
	}

	bool JudgeExist(const int & _row, const int & _col, const double & _val) //�жϾ����λ���Ƿ����Ԫ�أ��������Ԫ�صĻ�������TRUE
	{
		if (_row == -1 || _col == -1)
		{
			return false;
		}
		bool exist_element = false;
		if (_row == _col)
		{
			if (diagonal_val_[_row] != 0.0)
			{
				exist_element = true;
			}
		}
		else // ��������
		{
			int size = col_index_[_col].size();
			int i;
			for (i = 0; i < size; ++i)
			{
				if (col_index_[_col][i] == _row)
				{
					break;
				}
			}
			if (i < size)
			{
				exist_element = true;
			}
		}

		return exist_element;
	}
};

class SpMat // ����CRAM�㷨��ȼ�ľ���洢��ʽ
	/* �Խ�Ԫ�����洢��������Ԫ�ز�����ѹ���洢��������Ԫ�ز�����ѹ���洢*/
{
public:
	vector<Complex > diagonal_val_; // value of diagonal element
	vector<vector<Complex > > lower_val_; // value of lower element
	vector<vector<int> > lower_index_; // ������Ԫ�ص�λ��
	vector<vector<Complex > > upper_val_; // value of upper element
	vector<vector<int> > upper_index_; // ������Ԫ�ص�λ��

	vector<Complex > ce_; // �������Ԫ�������Խ�Ԫ����ֵ��Ϊ�˷������CRAM��⣬Ԫ�����������趨Ϊ������

	int element_amount_ = 0; // ����ǶԽǷ���Ԫ����
	int spmat_dimen_; // ����ά��

public:
	SpMat()
	{
		diagonal_val_.resize(cNuclNum);
		lower_val_.resize(cNuclNum);
		lower_index_.resize(cNuclNum);
		upper_val_.resize(cNuclNum);
		upper_index_.resize(cNuclNum);
		spmat_dimen_ = cNuclNum;
	};

	SpMat(int nuclide_number_)
	{
		diagonal_val_.resize(nuclide_number_);
		lower_val_.resize(nuclide_number_);
		lower_index_.resize(nuclide_number_);
		upper_val_.resize(nuclide_number_);
		upper_index_.resize(nuclide_number_);
		spmat_dimen_ = nuclide_number_;
	}

	void Reset()
	{
		element_amount_ = 0;
		ce_.resize(0);
		for (int i = 0; i < spmat_dimen_; ++i)
		{
			diagonal_val_[i] = 0.0;
			lower_val_[i].resize(0);
			lower_index_[i].resize(0);
			upper_val_[i].resize(0);
			upper_index_[i].resize(0);
		}
	}

	void Resize(int nuclide_number_)
	{
		diagonal_val_.resize(nuclide_number_);
		lower_val_.resize(nuclide_number_);
		lower_index_.resize(nuclide_number_);
		upper_val_.resize(nuclide_number_);
		upper_index_.resize(nuclide_number_);
		spmat_dimen_ = nuclide_number_;
	}

	void AddElement(const int &_row, const int & _col, const double &_val)
	{
		if (_row == -1 || _col == -1)
		{
			InfoMessage::ErrorMessage("Position: void SpMat::AddElement; \n Warning: index of row or column equal -1.", 0);
			return;
		}
		if (_row == _col)
		{
			diagonal_val_[_row] += _val;
			//element_amount_ ++;
		}
		else if (_row > _col) // �����Ǿ���
		{
			/* �����жϸ�λ���Ƿ��Ѿ�����Ԫ�أ������δ���ڣ�ֱ����ӣ�����Ѿ����ڣ����ӷ����� */
			int PL = -1; // ��ʶ���������ж��Ƿ����Ԫ��
			for (unsigned int i = 0; i < lower_index_[_col].size(); ++i)
			{
				if (lower_index_[_col][i] == _row)
				{ 
					PL = i; // ��ʾ�������к���
					break;
				}
			}
			if (PL == -1)
			{
				lower_index_[_col].push_back(_row);
				lower_val_[_col].push_back(_val);
				element_amount_++;
			}
			else
			{
				lower_val_[_col][PL] += _val;
			}
		}
		else // �����Ǿ���
		{
			/* �����жϸ�λ���Ƿ��Ѿ�����Ԫ�أ������δ���ڣ�ֱ����ӣ�����Ѿ����ڣ����ӷ����� */
			int PL = -1; // ��ʶ���������ж��Ƿ����Ԫ��
			for (unsigned int i = 0; i < upper_index_[_row].size(); ++i)
			{
				if (upper_index_[_row][i] == _col)
				{
					PL = i; // ��ʾ�������к���
					break;
				}
			}
			if (PL == -1)
			{
				upper_index_[_row].push_back(_col);
				upper_val_[_row].push_back(_val);
				element_amount_++;
			}
			else
			{
				upper_val_[_row][PL] += _val;
			}
		}
	}
	
	double Element(const int & _row, const int & _col)
	{
		if (_row == _col)
		{
			return diagonal_val_[_row].real_;
		}
		else if (_row > _col) // ��������
		{
			int size = lower_index_[_col].size();
			int i;
			for (i = 0; i < size; ++i)
			{
				if (lower_index_[_col][i] == _row)
				{
					break;
				}
			}
			if (i < size)
			{
				return lower_val_[_col][i].real_;
			}
			else
			{
				return 0.0;
			}
		}
		else // ��������
		{
			int size = upper_index_[_row].size();
			int i;
			for (i = 0; i < size; ++i)
			{
				if (upper_index_[_row][i] == _col)
				{
					break;
				}
			}
			if (i < size)
			{
				return upper_val_[_row][i].real_;
			}
			else
			{
				return 0.0;
			}
		}
	}

	bool JudgeExist(const int & _row, const int & _col) //�жϾ����λ���Ƿ����Ԫ�أ��������Ԫ�صĻ�������TRUE
	{
		bool exist_element = false;
		if (_row == _col)
		{
			if (diagonal_val_[_row] != 0.0)
			{
				exist_element = true;
			}
		}
		else if (_row > _col) // ��������
		{
			int size = lower_index_[_col].size();
			int i;
			for (i = 0; i < size; ++i)
			{
				if (lower_index_[_col][i] == _row)
				{
					break;
				}
			}
			if (i < size)
			{
				exist_element = true;
			}
		}
		else // ��������
		{
			int size = upper_index_[_row].size();
			int i;
			for (i = 0; i < size; ++i)
			{
				if (upper_index_[_row][i] == _col)
				{
					break;
				}
			}
			if (i < size)
			{
				exist_element = true;
			}
		}
		return exist_element;
	}

	bool JudgeExist(const int & _row, const int & _col,const double &_val) //�жϾ����λ���Ƿ����Ԫ�أ��������Ԫ�صĻ�������TRUE
	{
		if (_row == -1 || _col == -1)
		{
			return false;
		}
		bool exist_element = false;
		if (_row == _col)
		{
			if (diagonal_val_[_row] != 0.0)
			{
				exist_element = true;
			}
		}
		else if (_row > _col) // ��������
		{
			int size = lower_index_[_col].size();
			int i;
			for (i = 0; i < size; ++i)
			{
				if (lower_index_[_col][i] == _row)
				{
					break;
				}
			}
			if (i < size && lower_val_[_col][i] == _val)
			{
				exist_element = true;
			}
		}
		else // ��������
		{
			int size = upper_index_[_row].size();
			int i;
			for (i = 0; i < size; ++i)
			{
				if (upper_index_[_row][i] == _col)
				{
					break;
				}
			}
			if (i < size && upper_val_[_row][i] == _val)
			{
				exist_element = true;
			}
		}
		return exist_element;
	}

	SpMat Identity()
	{
		SpMat I(spmat_dimen_);
		for (int i=0;i<spmat_dimen_;++i)
		{
			I.diagonal_val_[i] = 1;
		}
		
		return I;
	}

	SpMat operator*(const Complex &time)
	{
		int size = diagonal_val_.size();
		SpMat matrix(*this); // ���ƹ���matrix
		for (int i = 0; i < size; ++i)
		{
			matrix.diagonal_val_[i] = diagonal_val_[i] * time;
		}
		for (int i = 0; i < size; ++i)
		{
			int size1 = lower_val_[i].size();
			int size2 = upper_val_[i].size();
			for (int j = 0; j < size1; ++j)
			{
				matrix.lower_val_[i][j] = lower_val_[i][j] * time;
			}
			for (int j = 0; j < size2; ++j)
			{
				matrix.upper_val_[i][j] = upper_val_[i][j] * time;
			}
		}
		return matrix;
	}

	SpMat operator+(const SpMat& matrix)
	{
		SpMat summatrix(matrix); // ��matrix������this��ֵ��ԭ���ǿ��ǵ��ѱ����Ҫ����˥�����
		int size = summatrix.diagonal_val_.size();
		for (int i = 0; i < size; ++i)
		{
			summatrix.diagonal_val_[i] += diagonal_val_[i];
		}
		for (int i = 0; i < size; ++i)
		{
			int size1 = summatrix.lower_val_[i].size();
			int size2 = lower_val_[i].size();
			for (int j = 0; j < size2; ++j)
			{
				int jj = 0;
				for (jj = 0; jj < size1; ++jj)
				{
					if (summatrix.lower_index_[i][jj] == lower_index_[i][j])
					{
						summatrix.lower_val_[i][jj] += lower_val_[i][j];
						break;
					}
				}
				if (jj >= size1)
				{
					summatrix.lower_index_[i].push_back(lower_index_[i][j]);
					summatrix.lower_val_[i].push_back(lower_val_[i][j]);
					summatrix.element_amount_++;
				}
			}

			int size3 = summatrix.upper_val_[i].size();
			int size4 = upper_val_[i].size();
			for (int j = 0; j < size4; ++j)
			{
				int jj = 0;
				for (jj = 0; jj < size3; ++jj)
				{
					if (summatrix.upper_index_[i][jj] == upper_index_[i][j])
					{
						summatrix.upper_val_[i][jj] += upper_val_[i][j];
						break;
					}
				}
				if (jj >= size3)
				{
					summatrix.upper_index_[i].push_back(upper_index_[i][j]);
					summatrix.upper_val_[i].push_back(upper_val_[i][j]);
					summatrix.element_amount_++;
				}
			}
		}
		return summatrix;
	}

	SpMat& operator=(const SpMat& matrix)
	{
		diagonal_val_ = matrix.diagonal_val_;
		lower_val_ = matrix.lower_val_;
		lower_index_ = matrix.lower_index_;
		upper_val_ = matrix.upper_val_;
		upper_index_ = matrix.upper_index_;
		spmat_dimen_ = matrix.spmat_dimen_;
		return *this;
	}

	SpMat ConvertCram(const Complex &theta,const double &time)
	{
		SpMat CramSpMat(*this);
		int size = CramSpMat.diagonal_val_.size();
		for (int i = 0; i < size; ++i)
		{
			CramSpMat.diagonal_val_[i] = CramSpMat.diagonal_val_[i] * time - theta;
		}
		
		for (int i = 0; i < size; ++i)
		{
			int size2 = CramSpMat.lower_val_[i].size();
			for (int j = 0; j < size2; ++j)
			{
				CramSpMat.lower_val_[i][j] = CramSpMat.lower_val_[i][j] * time;
			}
		}

		for (int i = 0; i < size; ++i)
		{
			int size3 = CramSpMat.upper_val_[i].size();
			for (int j = 0; j < size3; ++j)
			{
				CramSpMat.upper_val_[i][j] = CramSpMat.upper_val_[i][j] * time;
			}
		}
		return CramSpMat;
	}

	void AddSpMat(const SpMat& matrix)
	{
		int size = diagonal_val_.size();
		for (int i = 0; i < size; ++i)
		{
			diagonal_val_[i] += matrix.diagonal_val_[i];
		}
		for (int i = 0; i < size; ++i)
		{
			int size1 = lower_val_[i].size();
			int size2 = matrix.lower_val_[i].size();
			for (int j = 0; j < size2; ++j)
			{
				int jj = 0;
				for (jj = 0; jj < size1; ++jj)
				{
					if (lower_index_[i][jj] == matrix.lower_index_[i][j])
					{
						lower_val_[i][jj] += matrix.lower_val_[i][j];
						break;
					}
				}
				if (jj >= size1)
				{
					lower_index_[i].push_back(matrix.lower_index_[i][j]);
					lower_val_[i].push_back(matrix.lower_val_[i][j]);
					element_amount_++;
				}
			}

			int size3 = upper_val_[i].size();
			int size4 = matrix.upper_val_[i].size();
			for (int j = 0; j < size4; ++j)
			{
				int jj = 0;
				for (jj = 0; jj < size3; ++jj)
				{
					if (upper_index_[i][jj] == matrix.upper_index_[i][j])
					{
						upper_val_[i][jj] += matrix.upper_val_[i][j];
						break;
					}
				}
				if (jj >= size3)
				{
					upper_index_[i].push_back(matrix.upper_index_[i][j]);
					upper_val_[i].push_back(matrix.upper_val_[i][j]);
					element_amount_++;
				}
			}
		}
	}

	void AddSpMat(const SpMat& matrix,const Complex &flux)
	{
		int size = diagonal_val_.size();
		for (int i = 0; i < size; ++i)
		{
			diagonal_val_[i] += matrix.diagonal_val_[i] * flux;
		}
		for (int i = 0; i < size; ++i)
		{
			int size1 = lower_val_[i].size();
			int size2 = matrix.lower_val_[i].size();
			for (int j = 0; j < size2; ++j)
			{
				int jj = 0;
				for (jj = 0; jj < size1; ++jj)
				{
					if (lower_index_[i][jj] == matrix.lower_index_[i][j])
					{
						lower_val_[i][jj] += matrix.lower_val_[i][j] * flux;
						break;
					}
				}
				if (jj >= size1)
				{
					lower_index_[i].push_back(matrix.lower_index_[i][j]);
					lower_val_[i].push_back(matrix.lower_val_[i][j] * flux);
					element_amount_++;
				}
			}

			int size3 = upper_val_[i].size();
			int size4 = matrix.upper_val_[i].size();
			for (int j = 0; j < size4; ++j)
			{
				int jj = 0;
				for (jj = 0; jj < size3; ++jj)
				{
					if (upper_index_[i][jj] == matrix.upper_index_[i][j])
					{
						upper_val_[i][jj] += matrix.upper_val_[i][j] * flux;
						break;
					}
				}
				if (jj >= size3)
				{
					upper_index_[i].push_back(matrix.upper_index_[i][j]);
					upper_val_[i].push_back(matrix.upper_val_[i][j] * flux);
					element_amount_++;
				}
			}
		}
	}
/////////////////////////////////////////////////////////////////////////////////////
///////////////    ����LU�ֽ�������ϡ��洢����һ��Ŀ�����Ż������ٶ� /////////////////
public:
		/* ����LU�ֽ���Ҫ�ı��� */
		vector<int> IRC;
		vector<int> ICFR;
		vector<int> LUP;
		
public:

	void SymbolLUElimination() // ����LU�ֽ�,ͬʱ����ϡ������ce_
	{
		int nuclide_number_ = spmat_dimen_;
		vector<vector<int> > UpperIndex_temp(upper_index_), LowerIndex_temp(lower_index_);
		vector<vector<int> > LUP_Track(3, vector<int>(1)); // LUP�б��׷���б�
		for (int i = 0; i < nuclide_number_; ++i)
		{
			int rowsize = LowerIndex_temp[i].size();
			int colsize = UpperIndex_temp[i].size();
			for (int rr = 0; rr < rowsize; ++rr)
			{
				int rowIndex = LowerIndex_temp[i][rr];
				for (int cc = 0; cc < colsize; ++cc)
				{
					int colIndex = UpperIndex_temp[i][cc];
					if (rowIndex == colIndex)
					{
						LUP_Track[0].push_back(0);
						LUP_Track[1].push_back(rowIndex);
						LUP_Track[2].push_back(0);
					}
					else if (rowIndex > colIndex)
					{
						int size = LowerIndex_temp[colIndex].size();
						int check = 0;
						for (check = 0; check < size; ++check)
						{
							if (rowIndex == LowerIndex_temp[colIndex][check])
							{
								break;
							}
						}
						if (check >= size) // ����Ԫ
						{
							LowerIndex_temp[colIndex].push_back(rowIndex);

							LUP_Track[0].push_back(-1);
							LUP_Track[1].push_back(colIndex);
							LUP_Track[2].push_back(size);
						}
						else
						{
							LUP_Track[0].push_back(-1);
							LUP_Track[1].push_back(colIndex);
							LUP_Track[2].push_back(check);
						}
					}
					else
					{
						int size = UpperIndex_temp[rowIndex].size();
						int check = 0;
						for (check = 0; check < size; ++check)
						{
							if (colIndex == UpperIndex_temp[rowIndex][check])
							{
								break;
							}
						}
						if (check >= size) // ����Ԫ
						{
							UpperIndex_temp[rowIndex].push_back(colIndex);

							LUP_Track[0].push_back(1);
							LUP_Track[1].push_back(rowIndex);
							LUP_Track[2].push_back(size);
						}
						else
						{
							LUP_Track[0].push_back(1);
							LUP_Track[1].push_back(rowIndex);
							LUP_Track[2].push_back(check);
						}
					}
				}
			}
		}

		////////////////////// Generate IRC, ICFR, LUP////////////////////////

		IRC.resize(nuclide_number_);
		ICFR.resize(2 * nuclide_number_ + 1);

		for (int i = 0; i < nuclide_number_; ++i)
		{
			IRC[i] = i;
		}
		//ICFR[0] = 0;
		ICFR[0] = nuclide_number_;

		//  Low triangle pointers - column compressed
		for (int i = 0; i < nuclide_number_; ++i)
		{
			int size = LowerIndex_temp[i].size();
			for (int j = 0; j<size; ++j)  // 0 -> size
			{
				IRC.push_back(LowerIndex_temp[i][j]);
			}
			ICFR[i + 1] = ICFR[i] + size;
		}

		//  Upper triangle pointers - row compressed
		for (int i = 0; i < nuclide_number_; ++i)
		{
			int size2 = UpperIndex_temp[i].size();
			for (int j = 0; j<size2; ++j)  // 0 -> size
			{
				IRC.push_back(UpperIndex_temp[i][j]);
			}
			ICFR[nuclide_number_ + i + 1] = ICFR[nuclide_number_ + i] + size2;
		}

		int size = LUP_Track[0].size();
		LUP.resize(size);
		for (int i = 0; i<size; ++i)
		{
			//  Check Flag:  0 diagonal   -1 lower triangle  +1 upper triangle
			if (LUP_Track[0][i] == 0)            // diagonal
			{
				LUP[i] = LUP_Track[1][i];
			}
			else if (LUP_Track[0][i] == -1)    // lower triangle 
			{
				LUP[i] = ICFR[LUP_Track[1][i]] + LUP_Track[2][i];
			}
			else    // upper triangle 
			{
				LUP[i] = ICFR[nuclide_number_ + LUP_Track[1][i]] + LUP_Track[2][i];
			}
		}



		return;
	}

	void SymbolLUElimination(const int &nuclide_number_) // ����LU�ֽ�,ͬʱ����ϡ������ce_
	{
		//int nuclide_number_ = cNuclNum;
		vector<vector<int> > UpperIndex_temp(upper_index_), LowerIndex_temp(lower_index_);
		vector<vector<int> > LUP_Track(3, vector<int>(1)); // LUP�б��׷���б�
		for (int i = 0; i < nuclide_number_; ++i)
		{
			int rowsize = LowerIndex_temp[i].size();
			int colsize = UpperIndex_temp[i].size();
			if (i == 1675)
			{
				int jj = 0;
			}
			for (int rr = 0; rr < rowsize; ++rr)
			{
				int rowIndex = LowerIndex_temp[i][rr];
				for (int cc = 0; cc < colsize; ++cc)
				{
					int colIndex = UpperIndex_temp[i][cc];
					if (rowIndex == colIndex)
					{
						LUP_Track[0].push_back(0);
						LUP_Track[1].push_back(rowIndex);
						LUP_Track[2].push_back(0);
					}
					else if (rowIndex > colIndex)
					{
						int size = LowerIndex_temp[colIndex].size();
						int check = 0;
						for (check = 0; check < size; ++check)
						{
							if (rowIndex == LowerIndex_temp[colIndex][check])
							{
								break;
							}
						}
						if (check >= size) // ����Ԫ
						{
							LowerIndex_temp[colIndex].push_back(rowIndex);

							LUP_Track[0].push_back(-1);
							LUP_Track[1].push_back(colIndex);
							LUP_Track[2].push_back(size);
						}
						else
						{
							LUP_Track[0].push_back(-1);
							LUP_Track[1].push_back(colIndex);
							LUP_Track[2].push_back(check);
						}
					}
					else
					{
						int size = UpperIndex_temp[rowIndex].size();
						int check = 0;
						for (check = 0; check < size; ++check)
						{
							if (colIndex == UpperIndex_temp[rowIndex][check])
							{
								break;
							}
						}
						if (check >= size) // ����Ԫ
						{
							UpperIndex_temp[rowIndex].push_back(colIndex);

							LUP_Track[0].push_back(1);
							LUP_Track[1].push_back(rowIndex);
							LUP_Track[2].push_back(size);
						}
						else
						{
							LUP_Track[0].push_back(1);
							LUP_Track[1].push_back(rowIndex);
							LUP_Track[2].push_back(check);
						}
					}
				}
			}
		}

		////////////////////// Generate IRC, ICFR, LUP////////////////////////

		IRC.resize(nuclide_number_);
		ICFR.resize(2 * nuclide_number_ + 1);

		for (int i = 0; i < nuclide_number_; ++i)
		{
			IRC[i] = i;
		}
		//ICFR[0] = 0;
		ICFR[0] = nuclide_number_;

		//  Low triangle pointers - column compressed
		for (int i = 0; i < nuclide_number_; ++i)
		{
			int size = LowerIndex_temp[i].size();
			for (int j = 0; j<size; ++j)  // 0 -> size
			{
				IRC.push_back(LowerIndex_temp[i][j]);
			}
			ICFR[i + 1] = ICFR[i] + size;
		}

		//  Upper triangle pointers - row compressed
		for (int i = 0; i < nuclide_number_; ++i)
		{
			int size2 = UpperIndex_temp[i].size();
			for (int j = 0; j<size2; ++j)  // 0 -> size
			{
				IRC.push_back(UpperIndex_temp[i][j]);
			}
			ICFR[nuclide_number_ + i + 1] = ICFR[nuclide_number_ + i] + size2;
		}

		int size = LUP_Track[0].size();
		LUP.resize(size);
		for (int i = 0; i<size; ++i)
		{
			//  Check Flag:  0 diagonal   -1 lower triangle  +1 upper triangle
			if (LUP_Track[0][i] == 0)            // diagonal
			{
				LUP[i] = LUP_Track[1][i];
			}
			else if (LUP_Track[0][i] == -1)    // lower triangle 
			{
				LUP[i] = ICFR[LUP_Track[1][i]] + LUP_Track[2][i];
			}
			else    // upper triangle 
			{
				LUP[i] = ICFR[nuclide_number_ + LUP_Track[1][i]] + LUP_Track[2][i];
			}
		}



		return;
	}

	void LUElimination(vector< Complex > &VectorB)
	{
		int nuclide_number_ = spmat_dimen_;
		int LUPCount = 1;
		int LUPpos;

		ce_.resize(IRC.size());
		for (int i = 0; i < nuclide_number_; ++i)
		{
			ce_[i] = diagonal_val_[i];
			int RowIB = ICFR[i];
			//int RowIC = ICFR[i + 1];
			int size = lower_index_[i].size();
			for (int j = RowIB; j < RowIB + size; ++j)
			{
				ce_[j] = lower_val_[i][j - RowIB];
			}

			int ColIB = ICFR[nuclide_number_ + i];
			//int ColIC = ICFR[nuclide_number_ + i + 1];
			int size2 = upper_index_[i].size();
			for (int j = ColIB; j < ColIB + size2; ++j)
			{
				ce_[j] = upper_val_[i][j - ColIB];
			}
		}

		vector<Complex > MatrixA(ce_);
		Complex ElimiFactor;

		for (int k = 0; k<nuclide_number_ - 1; ++k)
		{
			////////////////////////// Forward step of Gauss Elimination////////////////////////////
			//if(abs(MatrixA[k])==0)
			//{
			//	Depth.ErrorWarning("Wrong Gauss Elimination����",1);
			//}

			int RowIB = ICFR[k];
			int RowIC = ICFR[k + 1];
			int ColIB = ICFR[nuclide_number_ + k];
			int ColIC = ICFR[nuclide_number_ + k + 1];
			for (int i = RowIB; i < RowIC; ++i)
			{
				int row = IRC[i];     /// Find the row need to be eliminated
				if (row <= k)
				{
					InfoMessage::ErrorMessage("Position: void SpMat::LUElimination; \n Error: the index of lower matrix L has errors.", 1);
				}

				ElimiFactor = MatrixA[i] / MatrixA[k];

				//cout<<MatrixA[i].real()<<'\n'<<MatrixA[i].imag()<<'\n';

				/////// Elimination and filling-in of  Matrix A and its pointers////////
				for (int j = ColIB; j < ColIC; ++j)    /// Eliminate the columns of the given row
				{
					int col = IRC[j];
					if (col <= k)
					{
						InfoMessage::ErrorMessage("Position: void SpMat::LUElimination; \n Error: the index of upper matrix U has errors.", 1);
					}

					LUPpos = LUP[LUPCount];

					//MatrixA[LUPpos] = MatrixA[LUPpos] - MatrixA[j] * ElimiFactor;
					MatrixA[LUPpos] -= MatrixA[j] * ElimiFactor;
					++LUPCount;

				} // end for(j=0;j<size2;++j)
				  /////// Elimination and filling-in of  vector b////////
				if (VectorB[k] != 0.0)
				{
					//VectorB[row] = VectorB[row] - VectorB[k] * ElimiFactor;
					VectorB[row] -= VectorB[k] * ElimiFactor;
				}
			}// end for(i=0;i<size1;++i)
		}// endffor(k=1;k<=nuclide_number_-1;++k)


		 ////////////////////////// Backward step of Gauss Elimination////////////////////////////
		for (int i = nuclide_number_ - 1; i >= 0; --i)
		{
			int ColIB = ICFR[nuclide_number_ + i];
			int ColIC = ICFR[nuclide_number_ + i + 1];
			//cout<<k<<"   "<<size1<<"   "<<size2<<'\n';
			for (int j = ColIB; j < ColIC; ++j)
			{
				int col = IRC[j];
				if (col <= i)
				{
					InfoMessage::ErrorMessage("Position: void SpMat::LUElimination; \n Error: the index of lower matrix L has errors.",1);
				}

				VectorB[i] -= MatrixA[j] * VectorB[col];
			}
			VectorB[i] = VectorB[i] / MatrixA[i];
		}

		return;
	}

	void LUElimination(vector< Complex > &VectorB, int &nuclide_number_)
	{
		//,4int nuclide_number_ = cNuclNum;
		int LUPCount = 1;
		int LUPpos;

		ce_.resize(IRC.size());
		for (int i = 0; i < nuclide_number_; ++i)
		{
			ce_[i] = diagonal_val_[i];
			int RowIB = ICFR[i];
			//int RowIC = ICFR[i + 1];
			int size = lower_index_[i].size();
			for (int j = RowIB; j < RowIB + size; ++j)
			{
				ce_[j] = lower_val_[i][j - RowIB];
			}

			int ColIB = ICFR[nuclide_number_ + i];
			//int ColIC = ICFR[nuclide_number_ + i + 1];
			int size2 = upper_index_[i].size();
			for (int j = ColIB; j < ColIB + size2; ++j)
			{
				ce_[j] = upper_val_[i][j - ColIB];
			}
		}

		vector<Complex > MatrixA(ce_);
		Complex ElimiFactor;

		for (int k = 0; k<nuclide_number_ - 1; ++k)
		{
			////////////////////////// Forward step of Gauss Elimination////////////////////////////
			//if(abs(MatrixA[k])==0)
			//{
			//	Depth.ErrorWarning("Wrong Gauss Elimination����",1);
			//}

			int RowIB = ICFR[k];
			int RowIC = ICFR[k + 1];
			int ColIB = ICFR[nuclide_number_ + k];
			int ColIC = ICFR[nuclide_number_ + k + 1];
			for (int i = RowIB; i < RowIC; ++i)
			{
				int row = IRC[i];     /// Find the row need to be eliminated
				if (row <= k)
				{
					InfoMessage::ErrorMessage("Position: void SpMat::LUElimination; \n Error: the index of lower matrix L has errors.", 1);
				}

				ElimiFactor = MatrixA[i] / MatrixA[k];

				//cout<<MatrixA[i].real()<<'\n'<<MatrixA[i].imag()<<'\n';

				/////// Elimination and filling-in of  Matrix A and its pointers////////
				for (int j = ColIB; j < ColIC; ++j)    /// Eliminate the columns of the given row
				{
					int col = IRC[j];
					if (col <= k)
					{
						InfoMessage::ErrorMessage("Position: void SpMat::LUElimination; \n Error: the index of upper matrix U has errors.", 1);
					}

					LUPpos = LUP[LUPCount];

					MatrixA[LUPpos] -= MatrixA[j] * ElimiFactor;
					++LUPCount;

				} // end for(j=0;j<size2;++j)
				  /////// Elimination and filling-in of  vector b////////
				if (VectorB[k] != 0.0)
				{
					VectorB[row] -= VectorB[k] * ElimiFactor;
				}
			}// end for(i=0;i<size1;++i)
		}// endffor(k=1;k<=nuclide_number_-1;++k)


		 ////////////////////////// Backward step of Gauss Elimination////////////////////////////
		for (int i = nuclide_number_ - 1; i >= 0; --i)
		{
			int ColIB = ICFR[nuclide_number_ + i];
			int ColIC = ICFR[nuclide_number_ + i + 1];
			//cout<<k<<"   "<<size1<<"   "<<size2<<'\n';
			for (int j = ColIB; j < ColIC; ++j)
			{
				int col = IRC[j];
				if (col <= i)
				{
					InfoMessage::ErrorMessage("Position: void SpMat::LUElimination; \n Error: the index of lower matrix L has errors.", 1);
				}

				VectorB[i] -= MatrixA[j] * VectorB[col];
			}
			VectorB[i] = VectorB[i] / MatrixA[i];
		}

		return;
	}

	void LUEliminationForCram(const Complex &theta, const double &time, vector< Complex > &VectorB)
	{
		int nuclide_number_ = spmat_dimen_;
		int LUPCount = 1;
		int LUPpos;

		//////////   ʹ��cram�����ľ���Ԥ���� ////////////
		ce_.resize(IRC.size());
		for (int i = 0; i < nuclide_number_; ++i)
		{
			ce_[i] = diagonal_val_[i] * time - theta;
			int RowIB = ICFR[i];
			//int RowIC = ICFR[i + 1];
			int size = lower_index_[i].size();
			for (int j = RowIB; j < RowIB + size; ++j)
			{
				ce_[j] = lower_val_[i][j - RowIB] * time;
			}

			int ColIB = ICFR[nuclide_number_ + i];
			//int ColIC = ICFR[nuclide_number_ + i + 1];
			int size2 = upper_index_[i].size();
			for (int j = ColIB; j < ColIB + size2; ++j)
			{
				ce_[j] = upper_val_[i][j - ColIB] * time;
			}
		}
		////////////////////////////////////////////////

		vector<Complex > MatrixA(ce_);
		Complex ElimiFactor;

		for (int k = 0; k<nuclide_number_ - 1; ++k)
		{
			////////////////////////// Forward step of Gauss Elimination////////////////////////////
			//if(abs(MatrixA[k])==0)
			//{
			//	Depth.ErrorWarning("Wrong Gauss Elimination����",1);
			//}

			int RowIB = ICFR[k];
			int RowIC = ICFR[k + 1];
			int ColIB = ICFR[nuclide_number_ + k];
			int ColIC = ICFR[nuclide_number_ + k + 1];
			for (int i = RowIB; i < RowIC; ++i)
			{
				int row = IRC[i];     /// Find the row need to be eliminated
				if (row <= k)
				{
					InfoMessage::ErrorMessage("Position: void SpMat::LUEliminationForCram; \n Error: the index of lower matrix L has errors.", 1);
				}

				ElimiFactor = MatrixA[i] / MatrixA[k];

				//cout<<MatrixA[i].real()<<'\n'<<MatrixA[i].imag()<<'\n';

				/////// Elimination and filling-in of  Matrix A and its pointers////////
				for (int j = ColIB; j < ColIC; ++j)    /// Eliminate the columns of the given row
				{
					int col = IRC[j];
					if (col <= k)
					{
						InfoMessage::ErrorMessage("Position: void SpMat::LUEliminationForCram; \n Error: the index of upper matrix U has errors.", 1);
					}

					LUPpos = LUP[LUPCount];

					MatrixA[LUPpos] -= MatrixA[j] * ElimiFactor;
					++LUPCount;

				} // end for(j=0;j<size2;++j)
				  /////// Elimination and filling-in of  vector b////////
				if (VectorB[k] != 0.0)
				{
					VectorB[row] -= VectorB[k] * ElimiFactor;
				}
			}// end for(i=0;i<size1;++i)
		}// endffor(k=1;k<=nuclide_number_-1;++k)


		 ////////////////////////// Backward step of Gauss Elimination////////////////////////////
		for (int i = nuclide_number_ - 1; i >= 0; --i)
		{
			int ColIB = ICFR[nuclide_number_ + i];
			int ColIC = ICFR[nuclide_number_ + i + 1];
			//cout<<k<<"   "<<size1<<"   "<<size2<<'\n';
			for (int j = ColIB; j < ColIC; ++j)
			{
				int col = IRC[j];
				if (col <= i)
				{
					InfoMessage::ErrorMessage("Position: void SpMat::LUEliminationForCram; \n Error: the index of lower matrix L has errors." ,1);
				}

				VectorB[i] -= MatrixA[j] * VectorB[col];
			}
			VectorB[i] = VectorB[i] / MatrixA[i];
		}

		return;
	}

	void LUEliminationForCram(const Complex &theta, const double &time, vector<Complex > &VectorB, const int &nuclide_number_)
	{
		//int nuclide_number_ = cNuclNum;
		int LUPCount = 1;
		int LUPpos;

		//////////   ʹ��cram�����ľ���Ԥ���� ////////////
		ce_.resize(IRC.size());
		for (int i = 0; i < nuclide_number_; ++i)
		{
			ce_[i] = diagonal_val_[i] * time - theta;
			int RowIB = ICFR[i];
			//int RowIC = ICFR[i + 1];
			int size = lower_index_[i].size();
			for (int j = RowIB; j < RowIB + size; ++j)
			{
				ce_[j] = lower_val_[i][j - RowIB] * time;
			}

			int ColIB = ICFR[nuclide_number_ + i];
			//int ColIC = ICFR[nuclide_number_ + i + 1];
			int size2 = upper_index_[i].size();
			for (int j = ColIB; j < ColIB + size2; ++j)
			{
				ce_[j] = upper_val_[i][j - ColIB] * time;
			}
		}
		////////////////////////////////////////////////

		vector<Complex > MatrixA(ce_);
		Complex ElimiFactor;

		for (int k = 0; k<nuclide_number_ - 1; ++k)
		{

			int RowIB = ICFR[k];
			int RowIC = ICFR[k + 1];
			int ColIB = ICFR[nuclide_number_ + k];
			int ColIC = ICFR[nuclide_number_ + k + 1];
			for (int i = RowIB; i < RowIC; ++i)
			{
				int row = IRC[i];     
				if (row <= k)
				{
					InfoMessage::ErrorMessage("Position: void SpMat::LUEliminationForCram; \n Error: the index of lower matrix L has errors.", 1);
				}

				ElimiFactor = MatrixA[i] / MatrixA[k];


				/////// Elimination and filling-in of  Matrix A and its pointers////////
				for (int j = ColIB; j < ColIC; ++j)    /// Eliminate the columns of the given row
				{
					int col = IRC[j];
					if (col <= k)
					{
						InfoMessage::ErrorMessage("Position: void SpMat::LUEliminationForCram; \n Error: the index of upper matrix U has errors.", 1);
					}

					LUPpos = LUP[LUPCount];

					MatrixA[LUPpos] -= MatrixA[j] * ElimiFactor;
					++LUPCount;

				} 
				  /////// Elimination and filling-in of  vector b////////
				if (VectorB[k] != 0.0)
				{
					VectorB[row] -= VectorB[k] * ElimiFactor;
				}
			}
		}


		 ////////////////////////// Backward step of Gauss Elimination////////////////////////////
		for (int i = nuclide_number_ - 1; i >= 0; --i)
		{
			int ColIB = ICFR[nuclide_number_ + i];
			int ColIC = ICFR[nuclide_number_ + i + 1];
			for (int j = ColIB; j < ColIC; ++j)
			{
				int col = IRC[j];
				if (col <= i)
				{
					InfoMessage::ErrorMessage("Position: void SpMat::LUEliminationForCram; \n Error: the index of lower matrix L has errors.", 1);
				}

				VectorB[i] -= MatrixA[j] * VectorB[col];
			}
			VectorB[i] = VectorB[i] / MatrixA[i];
		}

		return;
	}
};

#endif
