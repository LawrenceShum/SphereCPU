#pragma once
#include "TotalInclude.h"
#include <Eigen/core>
#include <Eigen/Dense>

class LinearSolver
{
private:
	//��
	int row;
	//��
	int column;
	//����A��Ax=b
	Eigen::MatrixXd	A;
	//����b
	Eigen::VectorXd b;

	//PCG�ⷨ
	void PCG();
	//Gauss-Seidel�ⷨ
	void GS();
	//jakobi�ⷨ
	void Jakobi();
	//�����ݷ���A������
	void fillA(double*);

public:
	LinearSolver(int,int);

	//��ȡA����
	Eigen::MatrixXd get_matrix_A();
	//��ȡ����b
	Eigen::VectorXd get_vector_b();

	~LinearSolver();
};