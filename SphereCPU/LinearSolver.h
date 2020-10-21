#pragma once
#include "TotalInclude.h"
#include <Eigen/core>
#include <Eigen/Dense>

class LinearSolver
{
private:
	//�ߴ�
	int size;
	//����A��Ax=b
	Eigen::MatrixXd	A;
	//����b
	Eigen::VectorXd b;
	//����x
	Eigen::VectorXd x;
	//��ʼ��A����
	void initialize_A();
	//��ʼ������b
	void initialize_b();
	//��ʼ������x
	void initialize_x();

	//PCG�ⷨ
	void PCG();
	//Gauss-Seidel�ⷨ
	void GS();
	//jakobi�ⷨ
	void Jakobi();
	//QR�ֽ�
	void householderQR();
	//colpivhouseholderQR
	void colPivHouseholderQR();
	//full piv householderQR
	void fullPivHouseholderQR();

	//��������
	int iteration;

public:
	LinearSolver(int,int);

	//��ȡA����
	Eigen::MatrixXd get_matrix_A();
	//��ȡ����b
	Eigen::VectorXd get_vector_b();
	//�趨A�����ڳ�Ա��ֵ
	void set_value_A(int,int,double);
	//�趨b������ֵ
	void set_value_b(int, double);
	//���A����
	void output_A();
	//���b����
	void output_b();
	//���A��������ֵ
	void output_eigenvalue();
	//�����Է�����
	void solveLinear(int);

	~LinearSolver();
};