#pragma once
#include "TotalInclude.h"
#include <Eigen/core>
#include <Eigen/Dense>

class LinearSolver
{
private:
	//行
	int row;
	//列
	int column;
	//矩阵A，Ax=b
	Eigen::MatrixXd	A;
	//向量b
	Eigen::VectorXd b;

	//PCG解法
	void PCG();
	//Gauss-Seidel解法
	void GS();
	//jakobi解法
	void Jakobi();
	//将数据放入A矩阵中
	void fillA(double*);

public:
	LinearSolver(int,int);

	//获取A矩阵
	Eigen::MatrixXd get_matrix_A();
	//获取向量b
	Eigen::VectorXd get_vector_b();

	~LinearSolver();
};