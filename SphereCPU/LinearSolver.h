#pragma once
#include "TotalInclude.h"
#include <Eigen/core>
#include <Eigen/Dense>

class LinearSolver
{
private:
	//尺寸
	int size;
	//矩阵A，Ax=b
	Eigen::MatrixXd	A;
	//向量b
	Eigen::VectorXd b;
	//向量x
	Eigen::VectorXd x;
	//初始化A矩阵
	void initialize_A();
	//初始化向量b
	void initialize_b();
	//初始化向量x
	void initialize_x();

	//PCG解法
	void PCG();
	//Gauss-Seidel解法
	void GS();
	//jakobi解法
	void Jakobi();
	//QR分解
	void householderQR();
	//colpivhouseholderQR
	void colPivHouseholderQR();
	//full piv householderQR
	void fullPivHouseholderQR();

	//迭代次数
	int iteration;

public:
	LinearSolver(int,int);

	//获取A矩阵
	Eigen::MatrixXd get_matrix_A();
	//获取向量b
	Eigen::VectorXd get_vector_b();
	//设定A矩阵内成员的值
	void set_value_A(int,int,double);
	//设定b向量的值
	void set_value_b(int, double);
	//输出A矩阵
	void output_A();
	//输出b向量
	void output_b();
	//输出A矩阵特征值
	void output_eigenvalue();
	//解线性方程组
	void solveLinear(int);

	~LinearSolver();
};