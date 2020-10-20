#include "LinearSolver.h"
#include "TotalInclude.h"
#include <iostream>

using namespace std;
using namespace Eigen;

LinearSolver::LinearSolver(int size, int iteration)
	:size(size), iteration(iteration)
{
	initialize_A();
	initialize_b();
	initialize_x();
}

LinearSolver::~LinearSolver()
{

}

//获得A矩阵
void LinearSolver::set_value_A(int row, int col, double value)
{
	A(row, col) = value;
}

//初始化A矩阵，初始为零矩阵
void LinearSolver::initialize_A()
{
	A = Eigen::MatrixXd::Zero(size,size);
}

void LinearSolver::initialize_b()
{
	b = Eigen::VectorXd::Zero(size);
}

void LinearSolver::initialize_x()
{
	x = Eigen::VectorXd::Zero(size);
}

//获得A矩阵
Eigen::MatrixXd LinearSolver::get_matrix_A()
{
	return A;
}

//获得b向量
Eigen::VectorXd LinearSolver::get_vector_b()
{
	return b;
}

void LinearSolver::output_A()
{
	cout << A << endl;
}

void LinearSolver::output_eigenvalue()
{
	EigenSolver<MatrixXd> es(A);
	cout << "矩阵A的特征值为" << es.eigenvalues() << endl;
}

void LinearSolver::Jakobi()
{

}

void LinearSolver::PCG()
{
}