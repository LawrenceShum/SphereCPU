#include "LinearSolver.h"
#include "TotalInclude.h"
#include <iostream>
#include <ctime>

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

//设置b向量
void LinearSolver::set_value_b(int row, double value)
{
	b(row) = value;
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
	cout << "矩阵A为：" << endl;
	cout << A << endl;
}

void LinearSolver::output_b()
{
	cout << "向量b为：" << endl;
	cout << b << endl;
}

void LinearSolver::output_eigenvalue()
{
	EigenSolver<MatrixXd> es(A);
	cout << "矩阵A的特征值为" << endl << es.eigenvalues() << endl;
}

double LinearSolver::get_value_b(int n)
{
	return b(n);
}

void LinearSolver::Jakobi()
{

}

void LinearSolver::PCG()
{

}

void LinearSolver::householderQR()
{
	clock_t time = clock();
	x = A.householderQr().solve(b);
	cout << "househoulderQR方法使用时间为：" << 1000 * (clock() - time) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << x << endl;
}

void LinearSolver::colPivHouseholderQR()
{
	clock_t time = clock();
	x = A.colPivHouseholderQr().solve(b);
	cout << "colPivHousehoulderQR方法使用时间为：" << 1000 * (clock() - time) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << x << endl;
}

void LinearSolver::fullPivHouseholderQR()
{
	clock_t time = clock();
	x = A.fullPivHouseholderQr().solve(b);
	cout << "fullPivHousehoulderQR方法使用时间为：" << 1000 * (clock() - time) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << x << endl;
}

void LinearSolver::solveLinear(int c)
{
	if (!c)
		householderQR();
	else if (c == 1)
		colPivHouseholderQR();
	else if (c == 2)
		fullPivHouseholderQR();
}