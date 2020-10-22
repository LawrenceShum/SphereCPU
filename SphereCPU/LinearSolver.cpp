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

//���A����
void LinearSolver::set_value_A(int row, int col, double value)
{
	A(row, col) = value;
}

//����b����
void LinearSolver::set_value_b(int row, double value)
{
	b(row) = value;
}

//��ʼ��A���󣬳�ʼΪ�����
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

//���A����
Eigen::MatrixXd LinearSolver::get_matrix_A()
{
	return A;
}

//���b����
Eigen::VectorXd LinearSolver::get_vector_b()
{
	return b;
}

void LinearSolver::output_A()
{
	cout << "����AΪ��" << endl;
	cout << A << endl;
}

void LinearSolver::output_b()
{
	cout << "����bΪ��" << endl;
	cout << b << endl;
}

void LinearSolver::output_eigenvalue()
{
	EigenSolver<MatrixXd> es(A);
	cout << "����A������ֵΪ" << endl << es.eigenvalues() << endl;
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
	cout << "househoulderQR����ʹ��ʱ��Ϊ��" << 1000 * (clock() - time) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << x << endl;
}

void LinearSolver::colPivHouseholderQR()
{
	clock_t time = clock();
	x = A.colPivHouseholderQr().solve(b);
	cout << "colPivHousehoulderQR����ʹ��ʱ��Ϊ��" << 1000 * (clock() - time) / (double)CLOCKS_PER_SEC << "ms" << endl;
	cout << x << endl;
}

void LinearSolver::fullPivHouseholderQR()
{
	clock_t time = clock();
	x = A.fullPivHouseholderQr().solve(b);
	cout << "fullPivHousehoulderQR����ʹ��ʱ��Ϊ��" << 1000 * (clock() - time) / (double)CLOCKS_PER_SEC << "ms" << endl;
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