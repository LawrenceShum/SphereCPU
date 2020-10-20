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

//���A����
void LinearSolver::set_value_A(int row, int col, double value)
{
	A(row, col) = value;
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
	cout << A << endl;
}

void LinearSolver::output_eigenvalue()
{
	EigenSolver<MatrixXd> es(A);
	cout << "����A������ֵΪ" << es.eigenvalues() << endl;
}

void LinearSolver::Jakobi()
{

}

void LinearSolver::PCG()
{
}