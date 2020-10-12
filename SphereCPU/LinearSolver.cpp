#include "LinearSolver.h"
#include "TotalInclude.h"

using namespace std;
using namespace Eigen;

LinearSolver::LinearSolver(int size, int iteration)
	:size(size), iteration(iteration)
{
	initialize_A();
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