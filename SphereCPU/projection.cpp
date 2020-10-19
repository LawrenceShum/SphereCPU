#include "SphereSolver.h"
#include "TotalInclude.h"
#include "LinearSolver.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <cmath>
#include <Eigen/IterativeLinearSolvers>
#include <unsupported/Eigen/IterativeSolvers>
//#include <unsupported/Eigen/FFT>
# include "mkl_service.h"
/* Include Poisson Library header files */
# include "mkl_poisson.h"
# include "mkl_dfti.h"

using namespace std;
using namespace Eigen;

int get_column_num(int i, int j, int pitch)
{
	return (j - 1)*pitch + i;
}

void SphereSolver::projection(int a)
{
	float invRadius = 1.0 / radius;
	float invDensity = 1.0 / density;

	int iteration = 20;

	//做projection的时候，极点被考虑成为一个点而不是之前的展开
	int size = n_phi*(n_theta - 1) + 2;

	//构建一个线性求解器
	//线性系统为Ax = b
	LinearSolver linear(size, iteration);

	//构建A矩阵
	//北极、南极两点的P放在矩阵A中每一行的最后两个
	//??problem??
	//还差两个方程式，北极和南极的projection
	//trick1:极点处的压力为极点周围压力的平均值
	int row = 0;
	//先考虑非极点情况
	for (int j = 1; j < n_theta; j++)
		{
			double scale1 = 0.5*dt*invDensity*invRadius*invGridLen*invGridLen;
			for (int i = 0; i < n_phi; i++)
			{
				int left_i, left_j,
					right_i, right_j,
					up_i, up_j,
					down_i, down_j;
				//判断左右脚标
				if (i == 0)
				{
					left_i = n_phi - 2;
					right_i = 2;
					left_j = j;
					right_j = j;
				}
				else if (i == 1)
				{
					left_i = n_phi - 1;
					right_i = 3;
					left_j = j;
					right_j = j;
				}
				else if (i == n_phi - 2)
				{
					left_i = n_phi - 4;
					right_i = 0;
					left_j = j;
					right_j = j;
				}
				else if (i == n_phi - 1)
				{
					left_i = n_phi - 3;
					right_i = 1;
					left_j = j;
					right_j = j;
				}
				else 
				{
					left_i = i - 2;
					right_i = i + 2;
					left_j = j;
					right_j = j;
				}

				//判断上下脚标
				if (j == 1)
				{
					up_i = (i + n_phi / 2) % n_phi;
					down_i = i;
					up_j = 1;
					down_j = 3;
				}
				else if (j == 2)
				{
					up_i = 0;
					down_i = i;
					up_j = n_theta;
					down_j = 4;
				}
				else if (j == n_theta - 1)
				{
					up_i = i;
					down_i = (i + n_phi / 2) % n_phi;
					up_j = j - 2;
					down_j = n_theta - 1;
				}
				else if (j == n_theta - 2)
				{
					up_i = i;
					down_i = 1;
					up_j = j - 2;
					down_j = n_theta;
				}
				else
				{
					up_i = i;
					down_i = i;
					up_j = j - 2;
					down_j = j + 2;
				}
				//将自身系数设为4
				linear.set_value_A(row, get_column_num(i, j, n_phi), 4);
				//将临近上下左右系数设为-1
				linear.set_value_A(row, get_column_num(left_i, left_j, n_phi), -1);
				linear.set_value_A(row, get_column_num(right_i, right_j, n_phi), -1);
				linear.set_value_A(row, get_column_num(up_i, up_j, n_phi), -1);
				linear.set_value_A(row, get_column_num(down_i, down_j, n_phi), -1);

				row++;
			}
		}
	//考虑北极南极点的情况
	//北极点
	for (int i = 0; i < n_phi; i++)
	{
		linear.set_value_A(row, get_column_num(i, 1, n_phi), 1);
		linear.set_value_A(row, get_column_num(0, n_theta, n_phi), -1);
	}
	//南极点
	row++;
	for (int i = 0; i < n_phi; i++)
	{
		linear.set_value_A(row, get_column_num(i, n_theta-1, n_phi), 1);
		linear.set_value_A(row, get_column_num(1, n_theta, n_phi), -1);
	}
	linear.output_A();
	//构建rhs，b


	//gradP储存压力P的梯度
	float* gradP = new float[n_phi*(n_theta + 1)];
	float scale = dt * invRadius * invDensity;

	//使用MKL计算projection
	if (!a)
	{
		//u为phi方向，v为theta方向
		float* u = vel_phi_this;
		float* v = vel_theta_this;
		//MKLj解泊松方程
		MKLSpherical();
		//更新phi方向的速度
		for (int y = 1; y < n_theta - 1; y++)
		{
			for (int x = 0; x < n_phi; x++)
			{
				float theta = y*gridLen;
				float phi = x*gridLen;
				float invSin = 1.0 / sin(theta);
				float factorPhi = -dt * invGridLen * invRadius * invSin * invDensity * 0.5;

				int gridLeftX = (x == 0 ? n_phi - 1 : x - 1);
				int gridRightX = (x == n_phi - 1 ? 0 : x + 1);

				float pressureGrad = presure_this[gridRightX + y*n_phi] - presure_this[gridLeftX + y*n_phi];
				float deltauPhi = factorPhi * pressureGrad;
				vel_phi_next[x + y*n_phi] = u[x + y*n_phi] + deltauPhi;
			}
		}
		//更新theta方向的速度
		for (int y = 1; y < n_theta - 1; y++)
		{
			for (int x = 0; x < n_phi; x++)
			{
				float theta = y*gridLen;
				float phi = x*gridLen;
				float factorTheta = -dt * invGridLen * invRadius * invDensity * 0.5;

				int gridLeftX = (x == 0 ? n_phi - 1 : x - 1);
				int gridRightX = (x == n_phi - 1 ? 0 : x + 1);

				float pressureGrad = presure_this[gridRightX + y*n_phi] - presure_this[gridLeftX + y*n_phi];
				float deltauTheta = factorTheta * pressureGrad;
				vel_theta_next[x + y*n_phi] = v[x + y*n_phi] + deltauTheta;
			}
		}
		//解决极点处的projection
		solvePolarProjection();
	}

	//使用PCG解决projection
	else if (a == 1)
	{
		
	}

	//使用GS解决projection
	else if (a == 2)
	{
		
	}

	//使用jacob解决projection
	else if (a == 3)
	{
		
	}
	else
		cout << "Something goes wrong in the projection step" << endl;

	swap_vel_phi();
	swap_vel_theta();
}


void SphereSolver::solvePolarProjection()
{
	float* v = vel_theta_this;
	float* p = presure_this;

	float invRadius = 1.0 / radius;
	float invDensity = 1.0 / density;

	//北极点处
	for (int k = 0; k < n_phi; k++)
	{
		float factor = -dt * invGridLen * invRadius * invDensity * 0.5;
		float u_n = v[k];
		float p_n = p[k];
		float p_down = p[k + 1 * n_phi];

		float deltaP = p_down - p_n;
		float deltaTheta = factor*deltaP;
		vel_theta_next[k] = u_n + deltaTheta;
	}
	//南极点处
	for (int k = 0; k < n_phi; k++)
	{
		float factor = -dt * invGridLen * invRadius * invDensity * 0.5;
		float u_n = v[k + n_phi*n_theta];
		float p_n = p[k + (n_phi-1)*n_theta];
		float p_down = p[k + n_phi*n_theta];

		float deltaP = p_down - p_n;
		float deltaTheta = factor*deltaP;
		vel_theta_next[k + n_phi*n_theta] = u_n + deltaTheta;
	}
	//共轭计算
	for (int k = 0; k < n_phi; k++)
	{
		int num = (k + static_cast<int>(n_phi / 4)) % n_phi;
		vel_phi_next[k] = vel_theta_next[num];
		vel_phi_next[k + n_phi*n_theta] = vel_theta_next[num + n_phi*n_theta];
	}
}

/*******************************************************************************
* Copyright 2006-2018 Intel Corporation.
*
* This software and the related documents are Intel copyrighted  materials,  and
* your use of  them is  governed by the  express license  under which  they were
* provided to you (License).  Unless the License provides otherwise, you may not
* use, modify, copy, publish, distribute,  disclose or transmit this software or
* the related documents without Intel's prior written permission.
*
* This software and the related documents  are provided as  is,  with no express
* or implied  warranties,  other  than those  that are  expressly stated  in the
* License.
*******************************************************************************/

/*
*  Content:
*  C double precision example of solving Helmholtz problem on a whole sphere
*  using Intel(R) MKL Poisson Library
*
*******************************************************************************/
void SphereSolver::MKLSpherical()
{
	float* uPhi = vel_phi_this;
	float* vTheta = vel_theta_this;
	float* p = presure_this;

	MKL_INT np = n_phi, nt = n_theta - 1;
	double invRad = 1.0 / radius;
	double pi = 3.14159265358979324;

	MKL_INT ip, it, i, stat;
	MKL_INT ipar[128];
	double ap, bp, at, bt, lp, lt, hp, ht, theta_i, theta_up, theta_down, invSin, sin_up, sin_down, ct, c1;
	double *dpar = NULL, *f = NULL, *u = NULL;
	double q;
	DFTI_DESCRIPTOR_HANDLE handle_s = 0;
	DFTI_DESCRIPTOR_HANDLE handle_c = 0;
	MKL_INT mem_error, error;

	error = 0;
	/* memory allocation */
	mem_error = 1;
	dpar = (double*)mkl_malloc((5 * np / 2 + nt + 10) * sizeof(double), 64);
	if (dpar == NULL) goto end;
	f = (double*)mkl_malloc((np + 1)*(nt + 1) * sizeof(double), 64);
	if (f == NULL) goto end;
	/* memory allocated correctly */
	mem_error = 0;

	/* Defining the rectangular domain on a sphere 0<p<2*pi, 0<t<pi
	for Helmholtz Solver on a sphere
	Poisson Library will automatically detect that this problem is on a whole sphere! */
	ap = 0.0E0;
	bp = 2 * pi;
	at = 0.0E0;
	bt = pi;

	/* Setting the coefficient q to 0.0E0 for Poisson problem
	If you like to solve Helmholtz problem, please set q to 1.0E0 */
	q = 0.0E0;

	/* Computing the mesh size hp in phi-direction */
	lp = bp - ap;
	hp = lp / np;
	/* Computing the mesh size ht in theta-direction */
	lt = bt - at;
	ht = lt / nt;

	// interior of grid
	for (it = 1; it < nt; it++)
	{
		for (ip = 1; ip < np - 1; ip++)
		{
			theta_i = ht*it;
			theta_up = ht*(it + 1);
			theta_down = ht*(it - 1);
			invSin = 1.0 / sin(theta_i);
			sin_up = sin(theta_up);
			sin_down = sin(theta_down);
			double factor = invRad * 0.5 * invGridLen * invSin * density;
			float u_left = uPhi[ip + it*np];
			//float u_left = uPhi->getValueAt(ip - 1, it);
			float u_right = uPhi[ip + (it + 1)*np];
			//float u_right = uPhi->getValueAt(ip + 1, it);
			float u_up = vTheta[ip + (it + 1)*np];
			//float u_up = vTheta->getValueAt(ip, it + 1);
			float u_down = vTheta[ip + (it - 1)*np];
			//float u_down = vTheta->getValueAt(ip, it - 1);
			f[ip + it*(np + 1)] = factor * (u_right - u_left + u_up * sin_up - u_down * sin_down);
		}
	}
	// second to last column
	for (it = 1; it < nt; it++)
	{
		theta_i = ht*it;
		theta_up = ht*(it + 1);
		theta_down = ht*(it - 1);
		invSin = 1.0 / sin(theta_i);
		sin_up = sin(theta_up);
		sin_down = sin(theta_down);
		double factor = invRad * 0.5 * invGridLen * invSin * density;
		float u_left = uPhi[(np - 2) + it*np];
		//float u_left = uPhi->getValueAt(np - 2, it);
		float u_right = uPhi[0 + it*np];
		//float u_right = uPhi->getValueAt(0, it);
		float u_up = vTheta[(np - 1) + (it + 1)*np];
		//float u_up = vTheta->getValueAt(np - 1, it + 1);
		float u_down = vTheta[(np - 1) + (it - 1)*np];
		//float u_down = vTheta->getValueAt(np - 1, it - 1);
		f[ip + it*(np + 1)] = factor * (u_right - u_left + u_up * sin_up - u_down * sin_down);
	}

	// phi = 0, 2Pi seam
	for (it = 1; it < nt; it++) {
		theta_i = ht*it;
		theta_up = ht*(it + 1);
		theta_down = ht*(it - 1);
		invSin = 1.0 / sin(theta_i);
		sin_up = sin(theta_up);
		sin_down = sin(theta_down);
		double factor = invRad * 0.5 * invGridLen * invSin;
		float u_left = uPhi[(np - 1) + it*np];
		//float u_left = uPhi->getValueAt(np - 1, it);
		float u_right = uPhi[1 + it*np];
		//float u_right = uPhi->getValueAt(1, it);
		float u_up = vTheta[0 + (it + 1)*np];
		//float u_up = vTheta->getValueAt(0, it + 1);
		float u_down = vTheta[0 + (it - 1)*np];
		//float u_down = vTheta->getValueAt(0, it - 1);
		float seamVal = factor * (u_right - u_left + u_up * sin_up - u_down * sin_down);
		f[it * (np + 1)] = seamVal;
		f[np + it * (np + 1)] = seamVal;
	}

	// poles
	for (ip = 0; ip <= np; ip++) {
		//gridShift = (ip + np / 2) % np;
		//fReal factor = invRad * 0.5 * invGridLen;
		//fReal u_down = vTheta->getValueAt(ip, 1);
		//fReal u_up = vTheta->getValueAt(gridShift, 1);
		//fReal NP_val = factor * (u_up - u_down);
		//f[ip] = NP_val;
		f[ip] = 0.0;
		//u_down = vTheta->getValueAt(ip, np - 2);
		//u_up = vTheta->getValueAt(gridShift, np - 2);
		//fReal SP_val = factor * (u_up - u_down);
		//f[ip + nt * (np + 1)] = SP_val;
		f[ip + nt * (np + 1)] = 0.0;
	}

	/* Initializing ipar array to make it free from garbage */
	for (i = 0; i<128; i++)
	{
		ipar[i] = 0;
	}

	/* Initializing simple data structures of Poisson Library
	for Poisson Solver on a sphere
	As we are looking for the solution on a whole sphere, this is a PERIODIC problem
	Therefore, the routines ending with "_p" are used to find the solution */
	d_init_sph_p(&ap, &bp, &at, &bt, &np, &nt, &q, ipar, dpar, &stat);
	if (stat < 0) {
		error = 1;
		goto end;
	}

	/* Initializing complex data structures of Poisson Library
	for Poisson Solver on a sphere
	NOTE: Right-hand side f may be altered after the Commit step. If you want
	to keep it, you should save it in another memory location! */
	d_commit_sph_p(f, &handle_s, &handle_c, ipar, dpar, &stat);
	if (stat < 0) {
		error = 1;
		goto end;
	}
	/* Computing the approximate solution of Poisson problem on a whole sphere */
	// note: solution is stored in f
	d_sph_p(f, &handle_s, &handle_c, ipar, dpar, &stat);
	if (stat < 0) {
		error = 1;
		goto end;
	}
	/* Cleaning the memory used by handle_s and handle_c */
	free_sph_p(&handle_s, &handle_c, ipar, &stat);
	if (stat < 0) {
		error = 1;
		goto end;
	}
	/* Now we can use handle_s and handle_c to solve another Poisson problem */
	/* after a proper initialization */

	/* Printing the results */
	//printf("The number of mesh intervals in phi-direction is np=%d\n", np);
	//printf("The number of mesh intervals in theta-direction is nt=%d\n\n", nt);

	for (it = 0; it <= nt; it++)
	{
		for (ip = 0; ip < np; ip++)
		{
			//printf("%10.3f", f[ip + it*(np + 1)]);
			double F = f[ip + it*(np + 1)];
			this->set_presure(ip, it, F);
			//p->writeValueTo(ip, it, f[ip + it*(np + 1)]);
		}
	}
	// TODO: check divergence free here

end:
	/* Free Intel(R) MKL memory if any was allocated */
	mkl_free(dpar);
	mkl_free(f);
	MKL_Free_Buffers();
	/* Failure message to print if something went wrong */
	if (mem_error == 1)
	{
		printf("| insufficient memory \n");
	}
	if (error != 0)
	{
		//printf("\nDouble precision Helmholtz example on a whole sphere has ");
		//printf("FAILED to compute the solution...\n");
	}
	/* Success message to print if everything is OK */
	//printf("\n Double precision Helmholtz example on a whole sphere has ");
	//printf("successfully PASSED\n through all steps of computation!\n");

	//printf("press any key to continue...");
	//getchar();
}

//计算projection中的矩阵A
void SphereSolver::computeA()
{
	
}


