#include "TotalInclude.h"
#include "SphereSolver.h"

using namespace std;

//spectual filter对极点上的速度场进行操作
void SphereSolver::spectual_filter()
{
	float ux_N = 0.0;
	float uy_N = 0.0;
	//北极点上
	for (int k = 0; k < n_phi; k++)
	{
		float u_theta = vel_theta_this[k];
		ux_N = u_theta*cos(k*gridLen);
		uy_N = u_theta*sin(k*gridLen);
	}
	ux_N = (2 / n_theta)*ux_N;
	uy_N = (2 / n_theta)*uy_N;

	//利用ux和uy来重构极点处的utheta和uphi
	for (int k = 0; k < n_phi; k++)
	{
		vel_theta_next[k] = ux_N*cos(k*gridLen) + uy_N*sin(k*gridLen);
	}
	for (int k = 0; k < n_phi; k++)
	{
		//求余操作
		int num = (k + static_cast<int>(n_phi / 2.0)) % n_phi;
		vel_phi_next[k] = vel_theta_next[num];
	}

	//南极点上
	float ux_S = 0.0;
	float uy_S = 0.0;
	//北极点上
	for (int k = 0; k < n_phi; k++)
	{
		int num = k + n_phi*n_theta;
		float u_theta = vel_theta_this[num];
		ux_N = u_theta*cos(k*gridLen);
		uy_N = u_theta*sin(k*gridLen);
	}
	ux_S = (2 / n_theta)*ux_S;
	uy_S = (2 / n_theta)*uy_S;

	//利用ux和uy来重构极点处的utheta和uphi
	for (int k = 0; k < n_phi; k++)
	{
		int num = k + n_phi*n_theta;
		vel_theta_next[num] = ux_S*cos(k*gridLen) + uy_S*sin(k*gridLen);
	}
	for (int k = 0; k < n_phi; k++)
	{
		//求余操作
		int num = (k + static_cast<int>(n_phi / 2.0)) % n_phi;
		vel_phi_next[k + n_phi*n_theta] = -vel_theta_next[num + n_phi*n_theta];
	}

	//将极点处next的值复制到this中，清空next
	for (int k = 0; k < n_phi; k++)
	{
		vel_phi_this[k] = vel_phi_next[k];
		vel_phi_this[k + n_phi*n_theta] = vel_phi_next[k + n_phi*n_theta];

		vel_phi_next[k] = 0.0;
		vel_phi_next[k + n_phi*n_theta];
	}
}