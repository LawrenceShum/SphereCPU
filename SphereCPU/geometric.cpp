#include "SphereSolver.h"
#include "TotalInclude.h"
#include "CubicSolver.h"

using namespace std;

void SphereSolver::geometric()
{
	for (int y = 1; y <= n_theta - 1; y++)
	{
		for (int x = 0; x < n_phi; x++)
		{
			velocity_geometric(x, y);
		}
	}

	swap_vel_phi();
	swap_vel_theta();

	//将next的速度清零
	for (int y = 0; y < n_theta + 1; y++)
	{
		for (int x = 0; x < n_phi; x++)
		{
			int num = n_phi*y + x;

			vel_phi_next[num] = 0.0;
			vel_theta_next[num] = 0.0;
		}
	}
}

void SphereSolver::velocity_geometric(int x, int y)
{
	int num = x + y*n_phi;
	float coTheta = y * gridLen;
	float coPhi = x * gridLen;

	//u是phi方向，v是theta方向
	float uPrev = get_vel_phi(y, x);
	float vPrev = get_vel_theta(y, x);

	float G = 0.0;
	float uNext = 0.0;
	if (abs(coTheta - M_PI / 2.0) < 1e-8)
	{
		G = 0.0;
		uNext = uPrev;
	}
	else
	{
		G = dt * cos(coTheta) / (radius * sin(coTheta));
		float cof = G * G;
		float A = 0.0;
		float B = (G * vPrev + 1.0) / cof;
		float C = -uPrev / cof;
		
		double solution[3];
		//解三次方程，x^3 + A*x^2 + B*x + C
		SolveP3(solution, A, B, C);

		uNext = static_cast<float>(solution[0]);
	}
	float vNext = vPrev + G * uNext * uNext;

	vel_phi_next[num] = uNext;
	vel_theta_next[num] = vNext;

	//极点处在这一步骤上没有做geometric
	//这里将极点处this的量复制到next处
	for (int x = 0; x < n_phi; x++)
	{
		vel_phi_next[x] = vel_phi_this[x];
		vel_theta_next[x] = vel_theta_this[x];

		vel_phi_next[x + n_theta*n_phi] = vel_phi_this[x + n_theta*n_phi];
		vel_theta_next[x + n_theta*n_phi] = vel_theta_this[x + n_theta*n_phi];
	}
}