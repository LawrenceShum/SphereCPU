#include "Sphere.h"
#include "SphereSolver.h"
#include "SphereQuantity.h"
#include <iostream>
#include <fstream>

using namespace std;

Sphere::Sphere(int n_theta, int n_phi, float dt, float DT, float radius, double A,
	double B, double C, double D)
	:n_theta(n_theta), n_phi(n_phi),dt(dt),DT(DT),radius(radius),A(A),B(B),C(C),D(D)
{
	simulation_time = 100;
}

Sphere::~Sphere()
{

}

//开始
void Sphere::start()
{
	SphereSolver solver(n_theta, n_phi, dt, DT, radius);

	float T = 0.0;

	ofstream out_phi("velocity_phi.txt");
	ofstream out_theta("velocity_theta.txt");
	ofstream presure("presure.txt");

	//开始仿真
	for (int i = 1; i <= simulation_time; i++)
	{
		//前进一个时间步
		while (T < i*DT)
		{
			solver.step(dt);
			T += dt;
		}
		solver.step(dt + i*DT - T);
		T = i*DT;

		//写入速度
		if (out_phi.is_open())
		{
			out_phi << "Frame " << i << endl;
			for (int y = 0; y < n_theta + 1; y++)
			{
				for (int x = 0; x < n_phi; x++)
				{
					float u = solver.get_vel_phi(y, x);
					out_phi << u << "    ";
				}
				out_phi << endl;
			}
		}
		else cout << "Failed to open velocity_phi.txt file" << endl;
		out_phi << endl;

		if (out_theta.is_open())
		{
			out_theta << "Frame " << i << endl;
			for (int y = 0; y < n_theta + 1; y++)
			{
				for (int x = 0; x < n_phi; x++)
				{
					float v = solver.get_vel_theta(y, x);
					out_theta << v << "    ";
				}
				out_theta << endl;
			}
		}
		else cout << "Failed to open velocity_theta.txt file" << endl;
		out_theta << endl;
	}
	out_phi.close();
	out_theta.close();
	solver.~SphereSolver();
}
