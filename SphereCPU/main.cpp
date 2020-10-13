#include "Sphere.h"
#include "TotalInclude.h"
#include <iostream>

using namespace std;

int main()
{
	//n_phi建议是4的倍数
	//n_theta必须是n_phi的1/2
	int n_theta = 6;
	int n_phi = 12;
	float dt = 0.01;
	float DT = 0.05;
	float radius = 10;
	double A = 0;
	double B = 0;
	double C = 0;
	double D = 0;
	Sphere sim(n_theta,n_phi,dt,DT,radius,A,B,C,D);

	sim.start();
	sim.~Sphere();
	return 0;
}