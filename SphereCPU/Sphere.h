#pragma once
#include "TotalInclude.h"

class Sphere
{
private:
	float radius;
	int n_theta;
	int n_phi;
	float gridLen;

	//timestep
	float dt;
	//帧与帧间隔时间
	float DT;

	int simulation_time;

	double A, B, C, D;

public:
	Sphere(int n_theta, int n_phi, float dt, float DT, float radius, double A,
		double B, double C, double D);
	~Sphere();
	void start();
};
