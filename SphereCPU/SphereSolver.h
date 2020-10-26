#pragma once
#include "TotalInclude.h"
#include "SphereQuantity.h"

class SphereSolver
{
private:
	//advection
	void advection();
	void velocity_advection(int, int);
	void polarAdvection();
	//geometric
	void geometric();
	void velocity_geometric(int, int);
	//spectual filter
	void spectual_filter();
	//projection
	void projection();
	void MKLSpherical();
	void solvePolarProjection();
	void computeA();
	//deffusion
	void diffusion();

	float sampleAt(float, float, float*);

	//timestep
	float dt;
	//每帧间隔时间
	float DT;
	//float theta;
	//float phi;
	//半径
	float radius;

	//仿真中规定theta方向的网格数量
	//一般，phi方向的数量是theta方向的2倍
	int n_theta;
	//仿真中规定phi方向的网格数量
	int n_phi;

	int simulation_time;

	//网格边长对应的角度
	float gridLen;
	float invGridLen;

	//当前计算步骤phi方向的速度
	float* vel_phi_this;
	//计算后phi方向的速度
	float* vel_phi_next;
	//当前计算步骤theta方向的速度
	float* vel_theta_this;
	//计算后theta方向的速度
	float* vel_theta_next;
	//当前计算步骤的压力
	float* presure_this;
	//计算后的压力
	float* presure_next;

	void initialize_velocity();
	void initialize_presure();

public:
	SphereSolver(int n_theta = 128, int n_phi = 256, float dt = 0.01, float DT = 0.05, float radius = 100);
	~SphereSolver();

	//前进一个timestep的计算
	//包括5个步骤
	//advection, projection, geometric, spectual filter, projection
	void step(float);

//	static void solvecubic(double, double, double, double);

	float get_vel_phi(int, int);
	float get_vel_theta(int, int);
	float get_presure(int, int);
	void set_vel_phi(int, int, float);
	void set_vel_theta(int, int, float);
	void set_presure(int, int, float);

	//this和next之间的互换
	void swap_vel_phi();
	void swap_vel_theta();
	void swap_presure();

	void write_vel();
	void write_presure();
	float Lerp(float, float, float);

	friend class Sphere;
	friend class LinearSolver;
	friend class Drawer;
};
