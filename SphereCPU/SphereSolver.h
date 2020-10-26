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
	//ÿ֡���ʱ��
	float DT;
	//float theta;
	//float phi;
	//�뾶
	float radius;

	//�����й涨theta�������������
	//һ�㣬phi�����������theta�����2��
	int n_theta;
	//�����й涨phi�������������
	int n_phi;

	int simulation_time;

	//����߳���Ӧ�ĽǶ�
	float gridLen;
	float invGridLen;

	//��ǰ���㲽��phi������ٶ�
	float* vel_phi_this;
	//�����phi������ٶ�
	float* vel_phi_next;
	//��ǰ���㲽��theta������ٶ�
	float* vel_theta_this;
	//�����theta������ٶ�
	float* vel_theta_next;
	//��ǰ���㲽���ѹ��
	float* presure_this;
	//������ѹ��
	float* presure_next;

	void initialize_velocity();
	void initialize_presure();

public:
	SphereSolver(int n_theta = 128, int n_phi = 256, float dt = 0.01, float DT = 0.05, float radius = 100);
	~SphereSolver();

	//ǰ��һ��timestep�ļ���
	//����5������
	//advection, projection, geometric, spectual filter, projection
	void step(float);

//	static void solvecubic(double, double, double, double);

	float get_vel_phi(int, int);
	float get_vel_theta(int, int);
	float get_presure(int, int);
	void set_vel_phi(int, int, float);
	void set_vel_theta(int, int, float);
	void set_presure(int, int, float);

	//this��next֮��Ļ���
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
