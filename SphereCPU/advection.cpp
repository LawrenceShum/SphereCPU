#include "SphereSolver.h"

#include "TotalInclude.h"

using namespace std;

void SphereSolver::advection()
{
	for (int y = 1; y <= n_theta - 1; y++)
	{
		for (int x = 0; x < n_phi; x++)
		{
			velocity_advection(x, y);
		}
	}

	polarAdvection();

	swap_vel_theta();
	swap_vel_phi();

	//��next���ٶ�����
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


void SphereSolver::velocity_advection(int x, int y)
{
	int num = y*n_phi + x;

	//ת����������ϵ�µ����꣬������
	float coPhi = x * this->gridLen;
	float coTheta = y * this->gridLen;

	//������ǰ�ٶ�
	float uPhi = sampleAt(coPhi, coTheta, vel_phi_this);
	float uTheta = sampleAt(coPhi, coTheta, vel_theta_this);

	float latRadius = this->radius * sin(coTheta);
	float cofPhi = dt / latRadius;
	float cofTheta = dt / radius;

	//����advection��deltaֵ
	float deltaPhi = cofPhi * uPhi;
	float deltaTheta = cofTheta * uTheta;

	//semi-lagrangian��advection
	float pPhi = coPhi - deltaPhi;
	float pTheta = coTheta - deltaTheta;

	float advectedVel_phi = sampleAt(pPhi, pTheta, vel_phi_this);
	float advectedVel_theta = sampleAt(pPhi, pTheta, vel_theta_this);

	this->vel_phi_next[num] = advectedVel_phi;
	this->vel_theta_next[num] = advectedVel_theta;
}

/*void SphereSolver::polarAdvection()
{
	//�ڱ�������theta = 0
	for (int x = 0; x < n_phi; x++)
	{
		float coPhi = x * this->gridLen;
		float coTheta = 0;

		if (x == 0)
		{
			float mid_vel = 0.5*(vel_theta_this[n_phi - 1] + vel_theta_this[1]);
			float derive_vel = 0.5*(vel_theta_this[n_phi - 1] - vel_theta_this[1]);
			vel_theta_next[x] = -(vel_theta_this[x] / radius)*derive_vel*(dt / gridLen) + mid_vel;
		}
		else if (x == n_phi - 1)
		{
			float mid_vel = 0.5*(vel_theta_this[n_phi - 2] + vel_theta_this[0]);
			float derive_vel = 0.5*(vel_theta_this[n_phi - 2] - vel_theta_this[0]);
			vel_theta_next[x] = -(vel_theta_this[x] / radius)*derive_vel*(dt / gridLen) + mid_vel;
		}
		else
		{
			float mid_vel = 0.5*(vel_theta_this[x + 1] + vel_theta_this[x - 1]);
			float derive_vel = 0.5*(vel_theta_this[x + 1] - vel_theta_this[x - 1]);
			vel_theta_next[x] = -(vel_theta_this[x] / radius)*derive_vel*(dt / gridLen) + mid_vel;
		}
	}
	//���ݼ���ı߽�����������������ٶ��ǹ���ģ���˼���phi������ٶ�
	for (int x = 0; x < n_phi; x++)
	{
		int conjugate = (x + static_cast<int>(n_phi / 4)) % n_phi;
		vel_phi_next[x] = vel_theta_next[conjugate];
	}

	//���ϼ�����theta = pi
	for (int x = 0; x < n_phi; x++)
	{
		float coPhi = x * this->gridLen;
		float coTheta = 0;

		if (x == 0)
		{
			float mid_vel = 0.5*(vel_theta_this[n_phi - 1 + n_theta*n_phi] + vel_theta_this[1 + n_theta*n_phi]);
			float derive_vel = 0.5*(vel_theta_this[n_phi - 1 + n_theta*n_phi] - vel_theta_this[1 + n_theta*n_phi]);
			vel_theta_next[x] = -(vel_theta_this[x] / radius)*derive_vel*(dt / gridLen) + mid_vel;
		}
		else if (x == n_phi - 1)
		{
			float mid_vel = 0.5*(vel_theta_this[n_phi - 2 + n_theta*n_phi] + vel_theta_this[0 + n_theta*n_phi]);
			float derive_vel = 0.5*(vel_theta_this[n_phi - 2 + n_theta*n_phi] - vel_theta_this[0 + n_theta*n_phi]);
			vel_theta_next[x] = -(vel_theta_this[x] / radius)*derive_vel*(dt / gridLen) + mid_vel;
		}
		else
		{
			float mid_vel = 0.5*(vel_theta_this[x + 1 + n_theta*n_phi] + vel_theta_this[x - 1 + n_theta*n_phi]);
			float derive_vel = 0.5*(vel_theta_this[x + 1 + n_theta*n_phi] - vel_theta_this[x - 1 + n_theta*n_phi]);
			vel_theta_next[x] = -(vel_theta_this[x] / radius)*derive_vel*(dt / gridLen) + mid_vel;
		}
	}
	//���ݹ����ϵ����phi�����ϵ��ٶ�
	for (int x = 0; x < n_phi; x++)
	{
		int conjugate = (x + static_cast<int>(n_phi / 4)) % n_phi;
		vel_phi_next[x + n_phi*n_theta] = -vel_theta_next[conjugate + n_phi*n_theta];
	}
}*/

//version 2 ���㴦����
void SphereSolver::polarAdvection()
{
	float invRadius = 1.0 / radius;
	//�ڱ����㴦
	for (int x = 0; x < n_phi; x++)
	{
		float coTheta = 0.0;
		float coPhi = x*gridLen;

		float uTheta = vel_theta_this[x];

		float pTheta = coTheta - dt*uTheta;
		float pPhi = coPhi;

		float midTheta = sampleAt(pPhi, pTheta, vel_theta_this);
		vel_theta_next[x] = uTheta + invRadius * (midTheta - uTheta);
	}
	//���ݹ����ϵ����phi������ٶ�
	for (int x = 0; x < n_phi; x++)
	{
		int conjugate = (x + static_cast<int>(n_phi / 4)) % n_phi;
		vel_phi_next[x] = vel_theta_next[conjugate];
	}

	//���ϼ��㴦
	for (int x = 0; x < n_phi; x++)
	{
		float coTheta = M_2PI;
		float coPhi = x*gridLen;

		float uTheta = vel_theta_this[x + n_theta*n_phi];

		float pTheta = coTheta - dt*uTheta;
		float pPhi = coPhi;

		float midTheta = sampleAt(pPhi, pTheta, vel_theta_this);
		vel_theta_next[x + n_theta*n_phi] = uTheta + invRadius * (midTheta - uTheta);
		//vel_theta_next[x + n_theta*n_phi] = sampleAt(pPhi, pTheta, vel_theta_this);
	}
	//���ݹ����ϵ����phi������ٶ�
	for (int x = 0; x < n_phi; x++)
	{
		int conjugate = (x + static_cast<int>(n_phi / 4)) % n_phi;
		vel_phi_next[x + n_theta*n_phi] = vel_theta_next[conjugate + n_theta*n_phi];
	}
}