#include "SphereSolver.h"
#include <iostream>
#include <fstream>

using namespace std;

SphereSolver::SphereSolver(int n_theta, int n_phi, float dt, float DT, float radius):
	n_theta(n_theta),n_phi(n_phi),dt(dt),DT(DT),radius(radius),gridLen(M_2PI/n_phi),invGridLen(1.0/gridLen)
{
	vel_phi_this = new float[n_phi*(n_theta+1)];
	vel_theta_this = new float[n_phi*(n_theta+1)];
	vel_phi_next = new float[n_phi*(n_theta+1)];
	vel_theta_next = new float[n_phi*(n_theta+1)];
	presure_this = new float[n_phi*(n_theta+1)];
	presure_next = new float[n_phi*(n_theta+1)];

	simulation_time = 100;

	initialize_velocity();
	initialize_presure();
}

SphereSolver::~SphereSolver()
{
	delete vel_phi_this;
	delete vel_theta_this;
	delete presure_this;
	delete vel_phi_next;
	delete vel_theta_next;
	delete presure_next;
}

//模拟步骤
void SphereSolver::step(float dt)
{
	float timestep = dt;

	//step 1
	advection();

	//cout << "hello" << endl;
	//steo 2
	geometric();

	//step 3 
	spectual_filter();

	//step 4
	projection();

	//将模拟结果以图片的形式输出
	//输出的图片格式为png

}

void SphereSolver::initialize_velocity()
{
	for (int j = 0; j < n_theta + 1; j++)
	{
		for (int i = 0; i < n_phi; i++)
		{
			this->set_vel_phi(j, i, 0.5);
			this->set_vel_theta(j, i, 1.0);
		}
	}
}

void SphereSolver::initialize_presure()
{
	for (int j = 0; j < n_theta+1; j++)
	{
		for (int i = 0; i < n_phi; i++)
		{
			this->set_presure(j, i, 0.0);
			this->set_presure(j, i, 0.0);
		}
	}
}

void SphereSolver::swap_vel_phi()
{
	float* mid = vel_phi_this;
	
	vel_phi_this = vel_phi_next;
	vel_phi_next = mid;
}

void SphereSolver::swap_vel_theta()
{
	float* mid = vel_theta_this;

	vel_theta_this = vel_theta_next;
	vel_theta_next = mid;
}

void SphereSolver::swap_presure()
{
	float* mid = presure_this;
	presure_this = presure_next;
	presure_next = mid;
}

float SphereSolver::get_vel_phi(int theta_index, int phi_index)
{
	return vel_phi_this[n_phi*theta_index + phi_index];
}

float SphereSolver::get_vel_theta(int theta_index, int phi_index)
{
	return vel_theta_this[n_phi*theta_index + phi_index];
}

float SphereSolver::get_presure(int theta_index, int phi_index)
{
	return presure_this[n_phi*theta_index + phi_index];
}

void SphereSolver::set_vel_phi(int theta_index, int phi_index, float value)
{
	vel_phi_this[n_phi*theta_index + phi_index] = value;
}

void SphereSolver::set_vel_theta(int theta_index, int phi_index, float value)
{
	vel_theta_this[n_phi*theta_index + phi_index] = value;
}

void SphereSolver::set_presure(int theta_index, int phi_index, float value)
{
	presure_this[n_phi*theta_index + phi_index] = value;
}

float SphereSolver::sampleAt(float coPhi, float coTheta, float* u)
{
	//首先将phi、theta都变成在0~2pi、0~pi的区间内
	int loops = static_cast<int>(std::floor(coTheta / M_2PI));
	coTheta = coTheta - loops * M_2PI;

	bool isFlipped = false;
	//穿过南极，phi向前推进pi
	if (coTheta > M_PI)
	{
		coTheta = M_2PI - coTheta;
		coPhi += M_PI;
		isFlipped = true;
	}
	loops = static_cast<int>(std::floor(coPhi / M_2PI));
	coPhi = coPhi - loops * M_2PI;

	float normedPhi = coPhi * invGridLen;
	float normedTheta = coTheta * invGridLen;

	int phiIndex = static_cast<int>(std::floor(normedPhi));
	int thetaIndex = static_cast<int>(std::floor(normedTheta));

	float alphaPhi = normedPhi - static_cast<float>(phiIndex);
	float alphaTheta = normedTheta - static_cast<float>(thetaIndex);

	int phiLower = phiIndex % n_phi;
	int phiHigher = (phiLower + 1) % n_phi;;
	int thetaLower = thetaIndex;
	int thetaHigher = thetaIndex + 1;

	float lowerBelt = Lerp(u[phiLower+thetaLower*n_phi], u[phiHigher+thetaLower*n_phi], alphaPhi);
	//float lowerBelt = Lerp<float>(get_vel_phi(phiLower, thetaLower), get_vel_phi(phiHigher, thetaLower), alphaPhi);
	float higherBelt = Lerp(u[phiLower + thetaHigher*n_phi], u[phiHigher + thetaHigher*n_phi], alphaPhi);
	//float higherBelt = Lerp<float>(get_vel_phi(phiLower, thetaHigher), get_vel_phi(phiHigher, thetaHigher), alphaPhi);

	return Lerp(lowerBelt, higherBelt, alphaTheta);
}

//插值
float SphereSolver::Lerp(float fromEndPoint, float toEndPoint, float factor)
{
	return (1.0 - factor) * fromEndPoint + factor * toEndPoint;
}

//将计算得到的速度写入txt文件内
void SphereSolver::write_vel()
{
	ofstream out_phi("velocity_phi.txt");
	if (out_phi.is_open())
	{
		for (int y = 0; y < n_theta+1; y++)
		{
			for (int x = 0; x < n_phi; x++)
			{
				float u = get_vel_phi(y, x);
				out_phi << u << "    ";
			}
			out_phi << endl;
		}
		out_phi.close();
	}
	else cout << "Failed to open the velocity_phi.txt file" << endl;

	ofstream out_theta("velocity_theta.txt");
	if (out_theta.is_open())
	{
		for (int y = 0; y < n_theta + 1; y++)
		{
			for (int x = 0; x < n_phi; x++)
			{
				float v = get_vel_theta(y, x);
				out_phi << v << "    ";
			}
			out_theta << endl;
		}
		out_theta.close();
	}
	else cout << "Failed to open the velocity_theta.txt file" << endl;
}

//将计算得到的presure写入txt文件内
void SphereSolver::write_presure()
{
	ofstream out_presure("presure.txt");
	if (out_presure.is_open())
	{
		for (int y = 0; y < n_theta + 1; y++)
		{
			for (int x = 0; x < n_phi; x++)
			{
				float presure = get_presure(y, x);
				out_presure << presure << "    ";
			}
			out_presure << endl;
		}
		out_presure.close();
	}
	else cout << "Failed to open the velocity_phi.txt file" << endl;
}