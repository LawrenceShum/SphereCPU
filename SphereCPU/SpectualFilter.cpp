#include "TotalInclude.h"
#include "SphereSolver.h"

using namespace std;

//spectual filter�Լ����ϵ��ٶȳ����в���
void SphereSolver::spectual_filter()
{
	float ux_N = 0.0;
	float uy_N = 0.0;
	//��������
	for (int k = 0; k < n_phi; k++)
	{
		float u_theta = vel_theta_this[k];

		ux_N += u_theta*cos(k*gridLen);
		uy_N += u_theta*sin(k*gridLen);
	}
	ux_N = (2 / n_theta) * ux_N;
	uy_N = (2 / n_theta) * uy_N;

	//����ux��uy���ع����㴦��utheta��uphi
	for (int k = 0; k < n_phi; k++)
	{
		float c = cos(k*gridLen);
		float s = sin(k*gridLen);

		vel_theta_next[k] = ux_N * c + uy_N * s;
		//�������
		vel_phi_next[k] = -ux_N * s + uy_N * c;
	}

	/*for (int k = 0; k < n_phi; k++)
	{
		//�������
		int num = (k + static_cast<int>(n_phi / 2.0)) % n_phi;
		vel_phi_next[k] = vel_theta_next[num];
	}*/

	//�ϼ�����
	float ux_S = 0.0;
	float uy_S = 0.0;
	//��������
	for (int k = 0; k < n_phi; k++)
	{
		int num = k + n_phi*n_theta;
		float u_theta = vel_theta_this[num];

		ux_S += u_theta*cos(k*gridLen);
		uy_S += u_theta*sin(k*gridLen);
	}

	ux_S = (2 / n_theta) * ux_S;
	uy_S = (2 / n_theta) * uy_S;

	//����ux��uy���ع����㴦��utheta��uphi
	for (int k = 0; k < n_phi; k++)
	{
		int num = k + n_phi*n_theta;
		float c = cos(k*gridLen);
		float s = sin(k*gridLen);

		vel_theta_next[num] = -ux_S * c + -uy_S * s;
		//�������
		vel_phi_next[num] = -ux_S * s + uy_S * c;
	}

	/*for (int k = 0; k < n_phi; k++)
	{
		//�������
		int num = (k + static_cast<int>(n_phi / 2.0)) % n_phi;
		vel_phi_next[k + n_phi*n_theta] = -vel_theta_next[num + n_phi*n_theta];
	}
	*/
	//�����㴦next��ֵ���Ƶ�this�У����next

	/*for (int k = 0; k < n_phi; k++)
	{
		vel_phi_this[k] = vel_phi_next[k];
		vel_phi_this[k + n_phi*n_theta] = vel_phi_next[k + n_phi*n_theta];

		vel_phi_next[k] = 0.0;
		vel_phi_next[k + n_phi*n_theta];
	}*/
}