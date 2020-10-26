#include "Sphere.h"
#include "SphereSolver.h"
#include "SphereQuantity.h"
#include <iostream>
#include <fstream>
#include "Drawer.h"

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

//��ʼ
void Sphere::start()
{
	int n_image = 1;

	SphereSolver solver(n_theta, n_phi, dt, DT, radius);
	float T = 0.0;

	//width, height
	Drawer draw(800, 400, &solver);

	ofstream out_phi("velocity_phi.txt");
	ofstream out_theta("velocity_theta.txt");
	ofstream presure("presure.txt");

	/*char p[] = "0000.png";
	unsigned thousand = (unsigned)(n_image / 1000);
	p[0] = (char)('0' + thousand);
	unsigned hundred = (unsigned)((n_image - thousand * 1000) / 100);
	p[1] = (char)('0' + hundred);
	unsigned ten = (unsigned)((n_image - thousand * 1000 - hundred * 100) / 10);
	p[2] = (char)('0' + ten);
	unsigned unit = (unsigned)(n_image - thousand * 1000 - hundred * 100 - ten * 10);
	p[3] = (char)('0' + unit);

	const char* filename = p;*/

	//��ʼ����
	for (int i = 1; i <= simulation_time; i++)
	{
		//ǰ��һ��ʱ�䲽
		while (T < i*DT)
		{
			solver.step(dt);
			T += dt;
			
			//��ʼ��ͼ
			char p[] = "0000.png";
			unsigned thousand = (unsigned)(n_image / 1000);
			p[0] = (char)('0' + thousand);
			unsigned hundred = (unsigned)((n_image - thousand * 1000) / 100);
			p[1] = (char)('0' + hundred);
			unsigned ten = (unsigned)((n_image - thousand * 1000 - hundred * 100) / 10);
			p[2] = (char)('0' + ten);
			unsigned unit = (unsigned)(n_image - thousand * 1000 - hundred * 100 - ten * 10);
			p[3] = (char)('0' + unit);
			const char* filename = p;

			//�������ӵ�λ��
			draw.calculate_particles();
			//�����ӻ��ϻ�����
			draw.draw_particles();

			//��ɫ����
			//draw.color_advect();
			//��ʼ����
			//draw.draw();
			//���pngͼƬ
			draw.output_png(filename);
			n_image++;
		}
		solver.step(dt + i*DT - T);
		T = i*DT;

		//��ʼ��ͼ
		char p[] = "0000.png";
		unsigned thousand = (unsigned)(n_image / 1000);
		p[0] = (char)('0' + thousand);
		unsigned hundred = (unsigned)((n_image - thousand * 1000) / 100);
		p[1] = (char)('0' + hundred);
		unsigned ten = (unsigned)((n_image - thousand * 1000 - hundred * 100) / 10);
		p[2] = (char)('0' + ten);
		unsigned unit = (unsigned)(n_image - thousand * 1000 - hundred * 100 - ten * 10);
		p[3] = (char)('0' + unit);
		const char* filename = p;

		//�������ӵ�λ��
		draw.calculate_particles();
		//�����ӻ��ϻ�����
		draw.draw_particles();

		//��ɫ����
		//draw.color_advect();
		//��ʼ����
		//draw.draw();
		//���pngͼƬ
		draw.output_png(filename);
		n_image++;


//*************д���ٶ�**************//
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
	draw.~Drawer();
}
