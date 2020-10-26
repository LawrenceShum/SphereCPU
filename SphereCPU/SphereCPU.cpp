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

//开始
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

	//开始仿真
	for (int i = 1; i <= simulation_time; i++)
	{
		//前进一个时间步
		while (T < i*DT)
		{
			solver.step(dt);
			T += dt;
			
			//开始画图
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

			//计算粒子的位置
			draw.calculate_particles();
			//将粒子画上画布上
			draw.draw_particles();

			//颜色对流
			//draw.color_advect();
			//开始画画
			//draw.draw();
			//输出png图片
			draw.output_png(filename);
			n_image++;
		}
		solver.step(dt + i*DT - T);
		T = i*DT;

		//开始画图
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

		//计算粒子的位置
		draw.calculate_particles();
		//将粒子画上画布上
		draw.draw_particles();

		//颜色对流
		//draw.color_advect();
		//开始画画
		//draw.draw();
		//输出png图片
		draw.output_png(filename);
		n_image++;


//*************写入速度**************//
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
