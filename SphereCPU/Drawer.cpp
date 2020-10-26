#include "Drawer.h"
#include "TotalInclude.h"

using namespace std;

Drawer::Drawer(unsigned width, unsigned height, SphereSolver* Psolver) :
	width(width), height(height), Psolver(Psolver),Angle2Pixel(width/M_2PI)
{
	image.resize(width * height * 4);
	particles.resize(num_particles);
	color_old.resize(width*height);
	color_new.resize(width*height);
	//initial_particles();
	initial(image);
	//draw_particles();
	
	output_png("0000.png");
}

Drawer::~Drawer()
{

}

void Drawer::initial(vector<unsigned char>& image)
{
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			if (y < (height / 16))
				color_old[y*width + x] = red;
			//set_point(image, x, y, red);
			else if (y < (height / 8))
				color_old[y*width + x] = white;
			//set_point(image, x, y, white);
			else if (y < (height * 3 / 16))
				color_old[y*width + x] = red;
			//set_point(image, x, y, red);
			else if (y < (height / 4))
				color_old[y*width + x] = white;
			//set_point(image, x, y, white);
			else if (y < (height * 5 / 16))
				color_old[y*width + x] = red;
			//set_point(image, x, y, red);
			else if (y < (height * 3 / 8))
				color_old[y*width + x] = white;
			//set_point(image, x, y, white);
			else if (y < (height * 7 / 16))
				color_old[y*width + x] = red;
			//set_point(image, x, y, red);
			else if (y < (height / 2))
				color_old[y*width + x] = white;
			else if (y < (height * 9 / 16))
				color_old[y*width + x] = red;
			else if (y < (height * 5 / 8))
				color_old[y*width + x] = white;
			else if (y < (height * 11 / 16))
				color_old[y*width + x] = red;
			else if (y < (height * 3 / 4))
				color_old[y*width + x] = white;
			else if (y < (height * 13 / 16))
				color_old[y*width + x] = red;
			else if (y < (height * 7 / 8))
				color_old[y*width + x] = white;
			else if (y < (height * 15 / 16))
				color_old[y*width + x] = red;
			else
				color_old[y*width + x] = white;
				//set_point(image, x, y, white);
		}
	}

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			set_point(image, x, y, color_old[y*width + x]);
		}
	}
}

void Drawer::initial_particles()
{
	for (unsigned long i = 0; i < num_particles; i++)
	{
		particles[i].x = 200 + (0 + i)%100;
		particles[i].y = 200 + std::floor(i / 100);
	}
}

//计算粒子运动
void Drawer::calculate_particles()
{
	double Pixel2Angle = 1.0 / Angle2Pixel;

	for (unsigned long i = 0; i < num_particles; i++)
	{
		//获取粒子当前所在位置的速度，通过插值的方法
		//并计算下一时间步，粒子的位置
		float phi_this, phi_next, theta_this, theta_next;

		phi_this = particles[i].x * Pixel2Angle;
		theta_this = particles[i].y * Pixel2Angle;

		//cout << phi_this << "   " << theta_this << endl;

		float* u = Psolver->vel_phi_this;
		float* v = Psolver->vel_theta_this;

		float velPhi = Psolver->sampleAt(phi_this, theta_this, u);
		float velTheta = Psolver->sampleAt(phi_this, theta_this, v);

		//更新粒子的位置
		//phi_next = velPhi * Psolver->dt + phi_this;
		//theta_next = velTheta * Psolver->dt + theta_this;
		phi_next = phi_this + velPhi * Psolver->dt;
		theta_next = theta_this + velTheta * Psolver->dt;

		//cout << phi_next << "   " << theta_next << endl;

		//检测坐标是否超出球坐标系的定义域
		//首先将phi、theta都变成在0~2pi、0~pi的区间内
		int loops = static_cast<int>(std::floor(theta_next / M_2PI));
		theta_next = theta_next - loops * M_2PI;

		bool isFlipped = false;
		//穿过南极，phi向前推进pi
		if (theta_next > M_PI)
		{
			theta_next = M_2PI - theta_next;
			phi_next += M_PI;
			isFlipped = true;
		}
		else if (theta_next < 0)
		{
			theta_next = theta_next - 2 * theta_next;
			phi_next += M_PI;
			isFlipped = true;
		}
		loops = static_cast<int>(std::floor(phi_next / M_2PI));
		phi_next = phi_next - loops * M_2PI;

		//cout << particles[i].x << "   " << particles[i].y << endl;

		//更新粒子的位置
		particles[i].x = std::floor(phi_next * Angle2Pixel);
		particles[i].y = std::floor(theta_next * Angle2Pixel);

		//cout << particles[i].x << "   " << particles[i].y << endl;

	}
}

//画粒子
void Drawer::draw_particles()
{
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			set_point(image, x, y, black);
		}
	}

	for (unsigned long i = 0; i < num_particles; i++)
	{
		set_point(image, particles[i].x, particles[i].y, red);
	}
}

//画画
void Drawer::draw()
{		
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			set_point(image, x, y, color_new[y*width + x]);
		}
	}
}

//更新颜色
void Drawer::update_color()
{
	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{
			color_old[y*width + x] = color_new[y * width + x];
		}
	}
}

//颜色对流
void Drawer::color_advect()
{
	double Pixel2Angle = 1.0 / Angle2Pixel;

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{

			/********************向前型对流******************/
			float phi_this, phi_old, theta_this, theta_old;
			phi_this = x * Pixel2Angle;
			theta_this = y * Pixel2Angle;

			//在当前像素位置上采样速度场的速度
			float* u = Psolver->vel_phi_this;
			float* v = Psolver->vel_theta_this;

			float vel_phi = Psolver->sampleAt(phi_this, theta_this, u);
			float vel_theta = Psolver->sampleAt(phi_this, theta_this, v);

			//根据采样得到的速度进行回溯，得到上一时间步该点的颜色，更新颜色
			/**********semi-lagrangian advection*********/
			//phi_old = phi_this - vel_phi * Psolver->dt;
			//theta_old = theta_this - vel_theta * Psolver->dt;

			/**********RK2 advection*********/
			float phi_mid = phi_this - vel_phi * Psolver->dt * 0.5;
			float theta_mid = theta_this - vel_phi * Psolver->dt * 0.5;

			phi_old = phi_this - Psolver->sampleAt(phi_mid, theta_mid, u) * Psolver->dt;
			theta_old = theta_this - Psolver->sampleAt(phi_mid, theta_mid, v) * Psolver->dt;

			//检测坐标是否超出球坐标系的定义域
			//首先将phi、theta都变成在0~2pi、0~pi的区间内
			int loops = static_cast<int>(std::floor(theta_old / M_2PI));
			theta_old = theta_old - loops * M_2PI;

			bool isFlipped = false;
			//穿过南极，phi向前推进pi
			if (theta_old > M_PI)
			{
				theta_old = M_2PI - theta_old;
				phi_old += M_PI;
				isFlipped = true;
			}
			/*else if (theta_old < 0)
			{
				theta_old = theta_old - 2 * theta_old;
				phi_old += M_PI;
				isFlipped = true;
			}*/
			loops = static_cast<int>(std::floor(phi_old / M_2PI));
			phi_old = phi_old - loops * M_2PI;

			//转换成像素坐标
			//int x_old = static_cast<int>(phi_old * Angle2Pixel);
			//int y_old = static_cast<int>(theta_old * Angle2Pixel);
			int x_old = std::floor(phi_old * Angle2Pixel);
			int y_old = std::floor(theta_old * Angle2Pixel);

			//cout << phi_old << "   " << theta_old << endl;
			//color_new[y*width + x] = color_old[y_old*width + x_old];
			color_new[y * width + x] = color_old[y_old * width + x_old];
			

			/*******************向后型对流*********************/
			/*float phi_this, phi_next, theta_this, theta_next;
			phi_this = x * Pixel2Angle;
			theta_this = y * Pixel2Angle;

			//在当前像素位置上采样速度场的速度
			float* u = Psolver->vel_phi_this;
			float* v = Psolver->vel_theta_this;

			float vel_phi = Psolver->sampleAt(phi_this, theta_this, u);
			float vel_theta = Psolver->sampleAt(phi_this, theta_this, v);

			//根据采样得到的速度进行回溯，得到上一时间步该点的颜色，更新颜色
			phi_next = phi_this + vel_phi * Psolver->dt;
			theta_next = theta_this + vel_theta * Psolver->dt;

			//检测坐标是否超出球坐标系的定义域
			//首先将phi、theta都变成在0~2pi、0~pi的区间内
			int loops = static_cast<int>(std::floor(theta_next / M_2PI));
			theta_next = theta_next - loops * M_2PI;

			bool isFlipped = false;
			//穿过南极，phi向前推进pi
			if (theta_next > M_PI)
			{
				theta_next = M_2PI - theta_next;
				phi_next += M_PI;
				isFlipped = true;
			}
			else if (theta_next < 0)
			{
				theta_next = -theta_next;
				phi_next += M_PI;
				isFlipped = true;
			}
			loops = static_cast<int>(phi_next / M_2PI);
			phi_next = phi_next - loops * M_2PI;

			//转换成像素坐标
			//int x_old = static_cast<int>(phi_old * Angle2Pixel);
			//int y_old = static_cast<int>(theta_old * Angle2Pixel);
			int x_next = static_cast<int>(phi_next * Angle2Pixel);
			int y_next = static_cast<int>(theta_next * Angle2Pixel);

			set_point(image, x_next, y_next, color_old[y*width + x]);*/
		}
	}
	update_color();
}

void Drawer::set_point(vector<unsigned char>& image, int x, int y, color color)
{
	//int x = x_i;
	//int y = height - y_i;
	image[4 * width * y + 4 * x + 0] = color.R;
	image[4 * width * y + 4 * x + 1] = color.G;
	image[4 * width * y + 4 * x + 2] = color.B;
	image[4 * width * y + 4 * x + 3] = color.A;
}


vector<unsigned char> Drawer::get_image() 
{
	return image;
}


//输出png图片
void Drawer::output_png(const char* filename)
{
	//Encode the image
	unsigned error = lodepng::encode(filename, image, width, height);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}