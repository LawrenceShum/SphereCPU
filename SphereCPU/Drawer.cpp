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

//���������˶�
void Drawer::calculate_particles()
{
	double Pixel2Angle = 1.0 / Angle2Pixel;

	for (unsigned long i = 0; i < num_particles; i++)
	{
		//��ȡ���ӵ�ǰ����λ�õ��ٶȣ�ͨ����ֵ�ķ���
		//��������һʱ�䲽�����ӵ�λ��
		float phi_this, phi_next, theta_this, theta_next;

		phi_this = particles[i].x * Pixel2Angle;
		theta_this = particles[i].y * Pixel2Angle;

		//cout << phi_this << "   " << theta_this << endl;

		float* u = Psolver->vel_phi_this;
		float* v = Psolver->vel_theta_this;

		float velPhi = Psolver->sampleAt(phi_this, theta_this, u);
		float velTheta = Psolver->sampleAt(phi_this, theta_this, v);

		//�������ӵ�λ��
		//phi_next = velPhi * Psolver->dt + phi_this;
		//theta_next = velTheta * Psolver->dt + theta_this;
		phi_next = phi_this + velPhi * Psolver->dt;
		theta_next = theta_this + velTheta * Psolver->dt;

		//cout << phi_next << "   " << theta_next << endl;

		//��������Ƿ񳬳�������ϵ�Ķ�����
		//���Ƚ�phi��theta�������0~2pi��0~pi��������
		int loops = static_cast<int>(std::floor(theta_next / M_2PI));
		theta_next = theta_next - loops * M_2PI;

		bool isFlipped = false;
		//�����ϼ���phi��ǰ�ƽ�pi
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

		//�������ӵ�λ��
		particles[i].x = std::floor(phi_next * Angle2Pixel);
		particles[i].y = std::floor(theta_next * Angle2Pixel);

		//cout << particles[i].x << "   " << particles[i].y << endl;

	}
}

//������
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

//����
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

//������ɫ
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

//��ɫ����
void Drawer::color_advect()
{
	double Pixel2Angle = 1.0 / Angle2Pixel;

	for (int y = 0; y < height; y++)
	{
		for (int x = 0; x < width; x++)
		{

			/********************��ǰ�Ͷ���******************/
			float phi_this, phi_old, theta_this, theta_old;
			phi_this = x * Pixel2Angle;
			theta_this = y * Pixel2Angle;

			//�ڵ�ǰ����λ���ϲ����ٶȳ����ٶ�
			float* u = Psolver->vel_phi_this;
			float* v = Psolver->vel_theta_this;

			float vel_phi = Psolver->sampleAt(phi_this, theta_this, u);
			float vel_theta = Psolver->sampleAt(phi_this, theta_this, v);

			//���ݲ����õ����ٶȽ��л��ݣ��õ���һʱ�䲽�õ����ɫ��������ɫ
			/**********semi-lagrangian advection*********/
			//phi_old = phi_this - vel_phi * Psolver->dt;
			//theta_old = theta_this - vel_theta * Psolver->dt;

			/**********RK2 advection*********/
			float phi_mid = phi_this - vel_phi * Psolver->dt * 0.5;
			float theta_mid = theta_this - vel_phi * Psolver->dt * 0.5;

			phi_old = phi_this - Psolver->sampleAt(phi_mid, theta_mid, u) * Psolver->dt;
			theta_old = theta_this - Psolver->sampleAt(phi_mid, theta_mid, v) * Psolver->dt;

			//��������Ƿ񳬳�������ϵ�Ķ�����
			//���Ƚ�phi��theta�������0~2pi��0~pi��������
			int loops = static_cast<int>(std::floor(theta_old / M_2PI));
			theta_old = theta_old - loops * M_2PI;

			bool isFlipped = false;
			//�����ϼ���phi��ǰ�ƽ�pi
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

			//ת������������
			//int x_old = static_cast<int>(phi_old * Angle2Pixel);
			//int y_old = static_cast<int>(theta_old * Angle2Pixel);
			int x_old = std::floor(phi_old * Angle2Pixel);
			int y_old = std::floor(theta_old * Angle2Pixel);

			//cout << phi_old << "   " << theta_old << endl;
			//color_new[y*width + x] = color_old[y_old*width + x_old];
			color_new[y * width + x] = color_old[y_old * width + x_old];
			

			/*******************����Ͷ���*********************/
			/*float phi_this, phi_next, theta_this, theta_next;
			phi_this = x * Pixel2Angle;
			theta_this = y * Pixel2Angle;

			//�ڵ�ǰ����λ���ϲ����ٶȳ����ٶ�
			float* u = Psolver->vel_phi_this;
			float* v = Psolver->vel_theta_this;

			float vel_phi = Psolver->sampleAt(phi_this, theta_this, u);
			float vel_theta = Psolver->sampleAt(phi_this, theta_this, v);

			//���ݲ����õ����ٶȽ��л��ݣ��õ���һʱ�䲽�õ����ɫ��������ɫ
			phi_next = phi_this + vel_phi * Psolver->dt;
			theta_next = theta_this + vel_theta * Psolver->dt;

			//��������Ƿ񳬳�������ϵ�Ķ�����
			//���Ƚ�phi��theta�������0~2pi��0~pi��������
			int loops = static_cast<int>(std::floor(theta_next / M_2PI));
			theta_next = theta_next - loops * M_2PI;

			bool isFlipped = false;
			//�����ϼ���phi��ǰ�ƽ�pi
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

			//ת������������
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


//���pngͼƬ
void Drawer::output_png(const char* filename)
{
	//Encode the image
	unsigned error = lodepng::encode(filename, image, width, height);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}