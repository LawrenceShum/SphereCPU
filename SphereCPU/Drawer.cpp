#include "Drawer.h"

using namespace std;

Drawer::Drawer(unsigned width, unsigned height) :
	width(width), height(height)
{
	image.resize(width * height * 4);
	particles.resize(num_particles);
	color_old.resize(width*height);
	color_new.resize(width*height);
	initial_particles();
	initial(image);
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
			if (y < (height / 8))
				color_old[y*width + x] = red;
			//set_point(image, x, y, red);
			else if (y < (height / 4))
				color_old[y*width + x] = white;
				//set_point(image, x, y, white);
			else if (y < (height * 3 / 8))
				color_old[y*width + x] = red;
				//set_point(image, x, y, red);
			else if (y < (height / 2))
				color_old[y*width + x] = white;
				//set_point(image, x, y, white);
			else if (y < (height * 5 / 8))
				color_old[y*width + x] = red;
				//set_point(image, x, y, red);
			else if (y < (height * 3 / 4))
				color_old[y*width + x] = white;
				//set_point(image, x, y, white);
			else if (y < (height * 7 / 8))
				color_old[y*width + x] = red;
				//set_point(image, x, y, red);
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
			/*
			if (y < (height / 8))
				set_point(image, x, y, red);
			else if (y < (height / 4))
				set_point(image, x, y, white);
			else if (y < (height * 3 / 8))
				set_point(image, x, y, red);
			else if (y < (height / 2))
				set_point(image, x, y, white);
			else if (y < (height * 5 / 8))
				set_point(image, x, y, red);
			else if (y < (height * 3 / 4))
				set_point(image, x, y, white);
			else if (y < (height * 7 / 8))
				set_point(image, x, y, red);
			else
				set_point(image, x, y, white);
				*/
		}
	}
}

void Drawer::initial_particles()
{
	for (unsigned long i = 0; i < num_particles; i++)
	{
		particles[i].x = (0 + i)%width;
		particles[i].y = 100 + static_cast<int>(i/width);
	}
}

//计算粒子运动
void Drawer::calculate_particles()
{
	for (unsigned long i = 0; i < num_particles; i++)
	{
		//获取粒子当前所在位置的速度，通过插值的方法
		//并计算下一时间步，粒子的位置
	}
}

//画出粒子
void Drawer::draw_particles()
{		

}

//颜色对流
void Drawer::color_advect()
{

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