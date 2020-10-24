#include "Drawer.h"

using namespace std;

Drawer::Drawer(unsigned width, unsigned height) :
	width(width), height(height)
{
	image.resize(width * height * 4);
	particles.resize(num_particles);
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
			//将画布变成黑色
			image[4 * y*width + 4 * x + 3] = 255;
		}
	}
}

void Drawer::initial_particles()
{
	for (unsigned i = 0; i < num_particles; i++)
	{
		particles[i].x = 0 + i;
		particles[i].y = 0 + i;
	}
}

void Drawer::set_point(vector<unsigned char>& image, int x_i, int y_i, color color)
{
	int x = x_i;
	int y = height - y_i;
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
void Drawer::output_png(const char* filename, std::vector<unsigned char>& image, unsigned width, unsigned height)
{
	//Encode the image
	unsigned error = lodepng::encode(filename, image, width, height);

	//if there's an error, display it
	if (error) std::cout << "encoder error " << error << ": " << lodepng_error_text(error) << std::endl;
}