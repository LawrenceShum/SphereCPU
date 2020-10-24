#pragma once
#include <iostream>
#include "lodepng.h"
#include <vector>

using namespace std;

struct color {
	int R;
	int G;
	int B;
	int A;
};

struct particle {
	double x;
	double y;
};

class Drawer 
{
private:
	//画布
	vector<unsigned char> image;
	//颜色定义
	color white = { 255,255,255,255 };
	color red = { 255,0,0,255 };
	color green = { 0,128,0,255 };
	color yellow = { 255,255,0,255 };
	//图像长宽
	unsigned width;
	unsigned height;
	//粒子的数量
	unsigned num_particles = 10000;
	vector<particle> particles;
	void initial_particles();
	
	//在png上画点
	void set_point(vector<unsigned char>&, int, int, color);
	void initial(vector<unsigned char>&);

public :
	Drawer(unsigned,unsigned);

	~Drawer();

	vector<unsigned char> get_image();
	void output_png(const char*, std::vector<unsigned char>&, unsigned, unsigned);
};