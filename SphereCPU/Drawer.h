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
	vector<color> color_old;
	vector<color> color_new;
	//颜色定义
	color white = { 255,255,255,255 };
	color red = { 255,0,0,255 };
	color green = { 0,128,0,255 };
	color yellow = { 255,255,0,255 };
	//图像长宽
	unsigned width;
	unsigned height;
	//粒子的数量
	unsigned long num_particles = 200*800;
	vector<particle> particles;
	void initial_particles();
	
	//在png上画点
	void set_point(vector<unsigned char>&, int, int, color);
	void initial(vector<unsigned char>&);
	
	//计算粒子的运动
	void calculate_particles();
	//画出粒子
	void draw_particles();
	//颜色的对流
	void color_advect();

public :
	Drawer(unsigned,unsigned);

	~Drawer();

	vector<unsigned char> get_image();
	void output_png(const char*);

	friend class SphereSolver;
};