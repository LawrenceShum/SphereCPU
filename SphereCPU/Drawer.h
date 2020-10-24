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
	vector<unsigned char> image;
	color white = { 255,255,255,255 };
	color red = { 255,0,0,255 };
	color green = { 0,128,0,255 };
	color yellow = { 255,255,0,255 };
	unsigned width;
	unsigned height;
	unsigned num_particles = 10000;
	vector<particle> particles;

	void set_point(vector<unsigned char>&, int, int, color);
	void initial(vector<unsigned char>&);

public :
	Drawer(unsigned,unsigned);

	~Drawer();

	vector<unsigned char> get_image();
};