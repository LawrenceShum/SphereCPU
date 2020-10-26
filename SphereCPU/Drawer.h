#pragma once
#include <iostream>
#include "lodepng.h"
#include "SphereSolver.h"
#include <vector>
#include "Sphere.h"

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
	//����
	vector<unsigned char> image;
	vector<color> color_old;
	vector<color> color_new;
	//�������ָ��
	SphereSolver* Psolver;
	//��ɫ����
	color white = { 255,255,255,255 };
	color red = { 255,0,0,255 };
	color green = { 0,128,0,255 };
	color yellow = { 255,255,0,255 };
	color black = { 0,0,0,255 };
	//ͼ�񳤿�
	unsigned width;
	unsigned height;
	//�����������صĻ���
	double Angle2Pixel;
	//���ӵ�����
	unsigned long num_particles = 100*100;

	vector<particle> particles;

	void initial_particles();

public :
	Drawer(unsigned,unsigned,SphereSolver*);

	~Drawer();

	vector<unsigned char> get_image();

	//��png�ϻ���
	void set_point(vector<unsigned char>&, int, int, color);
	void initial(vector<unsigned char>&);

	//�������ӵ��˶�
	void calculate_particles();
	void draw_particles();
	//������
	//����
	void draw();
	//��ɫ�Ķ���
	void color_advect();
	//������ɫ
	void update_color();

	void output_png(const char*);

	friend class SphereSolver;
};