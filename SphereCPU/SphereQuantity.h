#pragma once
#pragma once
#include "TotalInclude.h"

class SphereQuantity
{
private:
	int nPhi;
	int nTheta;
	float gridLen;
	float invGridLen;

	float* cpuBuffer;

public:
	SphereQuantity(size_t nPhi, size_t nTheta);

	~SphereQuantity();

	void swapGPUbuffer();
	void copyGPUtoCPU();
	void copyCPUtoGPU();

	size_t get_nPhi();
	size_t get_nTheta();

	float getCPUvalueAt(size_t x, size_t y);
	void setCPUvalueAt(size_t x, size_t y, float val);

	float* getGPUthisStep();
	float* getGPUnextStep();

	size_t getThisPitch();
	size_t getNextPitch();
};