#pragma once
#include <iostream>
#include <vector>


class lcg {

public:
	unsigned long long lamda;
	unsigned long long mu;
	unsigned long long ran;

	lcg(unsigned long long seed);
};

class Method {
public:
	double next(lcg& rn);
};