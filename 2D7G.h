#pragma once
#include "Assembly.h"
#include "Parsing.h"
#include "lcg2.h"
#include <queue>

struct Neutron {
public:
	double x = 0.0, y = 0.0, z = 0.0;
	double u = 0.0, v = 0.0, w = 0.0;
	int group = 0;
	double weight = 1.0;

	double center_x = 0.0, center_y = 0.0, center_z = 0.0;
	double local_x = 0.0, local_y = 0.0, local_z = 0.0;
	int current_index = 0;
	int current_cel_id = 0;
	int current_material = 0;

	double DTC = 0.0;
	double DTS = 0.0;
};

class Factory {
public:
	void ini_pos(lcg& rn, Neutron& neutron, Geometry& geometry);
	void ini_dir(lcg& rn, Neutron& neutron);
};

class Distance {
public:
	void distance(lcg& rn, Neutron& neutron, Material& material, Geometry& geometry, Assembly& assembly, Forming& forming, Coord& coord);
};


class Manager {
private:
	int NPS = 100000;
	int total_cycles = 300;
	int current_NPS = 0.0;
	double tally_sum = 0.0;
	double k = 1.31115;

	double generation_bank_count = 0.0;
	std::queue<Neutron> current_bank;
	std::queue<Neutron> next_bank;

	lcg rn = lcg(123456789ULL);
	Neutron neutron;
	Material material;
	Geometry geometry;
	Factory factory;
	Distance distance;
	Assembly assembly;
	Forming forming;
	Coord coord;
	Mat_Data mat;
	Cel_Data cel;
	Rep_Data rep;
	Sur_Data sur;

public:
	void initialize(int numNeutron);
	void cycle();
	void iteration();
};