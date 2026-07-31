#pragma once
#include "Parsing.h"
#include "lcg2.h"
#include <queue>

struct Neutron {
public:
	double x = 0.0, y = 0.0, z = 0.0;
	double u = 0.0, v = 0.0, w = 0.0;
	int group = 0;
	double weight = 0.0;

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
	void distance(lcg& rn, Neutron& neutron, Material& material, Geometry& geometry, Parsing& parsing, Coord& coord);
};


class Manager {
private:
	int NPS = 0;
	int total_cycles = 150;
	int inactive_cycles = 50;
	int current_NPS = 0.0;
	
	double tally_sum = 0.0;
	double k = 1.0;


	double generation_bank_count = 0.0;
	std::vector<Neutron> current_bank;
	std::vector<Neutron> next_bank;

	lcg rn = lcg(123456789ULL);
	Neutron neutron;
	Material material;
	Geometry geometry;
	Factory factory;
	Distance distance;
	Parsing parsing;
	Coord coord;
	Mat_Data mat;
	Cel_Data cel;
	Sur_Data sur;

public:
	void cycle();
	void iteration(int numNeutron);
	// Manager 클래스 (Parsing.h 또는 Simulator.h)
public:
	void set_data(const Material& mat, const Geometry& geo) {
		this->material = mat;
		this->geometry = geo;
	}
};
