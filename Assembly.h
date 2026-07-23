#pragma once
#include "Parsing.h"

class Assembly {
public:
	std::vector<int> distribution;
	int index = 0;
	int total_size = 0;

};

class Forming {
public:
	void distribution_forming(Assembly& assembly, Coord& coord, Material& material, Geometry& geometry, Cel_Data& cel, Rep_Data& rep, Sur_Data& sur);
	double evaluate_surface(double local_x, double local_y, double local_z, Sur_Data& sur) const;
	int determine_material(double local_x, double local_y, double local_z, int current_cell_id, const Geometry& geometry) const;
	void printing(Assembly& assembly, Material& material, Geometry& geometry);
};