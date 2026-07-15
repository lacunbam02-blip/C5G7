#pragma once
#include "Parsing.h"

class Assembly {
public:
	std::vector<int> distribution;
	int index = 0;
	int total_size = 0;

	double cen_x, cen_y, cen_z = 0.0;
};

class Forming {
public:
	void distribution_forming(Assembly& assembly, Coord& coord, Material& material, Geometry& geometry, Cel_Data& cel, Rep_Data rep, Sur_Data sur);
	void coordinate_forming(Assembly& assembly, Coord& coord, Material& material, Geometry& geometry, Cel_Data& cel, Rep_Data rep, Sur_Data sur);
};