#pragma once
#include "Parsing.h"

class Assembly {
public:
	std::vector<int> distribution;
	int index = 0;
	int total_size = 0;

	double cen_x= 0.0;
	double cen_y = 0.0;
	double cen_z = 0.0;
};

class Forming {
public:
	void distribution_forming(Assembly& assembly, Coord& coord, Material& material, Geometry& geometry, Cel_Data& cel, Rep_Data rep, Sur_Data sur);
	void coordinate_forming(Assembly& assembly, Coord& coord, Material& material, Geometry& geometry, Cel_Data& cel, Rep_Data rep, Sur_Data sur);
	double evaluate_surface(double lcoal_x, double local_y, double local_z, const Sur_Data& sur) const;
};