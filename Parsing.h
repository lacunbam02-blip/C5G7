#pragma once
#include <iostream>
#include <vector>
#include <string>
#include <fstream>
#include <numeric>

struct Mat_Data {
public:
	int id = 0;
	std::string name;
	std::vector<double> xs_t;
	std::vector<double> xs_a;
	std::vector<double> xs_c;
	std::vector<double> xs_f;
	std::vector<double> nu;
	std::vector<double> chi;
	std::vector<double> xs_s;
};

class Material {
public:
	std::vector<Mat_Data> materials;
};

//

struct Sur_Data {
	int id = 0;
	std::string type;
	std::vector<double> coefficient;
};

struct Sub_Cell {
	int material_id = 0;
	std::vector<int> boundary_condition;
};

struct Cel_Data {
	int id = 0;
	std::vector<Sub_Cell> sub_cells;
};

struct Coord {
	int i = 0;
	int j = 0;
	int k = 0;
};

struct Rep_Data {
	int id = 0;
	std::vector<Coord> replace_cell;
};


class Geometry {
public:
	std::vector<Sur_Data> surfaces;
	std::vector<Cel_Data> cells;

	int default_cell = 0;
	int size_i = 0;
	int size_j = 0;
	int size_k = 0;

	double pitch_r = 0.0;
	double pitch_i = 0.0;
	double pitch_j = 0.0;
	double pitch_k = 0.0;

	std::vector<Rep_Data> replace_cells;
};

//

class Parsing {
public:
	void parsing(const std::string& filename, Material& material, Geometry& geometry);
};
