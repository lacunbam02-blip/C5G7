#include "lcg2.h"
#include "Parsing.h"
#include "Assembly.h"
#include "2D7G.h"

int main() {
	Manager manager;

	Mat_Data mat_data;
	Material material;
	Geometry geometry;

	Parsing parsing;

	parsing.parsing("2D7G_Mat_Input.txt", material, geometry);

	Forming forming;

	forming.distribution_forming(assembly, coord, material, geometry, cel, rep, sur);
	forming.coordinate_forming(assembly, coord, material, geometry, cel, rep, sur);
}

