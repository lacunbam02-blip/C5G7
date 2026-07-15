#include "Parsing.h"
#include "Assembly.h"

int main() {
	Mat_Data mat_data;
	Material material;
	Geometry geometry;
	Cel_Data cel;
	Sur_Data sur;
	Rep_Data rep;
	Coord coord;

	Parsing parsing;

	parsing.parsing("2D7G_Mat_Input.txt", material, geometry);

	Assembly assembly;
	Forming forming;

	forming.distribution_forming(assembly, coord, material, geometry, cel, rep, sur);
}

