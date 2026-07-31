#include "../include/lcg2.h"
#include "../include/Parsing.h"
#include "../include/Simulator.h"
#include <iostream>
#include <filesystem>

int main() {

	Material material;
	Geometry geometry;
	Parsing parsing;

	parsing.parsing("../input/C5G7_Mat_Input.txt", material, geometry);
	parsing.printing(material, geometry);

	Manager manager;
	manager.set_data(material, geometry);
	manager.iteration(100000);
	return 0;
}

