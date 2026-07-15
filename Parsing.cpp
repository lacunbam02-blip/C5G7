#include "Parsing.h"

void Parsing::parsing(const std::string& filename, Material& material, Geometry& geometry) {
	std::ifstream input_file(filename);

	std::string word;
	while (input_file >> word) {
		if (word[0] == '#') continue;
		if (word == "Material") {
			Mat_Data mat;
			input_file >> mat.id;

			input_file >> word;
			if (word[0] == '#') mat.name = word.substr(1);     //check the role of substr later!

			for (int i = 0; i < 7; ++i) {
				double var1, var2, var3, var4, var5, var6;
				input_file >> var1 >> var2 >> var3 >> var4 >> var5 >> var6;
				mat.xs_t.push_back(var1);
				mat.xs_a.push_back(var2);
				mat.xs_c.push_back(var3);
				mat.xs_f.push_back(var4);
				mat.nu.push_back(var5);
				mat.chi.push_back(var6);
			}

			for (int i = 0; i < 49; ++i) {
				double var7;
				input_file >> var7;
				mat.xs_s.push_back(var7);
			}
			material.materials.push_back(mat);
		}
		if (word == "Surface") {
			Sur_Data sur;
			input_file >> sur.id;

			input_file >> sur.type;

			double var;
			input_file >> var;
			sur.coefficient.push_back(var);

			if (sur.type == "PX" && sur.coefficient.back() > 0) geometry.pitch_i = 2 * sur.coefficient.back();
			if (sur.type == "PY" && sur.coefficient.back() > 0) geometry.pitch_j = 2 * sur.coefficient.back();
			if (sur.type == "PZ" && sur.coefficient.back() > 0) geometry.pitch_k = 2 * sur.coefficient.back();

			geometry.surfaces.push_back(sur);
		}
		if (word == "Cell") {
			Cel_Data cel;
			input_file >> cel.id;

			std::string next_word;
			while (input_file >> next_word) {
				if (next_word[0] == 'C' || next_word == "size") {
					word = next_word;
					break;
				}
				cel.material_id = std::stoi(next_word);    // use (std::stoi) instead (static_cast)

				std::string condition_word;
				while (input_file >> condition_word) {
					if (condition_word == "/") break;
					cel.boundary_condition.push_back(std::stoi(condition_word));
				}
				geometry.cells.push_back(cel);
			}
		}
		if (word == "Size") {
			input_file >> geometry.size_i;
			input_file >> geometry.size_j;
			input_file >> geometry.size_k;
		}
		if (word == "Default_Cell") input_file >> geometry.default_cell;
		if (word == "Replace_Cell") {
			Rep_Data rep;
			input_file >> rep.id;
			
			std::string position_word;
			while (input_file >> position_word) {
				if (position_word[0] != '[') {
					word = position_word;
					break;
				}
				Coord coord;
				sscanf(position_word.c_str(), "[%d, %d]", &coord.i, &coord.j);

				coord.i = coord.i - 1;
				coord.j = coord.j - 1;
				coord.k = coord.k - 1;

				rep.replace_cell.push_back(coord);
			}
			geometry.replace_cells.push_back(rep);
		}
	}
}