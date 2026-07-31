#include "../include/Parsing.h"

void Parsing::parsing(const std::string& filename, Material& material, Geometry& geometry) {
	std::ifstream input_file(filename);

	if (!input_file.is_open()) {
		std::cerr << "Error: '" << filename << "' 파일을 찾을 수 없습니다!" << std::endl;
		exit(1);
	}   //파일 열리는지 확인



	std::string word;

	bool skip_read = false; //check

	while (skip_read || input_file >> word) {
		skip_read = false; //check

		if (word[0] == '#') continue;
		else if (word == "Material") {   //material 정보 파싱
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
		else if (word == "Surface") {
			Sur_Data sur;
			input_file >> sur.id;

			input_file >> sur.type;

			double var;
			input_file >> var;
			sur.coefficient.push_back(var);

			if (sur.type == "CZ" && sur.coefficient.back() > 0) geometry.pitch_r = sur.coefficient.back();
			if (sur.type == "PX" && sur.coefficient.back() > 0) geometry.pitch_i = 2 * sur.coefficient.back();
			if (sur.type == "PY" && sur.coefficient.back() > 0) geometry.pitch_j = 2 * sur.coefficient.back();
			if (sur.type == "PZ" && sur.coefficient.back() > 0) geometry.pitch_k = 2 * sur.coefficient.back();

			geometry.surfaces.push_back(sur);
		}
		else if (word == "Cell") {
			Cel_Data cel;
			input_file >> cel.id;

			std::string next_word;
			while (input_file >> next_word) {
				if (next_word[0] == '#') continue;
				if (next_word == "Size" || next_word == "size" || next_word == "Core" || next_word == "Cell") {
					word = next_word;
					skip_read = true;  //check
					break;
				}

				Sub_Cell sub_cell;
				sub_cell.material_id = std::stoi(next_word);    // use (std::stoi) instead (static_cast)

				std::string condition_word;
				while (input_file >> condition_word) {
					if (condition_word == "/") break;
					sub_cell.boundary_condition.push_back(std::stoi(condition_word));
				}
				cel.sub_cells.push_back(sub_cell);
			}
			geometry.cells.push_back(cel);
		}
		else if (word == "Size") {
			input_file >> geometry.size_i;
			input_file >> geometry.size_j;
			input_file >> geometry.size_k;
		}

		else if (word == "Core") {
			geometry.total_size = geometry.size_i * geometry.size_j * geometry.size_k;

			std::vector<int> extra_data;
			extra_data.reserve(geometry.size_i * geometry.size_j);  //임시 한 층 벡터

			for (int j = 0; j < geometry.size_i * geometry.size_j; j++) {
				int var8;
				input_file >> var8;
				extra_data.push_back(var8);
			}

			/*
			* 	for (int k = 0; k < geometry.size_k; k++) {
				if (k < 17) {
					// 0~16 층 상단 데이터
					for (int j = 0; j < geometry.size_i * geometry.size_j; j++) {
						geometry.distribution.push_back(0);
					}
				}
				else {
					// 17~50 층 하단 데이터
					for (int j = 0; j < geometry.size_i * geometry.size_j; j++) {
						geometry.distribution.push_back(extra_data[j]);
					}
				}
			}
			*/
			// Parsing_2.cpp 의 Core 파싱 부분
			for (int k = 0; k < geometry.size_k; k++) {
				if (k < 153) {
					// 0~33 층: Z축 바닥 (반사 경계면, 노심 중심부) -> 연료 배치
					for (int j = 0; j < geometry.size_i * geometry.size_j; j++) {
						geometry.distribution.push_back(extra_data[j]);
					}
				}
				else {
					// 34~50 층: Z축 천장 (누설 경계면, 반사체 덮개) -> 감속재(0) 배치
					for (int j = 0; j < geometry.size_i * geometry.size_j; j++) {
						geometry.distribution.push_back(0);
					}
				}
			}

		
			
	
		}
	}
}

double Parsing::evaluate_surface(double local_x, double local_y, double local_z, Sur_Data& sur) const {
	double value = 0.0;

	if (sur.type == "CZ") {
		value = (local_x * local_x) + (local_y * local_y) - (sur.coefficient[0] * sur.coefficient[0]);
	}
	else if (sur.type == "PX") {
		value = local_x - sur.coefficient[0];
	}
	else if (sur.type == "PY") {
		value = local_y - sur.coefficient[0];
	}
	else if (sur.type == "PZ") {
		value = local_z - sur.coefficient[0];
	}
	if (std::abs(value) < 1e-9) {
		value = 0.0;
	}

	return value;
}

int Parsing::determine_material(double local_x, double local_y, double local_z, int current_cell_id, const Geometry& geometry) const {
	for (const auto& cel : geometry.cells) {
		if (cel.id != current_cell_id) {     // find 'cell id' in cells vector
			continue;
		}

		for (const auto& sub : cel.sub_cells) {
			bool is_inside = true;            // asumption that neutron is in material

			for (int bc : sub.boundary_condition) {
				int sur_id = std::abs(bc);
				int sense = (bc > 0) ? 1 : -1;

				Sur_Data target_sur;
				for (const auto& s : geometry.surfaces) {
					if (s.id == sur_id) {
						target_sur = s;
						break;
					}
				}

				double val = evaluate_surface(local_x, local_y, local_z, target_sur);

				if (sense > 0 && val < 0.0) {
					is_inside = false;
					break;
				}
				else if (sense < 0 && val > 0.0) {
					is_inside = false;
					break;
				}
			}

			if (is_inside) {
				return sub.material_id;
			}
		}
		if (!cel.sub_cells.empty()) {
			return cel.sub_cells.back().material_id;
		}
		else {
			return -1;
		}  //check
	}
	return current_cell_id;  //check
}

//if value is 0?

void Parsing::printing(Material& material, Geometry& geometry) {
	std::cout << "==================================================\n";
	std::cout << "                  MATERIAL DATA                   \n";
	std::cout << "==================================================\n";
	for (const auto& mat : material.materials) {
		std::cout << "Material ID: " << mat.id << " (" << mat.name << ")\n";
		std::cout << "  Group |     xs_t     |     xs_a     |     xs_c     |     xs_f     |      nu      |     chi      \n";
		std::cout << "  ------|--------------|--------------|--------------|--------------|--------------|--------------\n";
		for (size_t g = 0; g < mat.xs_t.size(); ++g) {
			std::cout << "    " << g << "   | " << mat.xs_t[g] << " | " << mat.xs_a[g] << " | " << mat.xs_c[g] << " | " << mat.xs_f[g] << " | " << mat.nu[g] << " | " << mat.chi[g] << "\n";
		}
		std::cout << "\n  Scattering Matrix (xs_s):\n  ";
		for (size_t i = 0; i < mat.xs_s.size(); ++i) {
			std::cout << mat.xs_s[i] << " ";
			if ((i + 1) % 7 == 0) std::cout << "\n  ";
		}
		std::cout << "\n--------------------------------------------------\n";
	}

	std::cout << "\n==================================================\n";
	std::cout << "                  SURFACE DATA                    \n";
	std::cout << "==================================================\n";
	for (const auto& sur : geometry.surfaces) {
		std::cout << "Surface ID: " << sur.id << " | Type: " << sur.type << " | Coeff: ";
		for (double c : sur.coefficient) {
			std::cout << c << " ";
		}
		std::cout << "\n";
	}

	std::cout << "\n==================================================\n";
	std::cout << "                    CELL DATA                     \n";
	std::cout << "==================================================\n";
	for (const auto& cel : geometry.cells) {
		std::cout << "Cell ID: " << cel.id << "\n";
		for (const auto& sub : cel.sub_cells) {
			std::cout << "  Material ID: " << sub.material_id << " | BC: ";
			for (int bc : sub.boundary_condition) {
				std::cout << (bc > 0 ? "+" : "") << bc << " ";
			}
			std::cout << "/\n";
		}
	}

	std::cout << "\n==================================================\n";
	std::cout << "                  GEOMETRY SIZE                   \n";
	std::cout << "==================================================\n";
	std::cout << "Size (i, j, k): " << geometry.size_i << ", " << geometry.size_j << ", " << geometry.size_k << "\n";
	std::cout << "Pitch (i, j, k, r): " << geometry.pitch_i << ", " << geometry.pitch_j << ", " << geometry.pitch_k << ", " << geometry.pitch_r << "\n";

	std::cout << "\n=====================================================================================================\n";
	std::cout << "                                      CELL DISTRIBUTION (K = 0)                                     \n";
	std::cout << "=====================================================================================================\n";

	int layer_k = 17;
	for (int i = 0; i < geometry.size_i; ++i) {
		for (int j = 0; j < geometry.size_j; ++j) {
			int idx = layer_k * (geometry.size_i * geometry.size_j) + i * geometry.size_j + j;
			std::cout << geometry.distribution[idx] << " ";
		}
		std::cout << "\n";
	}
	std::cout << "=====================================================================================================\n";
}
