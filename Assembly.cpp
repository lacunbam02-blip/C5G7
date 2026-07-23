#include "Assembly.h"

void Forming::distribution_forming(Assembly& assembly, Coord& coord, Material& material, Geometry& geometry, Cel_Data& cel, Rep_Data& rep, Sur_Data& sur) {
	assembly.total_size = geometry.size_i * geometry.size_j * geometry.size_k;
	assembly.distribution.assign(assembly.total_size, geometry.default_cell);

	const int center_i = static_cast<int>(geometry.size_i / 2.0);
	const int center_j = static_cast<int>(geometry.size_j / 2.0);

	for (const auto& rep : geometry.replace_cells) {
		int change_cell_id = rep.id;
		for (const auto& change_coord : rep.replace_cell) {
			int di = change_coord.i - center_i;
			int dj = change_coord.j - center_j;
			
			std::pair <int, int> symmetry[8] = {
				{center_i + di, center_j + dj},
				{center_i - di, center_j + dj},
				{center_i + di, center_j - dj},
				{center_i - di, center_j - dj},
				{center_i + dj, center_j + di},
				{center_i - dj, center_j + di},
				{center_i + dj, center_j - di},
				{center_i - dj, center_j - di}
			};

			for (int i = 0; i < geometry.size_k; ++i) {
				for (int j = 0; j < 8; ++j) {
					int sym_x = symmetry[j].first;
					int sym_y = symmetry[j].second;
					if (sym_x < 0 || sym_x >= geometry.size_i || sym_y < 0 || sym_y >= geometry.size_j) {
						continue;
					}
					int idx = i * (geometry.size_i * geometry.size_j) + sym_x * geometry.size_i + sym_y;
					assembly.distribution[idx] = change_cell_id;
				}
			}
		}
	}
}

double Forming::evaluate_surface(double local_x, double local_y, double local_z, Sur_Data& sur) const {
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

int Forming::determine_material(double local_x, double local_y, double local_z, int current_cell_id, const Geometry& geometry) const {
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

void Forming::printing(Assembly& assembly, Material& material, Geometry& geometry) {
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
	std::cout << "Default Cell: " << geometry.default_cell << "\n";

	std::cout << "\n==================================================\n";
	std::cout << "             CELL DISTRIBUTION (K = 0)            \n";
	std::cout << "==================================================\n";

	int layer_k = 0;
	for (int i = 0; i < geometry.size_i; ++i) {
		for (int j = 0; j < geometry.size_j; ++j) {
			int idx = layer_k * (geometry.size_i * geometry.size_j) + i * geometry.size_j + j;
			std::cout << assembly.distribution[idx];
		}
		std::cout << "\n";
	}
	std::cout << "==================================================\n";
}