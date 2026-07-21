#include "Assembly.h"

void Forming::distribution_forming(Assembly& assembly, Coord& coord, Material& material, Geometry& geometry, Cel_Data& cel, Rep_Data rep, Sur_Data sur) {
	assembly.total_size = geometry.size_i * geometry.size_j * geometry.size_k;
	assembly.distribution.assign(assembly.total_size, geometry.default_cell);

	const int center_i = geometry.size_i/2;
	const int center_j = geometry.size_j/2;

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
		return cel.sub_cells.back().material_id;  //check
	}
	return -1;
}

//if value is 0?