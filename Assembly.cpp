#include "Assembly.h"

void Forming::distribution_forming(Assembly& assembly, Coord& coord, Material& material, Geometry& geometry, Cel_Data& cel, Rep_Data rep, Sur_Data sur) {
	assembly.index = coord.k * (geometry.size_i * geometry.size_j) + coord.j * geometry.size_i + coord.i;

	assembly.total_size = geometry.size_i * geometry.size_j * geometry.size_k;
	assembly.distribution.assign(assembly.total_size, geometry.default_cell);




	const int center_i = 8;
	const int center_j = 8;

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
	cel.id = assembly.distribution[assembly.index];
}

void Forming::coordinate_forming(Assembly& assembly, Coord& coord, Material& material, Geometry& geometry, Cel_Data& cel, Rep_Data rep, Sur_Data sur) {
	assembly.cen_z = geometry.pitch_k * (coord.k - ((assembly.total_size / 2) / (geometry.size_i * geometry.size_j)));
	assembly.cen_y = geometry.pitch_j * (coord.j - (((assembly.total_size / 2) % (geometry.size_i * geometry.size_j)))/ geometry.size_i);
	assembly.cen_x = geometry.pitch_i * (coord.i - (((assembly.total_size / 2) % (geometry.size_i * geometry.size_j))) % geometry.size_i);
	// conversion of coord(i, j, k) to center Cartesian(x, y, z)
	d
}