#include "Assembly.h"
#include "Parsing.h"
#include "lcg2.h"
#include "2D7G.h"

void Factory::ini_pos(lcg& rn, Neutron& neutron, Assembly& assembly, Geometry& geometry, Coord& coord) {
	Method method;
	neutron.x = method.next(rn) * geometry.pitch_i * geometry.size_i - geometry.pitch_i * geometry.size_i / 2;
	neutron.y = method.next(rn) * geometry.pitch_j * geometry.size_j - geometry.pitch_j * geometry.size_j / 2;
	neutron.z = method.next(rn) * geometry.pitch_k * geometry.size_k - geometry.pitch_k * geometry.size_k / 2;

	coord.i = (neutron.x + geometry.pitch_i * geometry.size_i / 2 ) / geometry.pitch_i;
	coord.j = (neutron.y + geometry.pitch_j * geometry.size_j / 2 ) / geometry.pitch_j;
	coord.k = (neutron.z + geometry.pitch_k * geometry.size_k / 2 ) / geometry.pitch_k;

	neutron.center_x = geometry.pitch_i * (0.5 + coord.i);
	neutron.center_y = geometry.pitch_j * (0.5 + coord.j);
	neutron.center_z = geometry.pitch_k * (0.5 + coord.k);

	neutron.local_x = neutron.x - neutron.center_x;
	neutron.local_y = neutron.y - neutron.center_y;
	neutron.local_z = neutron.z - neutron.center_z;
}

void Factory::ini_dir(lcg& rn, Neutron& neutron) {
	Method method;
	double x, y, z, r_2 = 0.0;
	for (int i = 0; i < 40000; ++i) {
		x = method.next(rn) * 2 - 1;
		y = method.next(rn) * 2 - 1;
		z = method.next(rn) * 2 - 1;
		r_2 = x * x + y * y + z * z;
		if (r_2 <= 1 && r_2 != 0) break;
	}
	double r = sqrt(r_2);
	neutron.u = (x / r);
	neutron.v = (y / r);
	neutron.w = (z / r);
}



void Distance::distance(lcg& rn, Neutron& neutron, Material& material, Geometry& geometry, Assembly& assembly) {
	Method method;
	neutron.DTC = -(log(method.next(rn)) / material.materials[geometry.cells[neutron.current_cel_id].material_id].xs_t[neutron.group]);

	// current location (i, j, k) or (x, y, z) -> use index and cel.id
	if (geometry.cells[neutron.current_cel_id].material_id == 0 || geometry.cells[neutron.current_cel_id].material_id == 2 || geometry.cells[neutron.current_cel_id].material_id == 3) {
		double dts_w = 0.0;
		double dts_r = 0.0;
		if (neutron.w > 0) dts_w = (geometry.pitch_k / 2 - neutron.local_z) / neutron.w;
		else if (neutron.w < 0) dts_w = (-1 * geometry.pitch_k / 2 - neutron.local_z) / neutron.w;
		else dts_w = 1e30;

		double a = -neutron.u * neutron.local_x - neutron.v * neutron.local_y;
		double b = (neutron.u * neutron.local_x + neutron.v * neutron.local_y) * (neutron.u * neutron.local_x + neutron.v * neutron.local_y) - (neutron.local_x * neutron.local_x + neutron.local_y * neutron.local_y - geometry.pitch_r * geometry.pitch_r) * (neutron.u * neutron.u + neutron.v * neutron.v);
		double c = neutron.u * neutron.u + neutron.v * neutron.v;
		if (b > 0 && c != 0) dts_r = (a + sqrt(b)) / c;
		else dts_r = 1e30;

		neutron.DTS = std::min({ dts_w, dts_r });
	}
	if (geometry.cells[neutron.current_cel_id].material_id == 1) {
		double dts_u = 0.0;
		double dts_v = 0.0;
		double dts_w = 0.0;
		double dts_r = 0.0;
		if (neutron.u > 0) dts_u = (geometry.pitch_i / 2 - neutron.local_x) / neutron.u;
		else if (neutron.u < 0) dts_u = (-1 * geometry.pitch_i / 2 - neutron.local_x) / neutron.u;
		else dts_u = 1e30;

		if (neutron.v > 0) dts_v = (geometry.pitch_j / 2 - neutron.local_y) / neutron.v;
		else if (neutron.v < 0) dts_v = (-1 * geometry.pitch_j / 2 - neutron.local_y) / neutron.v;
		else dts_v = 1e30;

		if (neutron.w > 0) dts_w = (geometry.pitch_k / 2 - neutron.local_z) / neutron.w;
		else if (neutron.w < 0) dts_w = (-1 * geometry.pitch_k / 2 - neutron.local_z) / neutron.w;
		else dts_w = 1e30;

		double a = -neutron.u * neutron.local_x - neutron.v * neutron.local_y;
		double b = (neutron.u * neutron.local_x + neutron.v * neutron.local_y) * (neutron.u * neutron.local_x + neutron.v * neutron.local_y) - (neutron.local_x * neutron.local_x + neutron.local_y * neutron.local_y - geometry.pitch_r * geometry.pitch_r) * (neutron.u * neutron.u + neutron.v * neutron.v);
		double c = neutron.u * neutron.u + neutron.v * neutron.v;
		if (b > 0 && c != 0) dts_r = (a - sqrt(b)) / c;
		else dts_r = 1e30;

		neutron.DTS = std::min({ dts_u, dts_v, dts_w, dts_r });
	}
}

void Manager::initialize() {
	Parsing parsing;

	parsing.parsing("2D7G_Mat_Input.txt", material, geometry);

	Forming forming;

	forming.distribution_forming(assembly, coord, material, geometry, cel, rep, sur);
}