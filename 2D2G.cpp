#include "Assembly.h"
#include "Parsing.h"
#include "lcg2.h"
#include "2D7G.h"
#include <algorithm>

void Factory::ini_pos(lcg& rn, Neutron& neutron, Geometry& geometry) {
	Method method;
	neutron.x = method.next(rn) * geometry.pitch_i * geometry.size_i - geometry.pitch_i * geometry.size_i / 2;
	neutron.y = method.next(rn) * geometry.pitch_j * geometry.size_j - geometry.pitch_j * geometry.size_j / 2;
	neutron.z = method.next(rn) * geometry.pitch_k * geometry.size_k - geometry.pitch_k * geometry.size_k / 2;
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



void Distance::distance(lcg& rn, Neutron& neutron, Material& material, Geometry& geometry, Assembly& assembly, Forming& forming, Coord& coord) {
	
	int calc_i = static_cast<int>((neutron.x + geometry.pitch_i * geometry.size_i / 2) / geometry.pitch_i);
	int calc_j = static_cast<int>((neutron.y + geometry.pitch_j * geometry.size_j / 2) / geometry.pitch_j);
	int calc_k = static_cast<int>((neutron.z + geometry.pitch_k * geometry.size_k / 2) / geometry.pitch_k);

	if (calc_i >= geometry.size_i) calc_i = geometry.size_i - 1;
	else if (calc_i < 0) calc_i = 0;

	if (calc_j >= geometry.size_j) calc_j = geometry.size_j - 1;
	else if (calc_j < 0) calc_j = 0;

	if (calc_k >= geometry.size_k) calc_k = geometry.size_k - 1;
	else if (calc_k < 0) calc_k = 0;

	coord.i = std::max(0, std::min(calc_i, geometry.size_i - 1));
	coord.j = std::max(0, std::min(calc_j, geometry.size_j - 1));
	coord.k = std::max(0, std::min(calc_k, geometry.size_k - 1));   //check

	neutron.center_x = geometry.pitch_i * (0.5 + coord.i);
	neutron.center_y = geometry.pitch_j * (0.5 + coord.j);
	neutron.center_z = geometry.pitch_k * (0.5 + coord.k);

	neutron.local_x = neutron.x - neutron.center_x;
	neutron.local_y = neutron.y - neutron.center_y;
	neutron.local_z = neutron.z - neutron.center_z;

	neutron.current_index = coord.k * (geometry.size_i * geometry.size_j) + coord.j * geometry.size_i + coord.i;
	neutron.current_cel_id = assembly.distribution[neutron.current_index];
	
	
	
	
	Method method;

	neutron.current_material = forming.determine_material(neutron.local_x, neutron.local_y, neutron.local_z, neutron.current_cel_id, geometry);

	if (neutron.current_material == -1) {
		return;
	}

	neutron.DTC = -(log(method.next(rn)) / material.materials[neutron.current_material].xs_t[neutron.group]);

	// current location (i, j, k) or (x, y, z) -> use index and cel.id
	if (neutron.current_material == 0 || neutron.current_material == 2 || neutron.current_material == 3) {
		double dts_w = 0.0;
		double dts_r = 0.0;
		if (neutron.w > 0) dts_w = (geometry.pitch_k / 2 - neutron.local_z) / neutron.w;
		else if (neutron.w < 0) dts_w = (-1 * geometry.pitch_k / 2 - neutron.local_z) / neutron.w;
		else dts_w = 1e30;

		if (dts_w < 1e-9) dts_w = 1e30; //check


		double a = -neutron.u * neutron.local_x - neutron.v * neutron.local_y;
		double b = (neutron.u * neutron.local_x + neutron.v * neutron.local_y) * (neutron.u * neutron.local_x + neutron.v * neutron.local_y) - (neutron.local_x * neutron.local_x + neutron.local_y * neutron.local_y - geometry.pitch_r * geometry.pitch_r) * (neutron.u * neutron.u + neutron.v * neutron.v);
		double c = neutron.u * neutron.u + neutron.v * neutron.v;
		if (b > 0 && c != 0) dts_r = (a + sqrt(b)) / c;
		else dts_r = 1e30;

		neutron.DTS = std::min({ dts_w, dts_r });
	}
	if (neutron.current_material == 1) {
		double dts_u = 0.0;
		double dts_v = 0.0;
		double dts_w = 0.0;
		double dts_r = 0.0;


		if (dts_u < 1e-9) dts_u = 1e30; //check
		if (dts_v < 1e-9) dts_v = 1e30;
		if (dts_w < 1e-9) dts_w = 1e30;


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

void Manager::cycle() {
	Method method;


	for (int i = 0; i < current_NPS; ++i) {
		if (current_bank.empty()) {
			factory.ini_pos(rn, neutron, geometry);
			neutron.group = 0;
		}
		else {
			neutron = current_bank.front();
			current_bank.pop();
		}
		factory.ini_dir(rn, neutron);
		bool loop_active = true;

		while (loop_active) {
			
			distance.distance(rn, neutron, material, geometry, assembly, forming, coord);

			if (neutron.current_material == -1) {
				loop_active = false;
				continue;
			}

			if (neutron.DTC < neutron.DTS) {
				double r = method.next(rn) * material.materials[neutron.current_material].xs_t[neutron.group];
				double sum_scattering = std::accumulate(material.materials[neutron.current_material].xs_s.begin() + neutron.group * 7, material.materials[neutron.current_material].xs_s.begin() + (neutron.group * 7 + 7), 0.0);

				if (r < sum_scattering) {
					double accumulated_xs = 0.0;
					for (int next_g = 0; next_g < 7; ++next_g) {
						accumulated_xs += material.materials[neutron.current_material].xs_s[neutron.group * 7 + next_g];

						if (r < accumulated_xs) {
							neutron.group = next_g;
							break;
						}
					}
					neutron.x = neutron.x + neutron.u * neutron.DTC;
					neutron.y = neutron.y + neutron.v * neutron.DTC;
					neutron.z = neutron.z + neutron.w * neutron.DTC;

					factory.ini_dir(rn, neutron);
				}

				else if (r < sum_scattering + material.materials[neutron.current_material].xs_c[neutron.group]) {
					loop_active = false;
				}

				else  {
					neutron.x = neutron.x + neutron.u * neutron.DTC;
					neutron.y = neutron.y + neutron.v * neutron.DTC;
					neutron.z = neutron.z + neutron.w * neutron.DTC;

					tally_sum += material.materials[neutron.current_material].nu[neutron.group] * material.materials[neutron.current_material].xs_f[neutron.group] / material.materials[neutron.current_material].xs_t[neutron.group];
					int nu_generated = static_cast<int>(material.materials[neutron.current_material].nu[neutron.group] *
						material.materials[neutron.current_material].xs_f[neutron.group] / material.materials[neutron.current_material].xs_t[neutron.group]
						/ k + method.next(rn));
					
					generation_bank_count += nu_generated;

					double accumulated_chi = 0.0;
					for (int g = 0; g < nu_generated; ++g) {
						double accumulated_chi = 0.0;

						for (int next_g = 0; next_g < 7; ++next_g) {
							accumulated_chi += material.materials[neutron.current_material].chi[next_g];

							if (method.next(rn) < accumulated_chi) {
								neutron.group = next_g;
								break;
							}
						}
						next_bank.push(neutron);
					}
					loop_active = false;
				}
			}

			else {
				neutron.x = neutron.x + neutron.u * neutron.DTS;
				neutron.y = neutron.y + neutron.v * neutron.DTS;
				neutron.z = neutron.z + neutron.w * neutron.DTS;

				if (neutron.x >= geometry.pitch_i * geometry.size_i / 2 || neutron.x <= -geometry.pitch_i * geometry.size_i / 2) {
					neutron.u *= -1.0;
				}

				if (neutron.y >= geometry.pitch_j * geometry.size_j / 2 || neutron.y <= -geometry.pitch_j * geometry.size_j / 2) {
					neutron.v *= -1.0;
				}

				if (neutron.z >= geometry.pitch_k * geometry.size_k / 2 || neutron.z <= -geometry.pitch_k * geometry.size_k / 2) {
					neutron.w *= -1.0;
				}
			}
			
		}

	}
}

void Manager::iteration() {

	double preceed_k = 0.0;
	int inactive_cycles = 50;
	double active_k_sum = 0.0;
	double active_k_sq_sum = 0.0;
	int active_count = 0;

	current_NPS = NPS;

	for (int i = 0; i < total_cycles; ++i) {
		generation_bank_count = 0.0;
		tally_sum = 0.0;

		cycle();

		preceed_k = k;
		k = tally_sum / current_NPS;

		current_bank = std::move(next_bank);
		current_NPS = current_bank.size();

		if (i >= inactive_cycles) {
			active_k_sum += k;
			active_k_sq_sum += (k * k);
			active_count++;
			std::cout << "Active Cycle " << i + 1 << " k_eff: " << k << "\n";
		}
		else {
			std::cout << "Inactive Cycle " << i + 1 << " k_eff: " << k << "\n";
		}
	}
	if (active_count > 0) {
		double avg_k = active_k_sum / active_count;
		double variance = (active_k_sq_sum / active_count) - (avg_k * avg_k);
		double std_dev = std::sqrt(variance / (active_count - 1));

		std::cout << avg_k << " ± " << std_dev << "\n";
	}
}