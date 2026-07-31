#include "../include/Parsing.h"
#include "../include/lcg2.h"
#include "../include/Simulator.h"
#include <algorithm>

void Factory::ini_pos(lcg& rn, Neutron& neutron, Geometry& geometry) {
	Method method;
	neutron.x = method.random_number_generator(rn) * geometry.pitch_i * geometry.size_i - geometry.pitch_i * geometry.size_i / 2.0;
	neutron.y = method.random_number_generator(rn) * geometry.pitch_j * geometry.size_j - geometry.pitch_j * geometry.size_j / 2.0;
	neutron.z = method.random_number_generator(rn) * geometry.pitch_k * geometry.size_k - geometry.pitch_k * geometry.size_k / 2.0;
}

void Factory::ini_dir(lcg& rn, Neutron& neutron) {
	Method method;
	double x, y, z, r_2 = 0.0;
	for (int i = 0; i < 40000; ++i) {
		x = method.random_number_generator(rn) * 2 - 1;
		y = method.random_number_generator(rn) * 2 - 1;
		z = method.random_number_generator(rn) * 2 - 1;
		r_2 = x * x + y * y + z * z;
		if (r_2 <= 1 && r_2 != 0) break;
	}
	double r = sqrt(r_2);
	neutron.u = (x / r);
	neutron.v = (y / r);
	neutron.w = (z / r);
}



void Distance::distance(lcg& rn, Neutron& neutron, Material& material, Geometry& geometry, Parsing& parsing, Coord& coord) {
	
	double logical_x = neutron.x;
	double logical_y = neutron.y;
	double logical_z = neutron.z;

	if (neutron.DTC >= neutron.DTS) {
		logical_x += neutron.u * 1e-6;
		logical_y += neutron.v * 1e-6;
		logical_z += neutron.w * 1e-6;
	}

	int calc_i = static_cast<int>((logical_x + geometry.pitch_i * geometry.size_i / 2.0) / geometry.pitch_i);
	int calc_j = static_cast<int>((logical_y + geometry.pitch_j * geometry.size_j / 2.0) / geometry.pitch_j);
	int calc_k = static_cast<int>((logical_z + geometry.pitch_k * geometry.size_k / 2.0) / geometry.pitch_k);

	if (calc_i >= geometry.size_i) calc_i = geometry.size_i - 1;
	else if (calc_i < 0) calc_i = 0;

	if (calc_j >= geometry.size_j) calc_j = geometry.size_j - 1;
	else if (calc_j < 0) calc_j = 0;

	if (calc_k >= geometry.size_k) calc_k = geometry.size_k - 1;
	else if (calc_k < 0) calc_k = 0;

	coord.i = std::max(0, std::min(calc_i, geometry.size_i - 1));
	coord.j = std::max(0, std::min(calc_j, geometry.size_j - 1));
	coord.k = std::max(0, std::min(calc_k, geometry.size_k - 1));   //check

	neutron.center_x = geometry.pitch_i * (0.5 + coord.i - geometry.size_i / 2.0);
	neutron.center_y = geometry.pitch_j * (0.5 + coord.j - geometry.size_j / 2.0);
	neutron.center_z = geometry.pitch_k * (0.5 + coord.k - geometry.size_k / 2.0);

	neutron.local_x = neutron.x - neutron.center_x;
	neutron.local_y = neutron.y - neutron.center_y;
	neutron.local_z = neutron.z - neutron.center_z;

	neutron.current_index = coord.k * (geometry.size_i * geometry.size_j) + coord.j * geometry.size_i + coord.i;
	neutron.current_cel_id = geometry.distribution[neutron.current_index];

	double logical_local_x = logical_x - neutron.center_x;
	double logical_local_y = logical_y - neutron.center_y;
	double logical_local_z = logical_z - neutron.center_z;


	Method method;

	neutron.current_material = parsing.determine_material(logical_local_x, logical_local_y, logical_local_z, neutron.current_cel_id, geometry);

	if (neutron.current_material == -1) {
		return;
	}

	neutron.DTC = -(log(method.random_number_generator(rn)) / material.materials[neutron.current_material].xs_t[neutron.group]);

	// current location (i, j, k) or (x, y, z) -> use index and cel.id
	if (neutron.current_material == 1 || neutron.current_material == 2 || neutron.current_material == 3 || neutron.current_material == 4 || neutron.current_material == 5 || neutron.current_material == 6) {
		double dts_w = 0.0;
		double dts_r = 0.0;
		if (neutron.w > 0) dts_w = (geometry.pitch_k / 2.0 - neutron.local_z) / neutron.w;
		else if (neutron.w < 0) dts_w = (-1 * geometry.pitch_k / 2.0 - neutron.local_z) / neutron.w;
		else dts_w = 1e30;

		double a = -neutron.u * neutron.local_x - neutron.v * neutron.local_y;
		double b = (neutron.u * neutron.local_x + neutron.v * neutron.local_y) * (neutron.u * neutron.local_x + neutron.v * neutron.local_y) - (neutron.local_x * neutron.local_x + neutron.local_y * neutron.local_y - geometry.pitch_r * geometry.pitch_r) * (neutron.u * neutron.u + neutron.v * neutron.v);
		double c = neutron.u * neutron.u + neutron.v * neutron.v;
		if (b > 0 && c != 0) dts_r = (a + sqrt(b)) / c;
		else dts_r = 1e30;

		if (dts_w < 1e-9) dts_w = 1e30; //check
		if (dts_r < 1e-9) dts_r = 1e30; //check

		neutron.DTS = std::min({ dts_w, dts_r });
	}
	if (neutron.current_material == 0) {
		double dts_u = 0.0;
		double dts_v = 0.0;
		double dts_w = 0.0;
		double dts_r = 0.0;


		if (neutron.u > 0) dts_u = (geometry.pitch_i / 2.0 - neutron.local_x) / neutron.u;
		else if (neutron.u < 0) dts_u = (-1 * geometry.pitch_i / 2.0 - neutron.local_x) / neutron.u;
		else dts_u = 1e30;

		if (neutron.v > 0) dts_v = (geometry.pitch_j / 2.0 - neutron.local_y) / neutron.v;
		else if (neutron.v < 0) dts_v = (-1 * geometry.pitch_j / 2.0 - neutron.local_y) / neutron.v;
		else dts_v = 1e30;

		if (neutron.w > 0) dts_w = (geometry.pitch_k / 2.0 - neutron.local_z) / neutron.w;
		else if (neutron.w < 0) dts_w = (-1 * geometry.pitch_k / 2.0 - neutron.local_z) / neutron.w;
		else dts_w = 1e30;

		double a = -neutron.u * neutron.local_x - neutron.v * neutron.local_y;
		double b = (neutron.u * neutron.local_x + neutron.v * neutron.local_y) * (neutron.u * neutron.local_x + neutron.v * neutron.local_y) - (neutron.local_x * neutron.local_x + neutron.local_y * neutron.local_y - geometry.pitch_r * geometry.pitch_r) * (neutron.u * neutron.u + neutron.v * neutron.v);
		double c = neutron.u * neutron.u + neutron.v * neutron.v;
		if (b > 0 && c != 0) dts_r = (a - sqrt(b)) / c;
		else dts_r = 1e30;

		if (dts_u < 1e-9) dts_u = 1e30; //check
		if (dts_v < 1e-9) dts_v = 1e30;
		if (dts_w < 1e-9) dts_w = 1e30;
		if (dts_r < 1e-9) dts_r = 1e30;

		neutron.DTS = std::min({ dts_u, dts_v, dts_w, dts_r });
	}

}

void Manager::cycle() {
	Method method;

	for (int i = 0; i < current_NPS; ++i) {

		neutron = current_bank[i];

		bool loop_active = true;

		while (loop_active) {

			distance.distance(rn, neutron, material, geometry, parsing, coord);

			if (neutron.current_material == -1) {
				loop_active = false;
				continue;
			}

			if (neutron.DTC < neutron.DTS) { // DTC가 DTS보다 짧으므로 Reaction

				tally_sum += material.materials[neutron.current_material].nu[neutron.group] *
					material.materials[neutron.current_material].xs_f[neutron.group] /
					material.materials[neutron.current_material].xs_t[neutron.group];

				double r = method.random_number_generator(rn) * material.materials[neutron.current_material].xs_t[neutron.group];
				double sum_scattering = std::accumulate(material.materials[neutron.current_material].xs_s.begin() + neutron.group * 7, material.materials[neutron.current_material].xs_s.begin() + (neutron.group * 7 + 7), 0.0);

				if (r < sum_scattering) {	// 반응 중 Scattering 선택
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

				else if (r < sum_scattering + material.materials[neutron.current_material].xs_c[neutron.group]) {	// 반응 중 Absorption
					// loop를 false로 만들어 흡수되어 반응하지 않는 걸 묘사
					loop_active = false;
				}

				else {		// fission
					neutron.x = neutron.x + neutron.u * neutron.DTC;
					neutron.y = neutron.y + neutron.v * neutron.DTC;
					neutron.z = neutron.z + neutron.w * neutron.DTC;

					double safe_k = (k > 0.0) ? k : 1.0;

					int nu_generated = static_cast<int>(material.materials[neutron.current_material].nu[neutron.group] / safe_k + method.random_number_generator(rn));

					generation_bank_count += nu_generated;

					for (int g = 0; g < nu_generated; ++g) {

						neutron.group = 0;
						double accumulated_chi = 0.0;
						double r_chi = method.random_number_generator(rn);

						for (int next_g = 0; next_g < 7; ++next_g) {
							accumulated_chi += material.materials[neutron.current_material].chi[next_g];

							if (r_chi < accumulated_chi) {
								neutron.group = next_g;
								break;
							}
						}
						factory.ini_dir(rn, neutron);
						next_bank.push_back(neutron);
					}
					loop_active = false;
				}
			}
			else {	// DTC>DTS : 반응 안하고 다음 surface로 이동
				neutron.x = neutron.x + neutron.u * (neutron.DTS);
				neutron.y = neutron.y + neutron.v * (neutron.DTS);
				neutron.z = neutron.z + neutron.w * (neutron.DTS);

				if ( neutron.x <= -geometry.pitch_i * geometry.size_i / 2.0 + 1e-9) {
					neutron.u *= -1.0;
				}
				else if (neutron.x >= geometry.pitch_i * geometry.size_i / 2.0 - 1e-9) {
					loop_active = false;
				}

				if (neutron.y <= -geometry.pitch_j * geometry.size_j / 2.0 + 1e-9) {
					neutron.v *= -1.0;
				}
				else if (neutron.y >= geometry.pitch_j * geometry.size_j / 2.0 - 1e-9) {
					loop_active = false;
				}

				if (neutron.z <= -geometry.pitch_k * geometry.size_k / 2.0 + 1e-9) {
					neutron.w *= -1.0;
				}
				else if (neutron.z >= geometry.pitch_k * geometry.size_k / 2.0 - 1e-9) {
					loop_active = false;
				}
			}
			/*
			else {	// DTC>DTS : 반응 안하고 다음 surface로 이동
				neutron.x = neutron.x + neutron.u * (neutron.DTS);
				neutron.y = neutron.y + neutron.v * (neutron.DTS);
				neutron.z = neutron.z + neutron.w * (neutron.DTS);

				double bound_x = geometry.pitch_i * geometry.size_i / 2.0;
				double bound_y = geometry.pitch_j * geometry.size_j / 2.0;
				double bound_z = geometry.pitch_k * geometry.size_k / 2.0;

				// 경계를 넘어가면 중성자가 시스템 밖으로 누설됨 (추적 종료)
				if (std::abs(neutron.x) >= bound_x - 1e-9 ||
					std::abs(neutron.y) >= bound_y - 1e-9 ||
					std::abs(neutron.z) >= bound_z - 1e-9) {
					loop_active = false;
				}
			}
			*/
		}

	}
}


void Manager::iteration(int numNeutron) {

	this -> NPS = numNeutron;

	double active_k_sum = 0.0;
	double active_k_sq_sum = 0.0;
	int active_count = 0;

	this->current_NPS = NPS;
	Method method;

	std::vector<int> accumulated_cell_count(geometry.size_i * geometry.size_j, 0);

	for (int i = 0; i < current_NPS; i++) {
		factory.ini_pos(rn, this->neutron, geometry);
		factory.ini_dir(rn, this->neutron);
		this->neutron.group = static_cast<int>(method.random_number_generator(rn) * 7);
		this->neutron.weight = 1.0;
		

		this->current_bank.push_back(this->neutron);

		if (i < 100) {
			Neutron localN = current_bank[i];
			std::cout << "Neutron IDX " << i << ", pos: [" << localN.x << ", " << localN.y << ", " << localN.z << "],\t";
			std::cout << "dir: [" << localN.u << ", " << localN.v << ", " << localN.w << "], Energy Group: " << localN.group << "\n";
		}

	}

	for (int i = 0; i < total_cycles; ++i) {
		this->generation_bank_count = 0.0;
		this->tally_sum = 0.0;

		this->cycle();


		this->k = tally_sum / static_cast<double>(current_NPS);

		//this->k = static_cast<double>(this->next_bank.size()) / static_cast<double>(current_NPS);


		current_bank = next_bank;
		next_bank.clear();
		next_bank.reserve(this->NPS);
		current_NPS = current_bank.size();



		if (i >= inactive_cycles) {
			active_k_sum += k;
			active_k_sq_sum += (k * k);
			active_count++;

			for (const auto& n : current_bank) {
				int cell_idx_d = n.current_index % (geometry.size_i * geometry.size_j);

				if (cell_idx_d >= 0 && cell_idx_d < accumulated_cell_count.size()) {
					accumulated_cell_count[cell_idx_d]++;
				}
			}

			std::cout << "Active Cycle " << i + 1 << " k_eff: " << k << "\n";
			std::cout << "                " << "current_NPS: " << current_NPS << "\n";
		}
		else {
			std::cout << "Inactive Cycle " << i + 1 << " k_eff: " << k << "\n";
			std::cout << "                " << "current_NPS: " << current_NPS << "\n";
		}
		if (i == inactive_cycles - 1) std::cout << "\n" << "--------------------------------" << "\n";
	}

	std::ofstream outFile("neutron_dist.txt");
	for (int row = 0; row < geometry.size_i; ++row) {
		for (int col = 0; col < geometry.size_j; ++col) {
			int idx = row * geometry.size_j + col;
			outFile << accumulated_cell_count[idx] << " ";
		}
		outFile << "\n";
	}
	outFile.close();
	std::cout << "\nAccumulated Neutron distribution saved to neutron_dist.txt\n";


	if (active_count > 0) {
		double avg_k = active_k_sum / active_count;
		double variance = (active_k_sq_sum / active_count) - (avg_k * avg_k);
		double std_dev = std::sqrt(variance / (active_count - 1));

		std::cout << avg_k << " standard deviation " << std_dev << "\n";
	}
}
