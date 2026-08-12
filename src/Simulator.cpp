#include "../include/Parsing.h"
#include "../include/lcg2.h"
#include "../include/Simulator.h"
#include <algorithm>

void Factory::ini_pos(lcg& rn, Neutron& neutron, Geometry& geometry) {
	Method method;
	neutron.x = method.random_number_generator(rn) * geometry.pitch_i * geometry.size_i * (5 / 3) - geometry.pitch_i * geometry.size_i / 2.0;
	neutron.y = method.random_number_generator(rn) * geometry.pitch_j * geometry.size_j * (5 / 3) - geometry.pitch_j * geometry.size_j / 2.0;
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

	//int lookup_k = (geometry.size_k - 1) - coord.k;
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

			if (neutron.current_material == -1) {     // neutron이 geometry 밖으로 나갈시 종료
				loop_active = false;
				continue;
			}

			tally.L = std::min(neutron.DTC, neutron.DTS);

			// 셀 내부 연료봉 영역 tally

			if (neutron.current_material == 1 || neutron.current_material == 2 || neutron.current_material == 3 || neutron.current_material == 4 || neutron.current_material == 5 || neutron.current_material == 6) {
				double lx = neutron.local_x;
				double ly = neutron.local_y;

				double u = neutron.u;
				double v = neutron.v;

				int rim_id = 0;
				std::vector<double> r = {0.0};
				for (int j = 0; j < 6; ++j) {
					r.push_back(geometry.pitch_r * (j+1) / 6.0 );
				}

				// x = lx + u * t;
				// y = ly + v * t;

				// 1단계

				double A = u * u + v * v;
				double B = 2 * (lx * u + ly * v);
				
				std::vector<double> t = {0.0, tally.L};

				for (int k = 1; k < r.size(); ++k) {
					double C = lx * lx + ly * ly - r[k] * r[k];
					double discriminant = B * B - 4 * A * C;

					if (discriminant > 0) {
						double t1 = (-B + sqrt(discriminant)) / (2 * A);
						double t2 = (-B - sqrt(discriminant)) / (2 * A);

						if (t1 > 0 && t1 < tally.L) t.push_back(t1);
						if (t2 > 0 && t2 < tally.L) t.push_back(t2);
					}
				}

				// 2단계

				std::sort(t.begin(), t.end());  // 오름차순 정렬

				auto it = std::unique(t.begin(), t.end(), [](double a, double b) {  // 중복 제거
					return std::abs(a - b) < 1e-9;
					});
				t.erase(it, t.end());

				std::vector<double> dif_L;
				for (int l = 0; l < t.size() - 1; ++l) {
					dif_L.push_back(t[l + 1] - t[l]);
				}

				std::vector<double> t_mid;
				for (int m = 0; m < t.size() - 1; ++m) {
					t_mid.push_back((t[m] + t[m + 1]) / 2.0);
				}

				// 3단계
				
				for (int n = 0; n < t_mid.size(); ++n) {
					double x_mid = lx + u * t_mid[n];
					double y_mid = ly + v * t_mid[n];
					double r_mid_sqr = x_mid * x_mid + y_mid * y_mid;

					for (int p = 0; p < r.size() - 1; ++p) {
						if (r_mid_sqr >= r[p] * r[p] && r_mid_sqr < r[p + 1] * r[p + 1]) {
							rim_id = p;

							tally.flux_tally[neutron.current_index * material.materials.size() * 6 + neutron.current_material * 6 + rim_id] += neutron.weight * dif_L[n];

							tally.fission_tally[neutron.current_index * material.materials.size() * 6 + neutron.current_material * 6 + rim_id] += neutron.weight
								* material.materials[neutron.current_material].xs_f[neutron.group] * dif_L[n];


							tally.fission_neutron_tally[neutron.current_index * material.materials.size() * 6 + neutron.current_material * 6 + rim_id] += neutron.weight
								* material.materials[neutron.current_material].xs_f[neutron.group]
								* dif_L[n] * material.materials[neutron.current_material].nu[neutron.group];
						}
					
					}

				}
			}

			else {
				tally.flux_tally[neutron.current_index * material.materials.size() * 6 + neutron.current_material * 6 + 0] += neutron.weight * tally.L;
			}

			neutron.x = neutron.x + neutron.u * tally.L;
			neutron.y = neutron.y + neutron.v * tally.L;
			neutron.z = neutron.z + neutron.w * tally.L;

			if (neutron.DTC < neutron.DTS) { // DTC가 DTS보다 짧으므로 Reaction

				int mat_idx = neutron.current_material;   // 현재 중성자가 위치한 물질의 인덱스 ( 효율을 위해, 보기 좋게)
				int g = neutron.group;
				double xs_t = material.materials[mat_idx].xs_t[g];
				double xs_f = material.materials[mat_idx].xs_f[g];
				double nu = material.materials[mat_idx].nu[g];

				// fission


				double safe_k = (tally.track_length_k > 0.0) ? tally.track_length_k : 1.0;
				double expected_neutron = neutron.weight * (nu * xs_f / xs_t) / safe_k;
				int nu_generated = static_cast<int>(expected_neutron + method.random_number_generator(rn));


				tally.collision_tally += neutron.weight * xs_f / xs_t * nu;

				for (int g = 0; g < nu_generated; ++g) {
					Neutron fission_neutron = neutron;
					fission_neutron.weight = 1.0;
					fission_neutron.group = 0;

					double accumulated_chi = 0.0;
					double r_chi = method.random_number_generator(rn);

					for (int next_g = 0; next_g < 7; ++next_g) {
						accumulated_chi += material.materials[neutron.current_material].chi[next_g];

						if (r_chi < accumulated_chi) {
							fission_neutron.group = next_g;
							break;
						}
					}
					factory.ini_dir(rn, fission_neutron);
					next_bank.push_back(fission_neutron);
				}

				// Implicit Capture

				double sum_scattering = std::accumulate(material.materials[neutron.current_material].xs_s.begin() + g * 7, material.materials[neutron.current_material].xs_s.begin() + (g * 7 + 7), 0.0);
				neutron.weight *= (sum_scattering / xs_t);

				// Russian Roulette

				if (neutron.weight < cutoff_weight) {
					double survival_probability = neutron.weight / target_weight;
					if (method.random_number_generator(rn) < survival_probability) {
						neutron.weight = target_weight;
					}
					else {
						loop_active = false;
						continue;
					}
				}

				// Scattering

				if (loop_active) {

					double r = method.random_number_generator(rn) * sum_scattering;
					double accumulated_xs = 0.0;
					for (int next_g = 0; next_g < 7; ++next_g) {
						accumulated_xs += material.materials[mat_idx].xs_s[g * 7 + next_g];

						if (r < accumulated_xs) {
							neutron.group = next_g;
							break;
						}
					}
					

					factory.ini_dir(rn, neutron);
				}

				
			}
			else {	// DTC>DTS : 반응 안하고 다음 surface로 이동

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
		}

	}
}


void Manager::iteration(int numNeutron) {
	Method method;

	this -> NPS = numNeutron;	
	this->current_NPS = NPS;





	// active tally 초기화
	this->tally.active_flux_tally.assign(geometry.total_size * material.total_materials * 6, 0.0);  //
	this->tally.active_fission_tally.assign(geometry.total_size * material.total_materials * 6, 0.0);

	




	// 초기 중성자 생성
	for (int i = 0; i < current_NPS; i++) {
		factory.ini_pos(rn, this->neutron, geometry);
		factory.ini_dir(rn, this->neutron);
		this->neutron.group = static_cast<int>(method.random_number_generator(rn) * 7);
		this->neutron.weight = 1.0;
		
		this->current_bank.push_back(this->neutron);

		if (i < 10) {
			Neutron localN = current_bank[i];
			std::cout << "Neutron IDX " << i << ", pos: [" << localN.x << ", " << localN.y << ", " << localN.z << "],\t";
			std::cout << "dir: [" << localN.u << ", " << localN.v << ", " << localN.w << "], Energy Group: " << localN.group << "\n";
		}

	}







	// 반복 계산 수행
	for (int i = 0; i < total_cycles; ++i) {

		// 사이클마다 초기화
		tally_manager.reset_cycle_tally(this->tally, this->geometry.total_size, this->material.total_materials);

		// 사이클 수행
		this->cycle();

		double tally_sum = std::accumulate(this->tally.fission_neutron_tally.begin(), this->tally.fission_neutron_tally.end(), 0.0);

		this->tally.track_length_k = tally_sum / static_cast<double>(current_NPS);  		//this->k = static_cast<double>(this->tally.next_bank.size()) / static_cast<double>(current_NPS);
		this->tally.collision_k = this->tally.collision_tally / static_cast<double>(current_NPS);



		this->current_bank = this->next_bank;
		this->next_bank.clear();
		this->next_bank.reserve(this->NPS);
		current_NPS = this->current_bank.size();


		// cycle 출력
		if (i >= inactive_cycles) {
			tally_manager.accumulate_active_tally(this->tally, current_NPS);

			std::cout << "Active Cycle " << i + 1 << " k_eff: " << tally.track_length_k << "\n";
			std::cout << "                " << "current_NPS: " << current_NPS << "\n";
		}
		else {
			std::cout << "Inactive Cycle " << i + 1 << " k_eff: " << tally.track_length_k << "\n";
			std::cout << "                " << "current_NPS: " << current_NPS << "\n";
		}
		if (i == inactive_cycles - 1) std::cout << "\n" << "--------------------------------" << "\n";
	}




	// active cycle 통계 계산
	tally_manager.statistics(this->tally, this->geometry, this->material);


	tally_manager.axial_distribution(this->tally, this->geometry, this->material);
	tally_manager.radial_distribution(this->tally, this->geometry, this->material);

	tally_manager.export_distributions(this->tally);
}