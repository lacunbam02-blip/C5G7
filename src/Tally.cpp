#include "../include/Tally.h"
#include <numeric>
#include <cmath>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <iomanip>


// 탈리 초기화
void Tally_Manager::reset_cycle_tally(Tally& tally, int total_cells, int total_materials) {

    tally.fission_neutron_tally.assign(total_cells * total_materials * 6, 0.0);
    tally.fission_tally.assign(total_cells * total_materials * 6, 0.0);
    tally.flux_tally.assign(total_cells * total_materials * 6, 0.0);

	tally.collision_tally = 0.0;
}


// k_eff 계산
void Tally_Manager::accumulate_active_tally(Tally& tally, int current_NPS) {
    tally.active_NPS_sum += current_NPS;

    tally.active_track_length_k_sum += tally.track_length_k;
    tally.active_track_length_k_sq_sum += (tally.track_length_k * tally.track_length_k);

    tally.active_collision_k_sum += tally.collision_k;
    tally.active_collision_k_sq_sum += (tally.collision_k * tally.collision_k);

    tally.active_count++;

    for (int idx = 0; idx < tally.flux_tally.size(); ++idx) {
        tally.active_flux_tally[idx] += tally.flux_tally[idx];
        tally.active_fission_tally[idx] += tally.fission_tally[idx];
    }
}

void Tally_Manager::statistics(Tally& tally, Geometry& geometry, Material& material) {
    if (tally.active_count <= 0) return;

	// 평균 및 표준편차 계산
    tally.avg_track_length_k = tally.active_track_length_k_sum / tally.active_count;
    tally.track_length_variance = (tally.active_track_length_k_sq_sum / tally.active_count) - (tally.avg_track_length_k * tally.avg_track_length_k);
    tally.std_dev_track_length_k = std::sqrt(tally.track_length_variance / (tally.active_count - 1));

    tally.avg_collision_k = tally.active_collision_k_sum / tally.active_count;
    tally.collision_variance = (tally.active_collision_k_sq_sum / tally.active_count) - (tally.avg_collision_k * tally.avg_collision_k);
    tally.std_dev_collision_k = std::sqrt(tally.collision_variance / (tally.active_count - 1));
    
    tally.avg_track_length_NPS = tally.active_NPS_sum / tally.active_count;

    double normalization_factor = static_cast<double>(tally.active_count) * tally.avg_track_length_NPS;
    tally.mean_flux.assign(geometry.total_size * material.total_materials * 6, 0.0);
    tally.mean_fission.assign(geometry.total_size * material.total_materials * 6, 0.0);

    for (int idx = 0; idx < tally.active_flux_tally.size(); ++idx) {
        tally.mean_flux[idx] = tally.active_flux_tally[idx] / normalization_factor;
        tally.mean_fission[idx] = tally.active_fission_tally[idx] / normalization_factor;
    }

    std::cout << "\n========================================\n";
    std::cout << "Average k_eff: " << tally.avg_track_length_k << " (std dev: " << tally.std_dev_track_length_k << ")\n";
    std::cout << "Average NPS: " << tally.avg_track_length_NPS << "\n";
	std::cout << "Average Collision k_eff: " << tally.avg_collision_k << " (std dev: " << tally.std_dev_collision_k << ")\n";
    std::cout << "========================================\n";
}



// Axial 분포 계산
void Tally_Manager::axial_distribution(Tally& tally, Geometry& geometry, Material& material) {
	tally.axial_flux_distribution.assign(geometry.size_k, 0.0);
	tally.axial_fission_distribution.assign(geometry.size_k, 0.0);

	std::vector<double> axial_flux_sum(geometry.size_k, 0.0);
	std::vector<double> axial_fission_sum(geometry.size_k, 0.0);

	int plane_size = geometry.size_i * geometry.size_j * material.total_materials * 6;
    int normalization_factor = 0.0;

    for (int k = 0; k < geometry.size_k; ++k) {
        for (int i = 0; i < plane_size; ++i) {
            if ((i/6) % material.total_materials != 0) {
                axial_flux_sum[k] += tally.mean_flux[k * plane_size + i];
                axial_fission_sum[k] += tally.mean_fission[k * plane_size + i];
                normalization_factor++;
            }
        }
        tally.axial_flux_distribution[k] = axial_flux_sum[k] / normalization_factor;
        tally.axial_fission_distribution[k] = axial_fission_sum[k] / normalization_factor;
        normalization_factor = 0;
    }
}


// Radial 분포 계산
void Tally_Manager::radial_distribution(Tally& tally, Geometry& geometry, Material& material) {
    int num_cells = geometry.size_i * geometry.size_j;
    int sub_bins = material.total_materials * 6;
    int plane_size = num_cells * sub_bins;


    tally.radial_flux_distribution.assign(num_cells * 7, 0.0);
    tally.radial_fission_distribution.assign(num_cells * 7, 0.0);

    tally.radial_fission_tally_pin.assign(num_cells * 2, 0.0);    // 각 셀의 핵분열 tally 합계 초기화

    for (int k = 0; k < geometry.size_k; ++k) {
        for (int c = 0; c < num_cells; ++c) {
            for (int mat = 0; mat < material.total_materials; ++mat) {
                if (mat != 0) {
                    for (int r = 0; r < 6; ++r) {
                        int idx = k * plane_size + c * sub_bins + mat * 6 + r;

                        tally.radial_flux_distribution[c * 7 + 1 + r] += tally.mean_flux[idx];
                        tally.radial_fission_distribution[c * 7 + 1 + r] += tally.mean_fission[idx];

                        tally.radial_fission_tally_pin[c * 2 + 1] += tally.mean_fission[idx];  // 각 셀의 핵분열 tally 합계 누적
                    }
                }
                else {
                    tally.radial_flux_distribution[c * 7] += tally.mean_flux[k * plane_size + c * sub_bins];
                    tally.radial_fission_distribution[c * 7] += tally.mean_fission[k * plane_size + c * sub_bins];

                    tally.radial_fission_tally_pin[c * 2] += tally.mean_fission[k * plane_size + c * sub_bins];
                }
            }   
        }
    }

    // Z축 층 수로 나누어 층당 평균값 계산
    for (int i = 0; i < num_cells * 7; ++i) {
        tally.radial_flux_distribution[i] /= geometry.size_k;
        tally.radial_fission_distribution[i] /= geometry.size_k;
    }
	for (int i = 0; i < num_cells * 2; ++i) {
		tally.radial_fission_tally_pin[i] /= geometry.size_k;
	}

    // 2. 정규화 (연료 핀 단위)
    double total_fission = 0.0;
    double total_flux = 0.0;

    int active_fuel_count = 0;

    for (int c = 0; c < num_cells; ++c) {
        double pin_fission_sum = 0.0;
        double pin_flux_sum = 0.0;

        // 정규화를 위해 "이 핀(셀)의 6개 껍질 핵분열 합"을 구함
        for (int r = 0; r < 7; ++r) {
            pin_fission_sum += tally.radial_fission_distribution[c * 7 + r];
            pin_flux_sum += tally.radial_flux_distribution[c * 7 + r];
        }

        // 핵분열이 일어난 셀(연료 핀)인 경우 카운트
        if (pin_fission_sum > 0.0) {
            total_fission += pin_fission_sum;
            total_flux += pin_flux_sum;

            active_fuel_count++;
        }
    }

    // 3. 연료 셀의 평균 계산 후 각 껍질 값들을 정규화
    if (active_fuel_count > 0 && total_fission > 0.0) {
        double avg_fission = total_fission / active_fuel_count; // 연료 핀 1개당 평균 출력
        double avg_flux = total_flux / active_fuel_count; // 연료 핀 1개당 평균 플럭스


        for (int i = 0; i < num_cells * 7; ++i) {
            tally.radial_fission_distribution[i] /= avg_fission;
            tally.radial_flux_distribution[i] /= avg_flux;
        }
		for (int i = 0; i < num_cells * 2; ++i) {
			tally.radial_fission_tally_pin[i] /= avg_fission;
		}
    }
}

void Tally_Manager::ring_distribution(Tally& tally, Geometry& geometry, Material& material) {
    // 연료봉 단위 평균 상대 출력

    int num_cells = geometry.size_i * geometry.size_j;
    std::vector<double> ring_flux_sum(7, 0.0);
    std::vector<double> ring_fission_sum(7, 0.0);
    int active_fuel_pin_count = 0;

    if (tally.radial_fission_distribution.size() >= static_cast<size_t>(num_cells * 7)) {
        for (int c = 0; c < num_cells; ++c) {
            double cell_fission_sum = 0.0;
            for (int r = 0; r < 7; ++r) {
                cell_fission_sum += tally.radial_fission_distribution[c * 7 + r];
            }
            if (cell_fission_sum > 0.0) {
                active_fuel_pin_count++;
                for (int r = 0; r < 7; ++r) {
                    ring_flux_sum[r] += tally.radial_flux_distribution[c * 7 + r];
                    ring_fission_sum[r] += tally.radial_fission_distribution[c * 7 + r];
                }
            }
        }
    }

    std::cout << "--------------------------------------------------------\n";
    std::cout << " Ring ID |   Avg Relative Flux   |  Avg Relative Fission (%) \n";
    std::cout << "---------|-----------------------|----------------------\n";

    if (active_fuel_pin_count > 0) {
        for (int r =1; r < 7; ++r) {
            double avg_flux = ring_flux_sum[r] / active_fuel_pin_count;
            double avg_fiss = ring_fission_sum[r] / active_fuel_pin_count;

            std::cout << "  Ring " << r << " |        " << std::fixed << std::setprecision(5) << avg_flux << "        |        " << avg_fiss << "\n";
        }
    }
}



void Tally_Manager::export_distributions(const Tally& tally) {
    // 1. Axial Flux
    std::ofstream f_ax_flux("axial_flux.txt");
    for (double val : tally.axial_flux_distribution) f_ax_flux << val << "\n";
    f_ax_flux.close();

    // 2. Axial Fission
    std::ofstream f_ax_fiss("axial_fission.txt");
    for (double val : tally.axial_fission_distribution) f_ax_fiss << val << "\n";
    f_ax_fiss.close();

    // 3. Radial Flux
    std::ofstream f_rad_flux("radial_flux.txt");
    for (double val : tally.radial_flux_distribution) f_rad_flux << val << " ";
    f_rad_flux.close();

    // 4. Radial Fission
    std::ofstream f_rad_fiss("radial_fission.txt");
    for (double val : tally.radial_fission_distribution) f_rad_fiss << val << " ";
    f_rad_fiss.close();

    // 5. Radial Fission (per pin)
    std::ofstream f_rad_fiss_pin("radial_fission_pin.txt");
    for (double val : tally.radial_fission_tally_pin) f_rad_fiss_pin << val << " ";
    f_rad_fiss_pin.close();

    std::cout << "\nSuccessfully exported 5 distribution text files!\n";
}