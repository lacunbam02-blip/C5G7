#pragma once
#include <iostream>
#include <vector>
#include "Parsing.h"


class Tally {
public:
	
	// 매 사이클 사용 데이터
	double L = 0.0;

	std::vector<double> flux_tally;
	std::vector<double> fission_tally;
	std::vector<double> fission_neutron_tally;

	double track_length_k = 1.0;    // k 값

	double collision_tally = 0.0;
	double collision_k = 1.0;




	// ative cycle 누적 통계
	double active_NPS_sum = 0.0;

	double active_track_length_k_sum = 0.0;
	double active_collision_k_sum = 0.0;

	double track_length_variance = 0.0;
	double collision_variance = 0.0;

	double active_track_length_k_sq_sum = 0.0;
	double active_collision_k_sq_sum = 0.0;

	int active_count = 0;
	std::vector<double> active_fission_tally;
	std::vector<double> active_flux_tally;




	//최종 결과
	double avg_track_length_k = 0.0;
	double avg_collision_k = 0.0;

	double std_dev_track_length_k = 0.0;
	double std_dev_collision_k = 0.0;

	double avg_track_length_NPS = 0.0;
	double avg_collision_NPS = 0.0;

	std::vector<double> mean_flux;
	std::vector<double> mean_fission;




	// axial 분포(UO2)
	std::vector<double> axial_flux_distribution;
	std::vector<double> axial_fission_distribution;

	// radial 분포
	std::vector<double> radial_flux_distribution;
	std::vector<double> radial_fission_distribution;
	std::vector<double> radial_fission_tally_pin; // 각 핀(셀) 단위로 핵분열 합계 저장

};

class Tally_Manager {
public:
	inline int get_idx(int cell_id, int material_id, int rim_id, int group) {
		return cell_id * (7 * 6 * 7) + material_id * (6 * 7) + rim_id * (7) + group;
	}
	void reset_cycle_tally(Tally& tally, int total_cells, int total_materials, int total_rims, int total_groups);
	void accumulate_active_tally(Tally& tally, int current_NPS);
	void statistics(Tally& tally, Geometry& geometry, Material& material);
	void axial_distribution(Tally& tally, Geometry& geometry, Material& material);
	void radial_distribution(Tally& tally, Geometry& geometry, Material& material);
	void ring_distribution(Tally& tally, Geometry& geometry, Material& material);
	void export_distributions(const Tally& tally);
};