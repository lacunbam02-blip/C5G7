#include "lcg2.h"

lcg::lcg(unsigned long long seed) {
	lamda = 6364136223846793005ULL;
	mu = 1ULL;
	ran = seed;
}

double Method::random_number_generator(lcg& rn) {
	rn.ran = (rn.lamda * rn.ran + rn.mu) & 0xFFFFFFFFFFFFULL;  // '0xFF' makes rest of value
	return static_cast<double>(rn.ran) / 281474976710656.0;   //281474976710656.0 : maximum value of 48bit  --> normalization
}