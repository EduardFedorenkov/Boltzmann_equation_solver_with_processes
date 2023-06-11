#pragma once

#include <armadillo>
#include <vector>
#include "medium.h"
#include "light_consts.h"

struct ExcToFinal{
	HeLvls final_lvl;
	double probability;
};

class Light{
public:
	Light();

	Light(const size_t z_);

	arma::cube GetExitationAtomsDF(const double ne, const double Te, const arma::cube& df_slice, const HeLvls ex_lvl) const;

private:

	double GetExcitationFrequency(const double ne, const double Te, const HeLvls Term) const;

	size_t z; // spectroscopic symbol (for neutral z = 1; for ions z = q + 1, where q - ion charge)
	std::vector<std::vector<ExcToFinal>> excitation_lvl_to_flvl;
	/*
	 * index of the vector means number in "HeLvls" enum class.
	 * vector<ExcToFinal> is a final lvls and probability to fall in it.
	 * example: excitation_lvl_to_flvl[7] = {{3, 0.7}, {0, 0.3}}
	 */
};
