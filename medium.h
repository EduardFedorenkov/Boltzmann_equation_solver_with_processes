#pragma once

#include <armadillo>
#include <vector>

using namespace arma;
using namespace std;

class Plasma{
public:
	Plasma(const double mass, const vector<double>& T, const vector<double>& n);

	Plasma(const double mass, const double T, const double n);

	Plasma(const double mass, const double T, const double n, const vec3 Vp);

	// this function makes 3.5 sigma maxwell distribution vel_1D size.
	vec MakeVel1DGrid(const size_t sg_idx, const size_t v_size) const;

	cube MakeSpitzerCondactDistr(const size_t sg_idx, const size_t v_size) const;

	cube MakeMaxwellDistr(const size_t sg_idx, const size_t v_size) const;

	double GetIonElectronCollisionsFrequency(const size_t sg_idx, const size_t v_size, const vec3 vi) const;

	size_t GetSpaceSize() const;

	double GetTemperature(const size_t idx) const;

	double GetDensity(const size_t idx) const;

	double GetIonElectronCrossSection(const size_t sg_idx, const double u) const;

	vec3 GetVel(const size_t idx) const;

	double GetTermalVel(const size_t idx) const;

	double GetIonMass() const;

private:
	vector<double> Tp;
	vector<double> np;
	vector<vec3> Vp;
	double ion_mass;
};
