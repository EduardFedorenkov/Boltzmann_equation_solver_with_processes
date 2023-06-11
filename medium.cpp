#include "medium.h"
#include "velocity_grid.h"
#include <cmath>

Plasma::Plasma(const double mass, const vector<double>& T, const vector<double>& n) : Tp(T), np(n),
Vp(vector<vec3>(Tp.size(), vec3({0, 0, 0}))), ion_mass(mass) {}

Plasma::Plasma(const double mass, const double T, const double n) : Tp(vector<double>{T}), np(vector<double>{n}),
Vp(vector<vec3>(1, vec3({0, 0, 0}))), ion_mass(mass){}

Plasma::Plasma(const double mass, const double T, const double n, const vec3 Vp) : Tp(vector<double>{T}), np(vector<double>{n}),
		Vp(vector<vec3>{Vp}), ion_mass(mass) {}

vec Plasma::MakeVel1DGrid(const size_t sg_idx, const size_t v_size) const{
	double Vmax = 3.5 * GetTermalVel(sg_idx) / datum::sqrt2;
	return vec(Make_1D_v_grid(v_size, Vmax));
}

cube Plasma::MakeSpitzerCondactDistr(const size_t sg_idx, const size_t v_size) const{
	vec vel_1D = MakeVel1DGrid(sg_idx, v_size);
	cube distr_maxwell = MakeMaxwellDistr(sg_idx, v_size);
	cube distr = distr_maxwell;
	double termal_vel = sqrt(2 * Tp[sg_idx] / ion_mass * datum::c_0 * datum::c_0 * 1e4);
	double factor = sqrt(datum::pi) * 0.25 / (termal_vel * termal_vel * termal_vel * termal_vel * termal_vel);
	for(size_t k = 0; k < v_size; k++){
		for(size_t l = 0; l < v_size; l++){
			for(size_t m = 0; m < v_size; m++){
				vec3 vp({vel_1D(m), vel_1D(l), vel_1D(k)});
				double norm_vp = norm(vp);
				double full_factor = dot(vp, Vp[sg_idx]) * norm_vp * norm_vp * norm_vp * factor;
				distr(m,l,k) += full_factor * distr_maxwell(m,l,k);
			}
		}
	}
	return distr;
}

cube Plasma::MakeMaxwellDistr(const size_t sg_idx, const size_t v_size) const{
	vec vel_1D = MakeVel1DGrid(sg_idx, v_size);
	cube distr(v_size,v_size,v_size,fill::zeros);
	double sqr_termal_vel = 2 * Tp[sg_idx] / ion_mass * datum::c_0 * datum::c_0 * 1e4;
	double factor = np[sg_idx] / pow(sqrt(datum::pi * sqr_termal_vel),3);
	for(size_t k = 0; k < v_size; k++){
		for(size_t l = 0; l < v_size; l++){
			for(size_t m = 0; m < v_size; m++){
				distr(m,l,k) = factor * exp(- (vel_1D(m)*vel_1D(m) + vel_1D(k)*vel_1D(k) + vel_1D(l)*vel_1D(l)) / sqr_termal_vel);
			}
		}
	}
	return distr;
}

double Plasma::GetIonElectronCollisionsFrequency(const size_t sg_idx, const size_t v_size, const vec3 vi) const{
	vec vel_1D = MakeVel1DGrid(sg_idx, v_size);
	cube plasma_distr = MakeMaxwellDistr(sg_idx, v_size);

	double sqr_termal_vel = 2 * Tp[sg_idx] / ion_mass * datum::c_0 * datum::c_0 * 1e4;
	double factor = np[sg_idx] / pow(sqrt(datum::pi * sqr_termal_vel),3);
	double frequency = 0;

	for(size_t k = 0; k < v_size; k++){
		for(size_t l = 0; l < v_size; l++){
			for(size_t m = 0; m < v_size; m++){
				vec3 vp = Vp[sg_idx] + vec3({vel_1D(m), vel_1D(l), vel_1D(k)});
				double u = norm(vp - vi);
				double plasmaDistr = factor * exp(- (vel_1D(m)*vel_1D(m) + vel_1D(k)*vel_1D(k) + vel_1D(l)*vel_1D(l)) / sqr_termal_vel);

				frequency += u * plasmaDistr * GetIonElectronCrossSection(sg_idx, u);
			}
		}
	}
	return frequency;
}

double Plasma::GetTemperature(const size_t idx) const{
	return Tp[idx];
}

double Plasma::GetDensity(const size_t idx) const{
	return np[idx];
}

double Plasma::GetIonElectronCrossSection(const size_t sg_idx, const double u) const{
	double mu = ion_mass * datum::m_e / (ion_mass + datum::m_e);
	double charge = datum::ec * 3 * 1e9; 					//[SGSE units]
	double crossSection = datum::pi * (Tp[sg_idx]/(8 * datum::pi * charge * charge * np[sg_idx])
			- charge * charge * charge * charge / (mu * mu * u * u * u * u));
	if(crossSection < 0)
		cerr << "Electron-Ion cross section is negative!!!";
	return crossSection;
}

vec3 Plasma::GetVel(const size_t idx) const{
	return Vp[idx];
}

double Plasma::GetTermalVel(const size_t idx) const{
	return sqrt(2 * Tp[idx] / ion_mass) * datum::c_0 * 100;
}

size_t Plasma::GetSpaceSize() const{
	return Tp.size();
}

double Plasma::GetIonMass() const{
	return ion_mass;
}
