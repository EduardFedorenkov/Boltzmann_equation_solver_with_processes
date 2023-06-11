//============================================================================
// Name        : Boltzmann_simple_transport.cpp
// Author      : Fedorenkov Eduard
// Version     : 1.0
// Copyright   : BINP Lab 9-0
// Description : Kinetic code for Gas-pasma interaction.
//============================================================================

#include "profile.h"
#include "test_runner.h"
#include "data_saver.h"
#include "distribution_func.h"
#include "elastic_collisions.h"
#include "collisions.h"
#include "processes.h"
#include "medium.h"
#include "gas.h"
#include <cmath>

void TestVelocityGrid(){
	VelocityGrid v(11,5);
	ASSERT_EQUAL(v.Get1DGrid(), vector<double>({-5, -4, -3, -2, -1, 0, 1, 2, 3, 4, 5}));
	ASSERT_EQUAL(v.GetGridStep(), 1);
	ASSERT_EQUAL(v.GetSize(), 11u);
	ASSERT_EQUAL(v.GetMax(), 5);
	VelocityGrid v2(vector<double>({-3.33, -2.22, -1.11, 0, 1.11, 2.22, 3.33}));
	ASSERT_EQUAL(v2.Get1DGrid(), vector<double>({-3.33, -2.22, -1.11, 0, 1.11, 2.22, 3.33}));
	DOUBLES_ASSERT_EQUAL_ACCURATE_TO(datum::eps, v2.GetGridStep(), 1.11);
	ASSERT_EQUAL(v2.GetSize(), 7u);
	ASSERT_EQUAL(v2.GetMax(), 3.33);
	VelocityGrid v3(size_t(7), 1.0, 2.0 * datum::c_0 * datum::c_0 * 1e4);
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps, vec(v3.Get1DGrid()), vec({-3,-2,-1,0,1,2,3}));
	VelocityGrid v4(v3.GetGridStep(), 1.0, 2.0 * datum::c_0 * datum::c_0 * 1e4);
	CONTAINERS_ASSERT_EQUAL_ACCURATE_TO(datum::eps, vec(v3.Get1DGrid()), vec(v4.Get1DGrid()));
}

void TestSpaceGrid(){
	SpaceGrid s;
	ASSERT_EQUAL(s.Get1DGrid(), vector<double>({0}));
	SpaceGrid x(10, 10, make_pair(BC_Type::PerfectReflection, BC_Type::ConstantTemperatureWall), make_pair(0.0, 300 * datum::k_evk));
	ASSERT_EQUAL(x.Get1DGrid(), vector<double>({0.5, 1.5, 2.5, 3.5, 4.5, 5.5, 6.5, 7.5, 8.5, 9.5}));
	ASSERT_EQUAL(x.GetDistance(), 9.5);
	ASSERT_EQUAL(x.GetGridStep(), 1);
	ASSERT_EQUAL(x.GetSize(), 10u);
	ASSERT(x.GetWalls().walls_BC.first == BC_Type::PerfectReflection);
	ASSERT(x.GetWalls().walls_BC.second == BC_Type::ConstantTemperatureWall);
	ASSERT_EQUAL(x.GetWalls().walls_T.first, 0.0);
	ASSERT_EQUAL(x.GetWalls().walls_T.second, 300 * datum::k_evk);
}

void SaveProcessData(const Plasma& p, const DistributionFunction& df, const size_t N_points, const double E){
	double Tg = *df.ComputeTemperatureStatic().begin();
	size_t vg_size = 11;
	VelocityGrid vg(vg_size, Tg, datum::m_p * datum::c_0 * datum::c_0 * 1e4);
	VelocityGrid vp(vg_size, Tg, datum::m_p * datum::c_0 * datum::c_0 * 1e4);
	Charge_exchange cx("ChargeExchange.txt", p, df);
	//He_ionization ioniz("ElectronIonization.txt", p);
	//HHplus_elastic hp_elastic("HpElastic.txt", vg, vp, p);
	//HH_elastic hh_elastic("HHElastic.txt", vg);
	cx.SaveCXRateCoeff(E, N_points);
	//ioniz.SaveIonizRateCoeff(N_points);
	//hp_elastic.SaveDiffCross(E, N_points);
	//hh_elastic.SaveHHDiffCross(E, N_points);
}

void CXCheck(const Plasma& p, const DistributionFunction& df){
	Charge_exchange cx("ChargeExchange.txt", p, df);
	auto rhs = cx.ComputePGRightHandSide(p, df);
	//size_t mid_point = (df.GetVelGrid().GetSize() - 1) / 2;
	cout << rhs[0] << endl;
	cout << p.GetDensity(0) << endl;
	cout << p.GetDensity(0);
	//auto plasma_distr = p.MakeMaxwellDistr(0, df.GetVelGrid().Get1DGrid());
	//cout << plasma_distr.slice(mid_point) << endl;
}

void HeElasticCheck(const Plasma& p){
	He_elastic Hee("HeElastic.txt", p.GetIonMass(), p);
	cout << datum::m_p * datum::c_0 * datum::c_0 / datum::eV / Hee.GetDiffusCoeff(0);
}

double TheoreticalSpitzerTestForce(const double n, const double T, const vec3& Ve, const bool IsVe){
	const double el_ch = datum::ec * 10 * datum::c_0;
	const double vel_sqr_for_lambda = IsVe ? norm(Ve) * norm(Ve) : 2 * T * 1.6e-12 / (datum::m_e * 1e3);
	const double Lambda_ee = std::log(0.5 * datum::m_e * 1e3 * vel_sqr_for_lambda *
			sqrt(T * 1.6e-12 / (8 * datum::pi * n)) / pow(el_ch, 3));
	const double tau_e = (3 * sqrt(datum::m_e * 1e3 * std::pow(T * 1.6e-12,3) ) ) /
			(4 * datum::sqrt2pi * n * std::pow(el_ch,4) * Lambda_ee);
	// Validate crocc section value:
	cout << "Lambda = " << Lambda_ee << endl;
	cout << "tau_e = " << tau_e << endl;
	//
	return datum::m_e * 1e3 * n * norm(Ve) / tau_e;
}

int main() {
	// Testing procedure
	TestRunner tr;
	RUN_TEST(tr, TestVelocityGrid);
	RUN_TEST(tr, TestSpaceGrid);

	// Gas params
	double Tg = 1;
	double ng = 1e14;
	double mg = datum::m_e * datum::c_0 * datum::c_0 / datum::eV;

	// Plasma params
	double Tp = 1;
	double np = 1e14;
	double mi = datum::m_e * datum::c_0 * datum::c_0 / datum::eV;
	vec3 Vp{1e6, 0, 0};
	Plasma p(mi, Tp, np, Vp);
	cout << "np = " << p.GetDensity(0) << endl;
	cout << "Vt = " << p.GetTermalVel(0) << endl;
	cout << "Tp = " << p.GetTemperature(0) << endl;

	// Velocity grid params
	size_t v_size = 11;
	VelocityGrid v(v_size, Tg, mg);

	// Time Evolution
	size_t N_steps = 1;
	cout << "Vp = " << Vp << ' ' <<  "V_grid = "  << *v.Get1DGrid().rbegin() << endl;
	cout << "vel_grid_step = " << v.GetGridStep() << endl;
	cout << "Vtp = " << p.GetTermalVel(0) << endl;

	std::cout << "Force_theory = " << TheoreticalSpitzerTestForce(np, Tp, Vp, false) << std::endl;

	DistributionFunction f_H(DistributionType::Maxwell, mg, v, ng, Tg);
	Gas H(f_H, vector<shared_ptr<PlasmaGasProcess>>({make_shared<GFastIons_elastic>(200, mg, p, v)}), vector<shared_ptr<GasGasProcess>>{});
	cout << "OK" << endl;
	H.SaveDistr(0, 0);
	for(size_t t = 0; t < N_steps; ++t){
		H.TimeEvolution_ConstTimeStep(p, 1e-7);
		H.SaveDistr(0, t+1);
	}
	auto F = H.SpitzerTestForce(p);
	std::cout << F.at(0) << std::endl;

	//SaveProcessData(p, f_H, 100, 50);
	//CXCheck(p, f_H);
	//HeElasticCheck(p);

	return 0;
}
