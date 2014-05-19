//============================================================================
// Name        : NOSE-rates-stochastic.cpp
// Author      : Jan Olsina
// Version     :
// Copyright   : 
// Description : Changing rates stochastic
//============================================================================

#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <random>
#include <chrono>

#include "NOSE-rates-stochastic.h"

using namespace std;

const double
		g_lambda = 100/ENERGY_INTERNAL_TO_CM,
		g_tauC = 100,
		g_temp = 300,
		g_step = 0.5,
		g_J = 50/ENERGY_INTERNAL_TO_CM;

const long g_sites = 3, g_steps = 1024*10, g_interpolation_steps = 1024, g_trajectories = 1000000;
site **mol;

int main(void) {

	init();
	process();
	return EXIT_SUCCESS;

	for(int i = 0; i < g_steps; i++) {
		double t = i*g_step;
		printf("%f %e %e\n", t,
				(mol[0]->m_gf)[i].real(),
				(mol[0]->m_gf)[i].imag());
	}

	return EXIT_SUCCESS;
}

dpc gfunctMatsubara(double t, double lambda, double tauC, double temp) {
	const double LL = 1.0/tauC;
	const double beta = 1.0 / KB / temp;
	const long n = 100;

	dpc out = 0;
	for(int i = n; i > 0; i--) { // from smaller terms for better convergence
		double freq = 2 * M_PI * i / beta;
		out += freq * exp(-freq * t) / (freq*freq - LL* LL);
	}
	out *= 4 * lambda * LL / beta;

	out += lambda * LL * exp(-t * LL) * ( 1/tan(LL * beta / 2) - dpc(0.0,1.0) );

	return out;
}

double basicRateElement(int donorN, int acceptorM, int donorTimeIndex, int transferTimeIndex) {

	const int 	n = donorN,										// shortcuts
				m = acceptorM,
				i = donorTimeIndex + transferTimeIndex,
				j = transferTimeIndex;

	dpc dummy = exp(
			-mol[n]->m_gf[j]-mol[m]->m_gf[j]					// basic g-functions
			-conj(mol[n]->m_gf[i]) + mol[n]->m_gf[i]									// reorganization energy g-funct
			-mol[n]->m_gf[i-j] + conj(mol[n]->m_gf[i-j])								// reorganization energy g-funct
			-(mol[m]->m_en - mol[n]->m_en)*g_step*j*dpc(0,1)	// energy difference
			);

	return 2 * dummy.real() * pow(mol[n]->m_J[m],2.0) * g_step;
}

double thirdOrderRateElement(int elderK, int donorN, int acceptorM, int elderTimeIndex, int donorTimeIndex, int transferTimeIndex) {

	const int 	n = donorN,										// shortcuts
				m = acceptorM,
				k = elderTimeIndex,
				i = donorTimeIndex + transferTimeIndex,
				j = transferTimeIndex;

	if(elderK == donorN)
		return basicRateElement(donorN, acceptorM, elderTimeIndex + donorTimeIndex, transferTimeIndex);

	if(elderK != acceptorM || k + i > g_steps)
		return basicRateElement(donorN, acceptorM, donorTimeIndex, transferTimeIndex);

	dpc dummy = exp(
			-mol[n]->m_gf[j]-mol[m]->m_gf[j]					// basic g-functions
			-conj(mol[n]->m_gf[i]) + mol[n]->m_gf[i]									// reorganization energy g-funct
			-mol[n]->m_gf[i-j] + conj(mol[n]->m_gf[i-j])								// reorganization energy g-funct
			-(mol[m]->m_en - mol[n]->m_en)*g_step*j*dpc(0,1)	// energy difference
			);

	dpc extension = exp(
			mol[m]->m_gf[k+i] - mol[m]->m_gf[i] - mol[m]->m_gf[i-j+k] + mol[m]->m_gf[i-j]
			);

	dpc extension2 = exp(-conj(
			mol[m]->m_gf[k+i] - mol[m]->m_gf[i] - mol[m]->m_gf[i-j+k] + mol[m]->m_gf[i-j]
			));

	dummy = dummy * extension * extension2;

	return 2 * dummy.real() * pow(mol[n]->m_J[m],2.0) * g_step;
}

void init() {
	mol = (site **)malloc(sizeof(site*)*g_sites);
	for(int i = 0; i < g_sites; i++) {
		mol[i] = new site(g_steps, g_sites);

		mol[i]->m_en = 0;

		if(i+1 < g_sites) 	mol[i]->m_J[(i+1)] = g_J;
		if(i-1 >= 0) 		mol[i]->m_J[(i-1)] = g_J;
	}

	#pragma omp parallel for
	for(int n = 0; n < g_sites; n++) {
		for(int m = 0; m < g_sites; m++) {
			if(m == n) continue;

			for(long i = 1; i < g_steps; i++) {
				for(long j = 1; j <= i; j++) {
					mol[n]->m_rates[m][i] += basicRateElement(n, m, i-j, j);
					//printf("%d %d %d %d %f %f\n", n,m,i,j,mol[n]->m_rates[m][i], dummy.real());
				}
			}
		}
	}  // parallel for

}

void process() {

	long **averaged_trajectories = (long**)malloc(sizeof(long*)*g_sites);
	for(int i = 0; i < g_sites; i++) {
		averaged_trajectories[i] = new long[g_steps];
	}
	for(int j = 0; j < g_sites; j++) {
	for(int i = 0; i < g_steps; i++) {
		averaged_trajectories[j][i] = 0;
	}
	}

	#pragma omp parallel for
	for(long traj = 0; traj < g_trajectories; traj++) {

		std::default_random_engine generator;
		generator.seed(traj + std::chrono::system_clock::now().time_since_epoch().count());
		std::uniform_real_distribution<double> distribution(0.0,1.0);

		int pos = 0;
		int since_last_jump = 0;
		for(long it = 0; it < g_steps; it++) {
			double gen = distribution(generator);

			double p = 0;
			for(int m = 0; m < g_sites; m++) {
				if(m == pos) continue;
				if(p + g_step * mol[pos]->m_rates[m][CAP(since_last_jump,g_steps)] >= gen) {
					// jumped
					pos = m;
					since_last_jump = 0;
					break;
				} else {
					// didn't jump'
					p += g_step * mol[pos]->m_rates[m][CAP(since_last_jump,g_steps)];
					since_last_jump++;
				}

				//printf("   %d %d %f\n", m, since_last_jump, mol[pos]->m_rates[m][CAP(since_last_jump,g_steps)]);
			}

			averaged_trajectories[pos][it]++;

		}


	} // parallel for

	for(int i = 0; i < g_steps; i++) {
		printf("%f ", i*g_step);
		for(int k = 0; k < g_sites; k++) {
			printf("%f ", ((double)averaged_trajectories[k][i])/g_trajectories);
		}
		printf("%f", mol[0]->m_rates[1][CAP(i,g_steps)],
				mol[0]->m_gf[CAP(i,g_steps)].real(),  mol[0]->m_gf[CAP(i,g_steps)].imag());
		printf("\n");
	}

	for(int i = 0; i < g_sites; i++) {
		delete [] averaged_trajectories[i];
		averaged_trajectories[i] = NULL;
	}
	free(averaged_trajectories);
}


site::site(int steps, int sites) {
	m_J = new double[sites];
	for(int i = 0; i < sites; i++)
		m_J[i] = 0;

	m_en = 0;
	m_lambda = g_lambda;
	m_tauC = g_tauC;

	dpc cfi, cfimin, hfi, hfimin, gfi, gfimin;
	m_gf = dpcArray(steps);
	const int stepFactor = 10;
	const double step = g_step / stepFactor;

	m_rates = (double**) malloc(sizeof(double*)*sites);
	for(int n = 0; n < sites; n++) {
		m_rates[n] = new double[steps];
		for(int i = 0; i < steps; i++)
			m_rates[n][i] = 0; // initiated when all sites are initiated and energies are known
	}

	cfi = gfunctMatsubara(0, m_lambda, m_tauC, g_temp);
	hfi = 0;
	gfi = 0;
	m_gf[0] = 0;

	for(long i = 1; i < steps * stepFactor; i++) {
		cfimin = cfi;
		hfimin = hfi;
		gfimin = gfi;

		cfi = gfunctMatsubara(step * i, m_lambda, m_tauC, g_temp);

		// Simpson integrator
		hfi = hfimin + step/6.0 * (cfimin + 4.0 * gfunctMatsubara(step * (i - 0.5), m_lambda, m_tauC, g_temp) + cfi);

		// common integrator here
		gfi = gfimin + step/2.0 * (hfimin + hfi);

		if(i % stepFactor == 0)
			m_gf[i / stepFactor] = gfi;
	}

}

site::~site() {
	delete m_J;
}
