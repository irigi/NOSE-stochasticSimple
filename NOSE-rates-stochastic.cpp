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

const long g_sites = 3,
		g_steps = 1024*10,
		g_ipl_steps = 1024,
		g_trajectories = 1000000;
site **mol;

const bool g_same_bath = true;
const double g_linearity_criterium = 1e-4;

int main(void) {

	init();
	output_rates();
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

	dpc dummy =
			-mol[n]->m_gf[j]-mol[m]->m_gf[j]					// basic g-functions
			-conj(mol[n]->m_gf[i]) + mol[n]->m_gf[i]									// reorganization energy g-funct
			-mol[n]->m_gf[i-j] + conj(mol[n]->m_gf[i-j])								// reorganization energy g-funct
			-(mol[m]->m_en - mol[n]->m_en)*g_step*j*dpc(0,1)	// energy difference
			;

	// complex exponential is costly, think twice before doing it
	if(dummy.real() > -14)
		dummy = exp(dummy);
	else
		dummy = dpc(0.0,0.0);

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

	if(elderK != acceptorM || k + i >= g_steps)
		return basicRateElement(donorN, acceptorM, donorTimeIndex, transferTimeIndex);

	dpc dummy =
			-mol[n]->m_gf[j]-mol[m]->m_gf[j]					// basic g-functions
			-conj(mol[n]->m_gf[i]) + mol[n]->m_gf[i]									// reorganization energy g-funct
			-mol[n]->m_gf[i-j] + conj(mol[n]->m_gf[i-j])								// reorganization energy g-funct
			-(mol[m]->m_en - mol[n]->m_en)*g_step*j*dpc(0,1)	// energy difference
			;

	dpc extension =
			mol[m]->m_gf[k+i] - mol[m]->m_gf[i] - mol[m]->m_gf[i-j+k] + mol[m]->m_gf[i-j]
			;

	dpc extension2 = -conj(
			mol[m]->m_gf[k+i] - mol[m]->m_gf[i] - mol[m]->m_gf[i-j+k] + mol[m]->m_gf[i-j]
			);

	// complex exponential is costly, think twice before doing it
	if((dummy + extension + extension2).real() > -14)
		dummy = exp(dummy + extension + extension2);
	else
		dummy = dpc(0.0,0.0);

	return 2 * dummy.real() * pow(mol[n]->m_J[m],2.0) * g_step;
}

void init() {
	mol = (site **)malloc(sizeof(site*)*g_sites);
	for(int i = 0; i < g_sites; i++) {
		mol[i] = new site(g_sites);

		mol[i]->m_en = 0;

		if(i+1 < g_sites) 	mol[i]->m_J[(i+1)] = g_J;
		if(i-1 >= 0) 		mol[i]->m_J[(i-1)] = g_J;
	}


	for(int n = 0; n < g_sites; n++) {
		for(int m = 0; m < g_sites; m++) {
			if(m == n) continue;

			#pragma omp parallel for
			for(long i = 1; i < g_steps; i++) {
				for(long j = 1; j <= i; j++) {
					mol[n]->m_rates[m][i] += basicRateElement(n, m, i-j, j);
					//printf("%d %d %d %d %f %f\n", n,m,i,j,mol[n]->m_rates[m][i], dummy.real());
				}
			} // parallel for
		}
	}

	for(int n = 0; n < g_sites; n++) {
		for(int m = 0; m < g_sites; m++) {
			if(m == n) continue;

			for(long k = 0; k < g_ipl_steps; k++) {
				#pragma omp parallel for
				for(long i = 1; i < g_steps; i++) {
					if(k * g_steps/g_ipl_steps > mol[n]->m_gLinInd) {
						// memory is lost, goes back to normal rate
						mol[n]->m_rates3[m][k][i] = mol[n]->m_rates3[m][k-1][i];
					} else {
						if((n > 0 || m > 1) && g_same_bath && mol[0]->m_J[1] != 0.0) {
							// if all baths are identical, we use the fact
							mol[n]->m_rates3[m][k][i] = mol[0]->m_rates3[1][k][i] * pow(mol[n]->m_J[m] / mol[0]->m_J[1],2.0);
						} else {
							for(long j = 1; j <= i; j++) {
								mol[n]->m_rates3[m][k][i] += thirdOrderRateElement(m, n, m, k * g_steps/g_ipl_steps, i-j, j);
							}
						}
					}
				} // parallel for
				printf("k = %ld out of %ld\n", k, g_ipl_steps);
			}
		}
	}

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


site::site(int sites) {
	m_J = new double[sites];
	for(int i = 0; i < sites; i++)
		m_J[i] = 0;

	m_en = 0;
	m_lambda = g_lambda;
	m_tauC = g_tauC;

	dpc cfi, cfimin, hfi, hfimin, gfi, gfimin;
	m_gf = dpcArray(g_steps);
	m_gLinInd = g_steps;
	const int stepFactor = 10;
	const double step = g_step / stepFactor;

	cfi = gfunctMatsubara(0, m_lambda, m_tauC, g_temp);
	hfi = 0;
	gfi = 0;
	m_gf[0] = 0;

	for(long i = 1; i < g_steps * stepFactor; i++) {
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

	m_gLinInd = g_linearity(m_gf);

	m_rates = (double**) malloc(sizeof(double*)*sites);
	for(int n = 0; n < sites; n++) {
		m_rates[n] = new double[g_steps];

		#pragma omp parallel for
		for(int i = 0; i < g_steps; i++)
			m_rates[n][i] = 0; // initiated when all sites are initiated and energies are known
	}

	m_rates3 = (double***) malloc(sizeof(double**)*sites);
	for(int n = 0; n < sites; n++) {
		m_rates3[n] = (double**) malloc(sizeof(double*)*g_ipl_steps);

		#pragma omp parallel for
		for(int k = 0; k < g_ipl_steps; k++) {
			m_rates3[n][k] = new double[g_steps];
			for(int i = 0; i < g_steps; i++)
				m_rates3[n][k][i] = 0; // initiated when all sites are initiated and energies are known
		} // parallel for
	}

}

site::~site() {
	delete m_J;
}

long g_linearity(dpcArray &gf) {
	// find the point of sufficient linearity
	double sum0 = 0;
	for(long i = 3; i < g_steps; i++) {
		double x0 =      i , x1 =      i-1 , x2 =      i-2 , x3 =      i-3 ;
		double y0 = gf[i].real(), y1 = gf[i-1].real(), y2 = gf[i-2].real(), y3 = gf[i-3].real();

		const double a = (-4*(x0 + x1 + x2 + x3)*(y0 + y1 + y2 + y3) +
		     16*(x0*y0 + x1*y1 + x2*y2 + x3*y3))/
		   (-4*pow(x0 + x1 + x2 + x3,2) +
		     16*(pow(x0,2) + pow(x1,2) + pow(x2,2) +
		    		 pow(x3,2)));
		const double b = (y0 + y1 + y2 + y3 +
		     ((x0 + x1 + x2 + x3)*
		        (4*(x0 + x1 + x2 + x3)*
		           (y0 + y1 + y2 + y3) -
		          16*(x0*y0 + x1*y1 + x2*y2 + x3*y3)))/
		      (-4*pow(x0 + x1 + x2 + x3,2) +
		        16*(pow(x0,2) + pow(x1,2) + pow(x2,2) +
		        		pow(x3,2))))/4.;

		double sum = sqrt(pow(a * x0 + b - y0, 2.0) +
			pow(a * x1 + b - y1, 2.0) +
			pow(a * x2 + b - y2, 2.0) +
			pow(a * x3 + b - y3, 2.0));

		if(i == 3)
			sum0 = sum;

		//printf("--> %e\n", sum/sum0);

		if(sum / sum0 < g_linearity_criterium) {
			printf("-->      %ld\n", i);
			return i;
		}
	}

	return g_steps;
}

void output_rates(void) {
	FILE *f;

	f = fopen("rate.dat","w");
	for(int j = 0; j < g_steps; j++) {
		fprintf(f, "%e %e\n", j*g_step, mol[0]->m_rates[1][CAP(j,g_steps)]);
	}
	fclose(f);

	for(int i = 0; i < g_ipl_steps; i += 10 ) {
		char buff[256];
		sprintf(buff, "rate%03d.dat", i);
		f = fopen(buff, "w");

		for(int j = 0; j < g_steps; j++) {
			fprintf(f, "%e %e\n", j*g_step, mol[0]->m_rates3[1][i][CAP(j,g_steps)]);
		}


		fclose(f);
	}
}

double rate3_quad(int donorMol, int acceptorMol, double elderTime, long timeIndex) {
	const double step = g_step * g_steps / g_ipl_steps;

	long i1 = (long)(elderTime/step + 0.5);

	const double
			y0 = mol[donorMol]->m_rates3[acceptorMol][CAP(i1-1,g_ipl_steps)][timeIndex],
			y1 = mol[donorMol]->m_rates3[acceptorMol][CAP(i1  ,g_ipl_steps)][timeIndex],
			y2 = mol[donorMol]->m_rates3[acceptorMol][CAP(i1+1,g_ipl_steps)][timeIndex],

			x0 = CAP(i1-1,g_ipl_steps) * step,
			x1 = CAP(i1  ,g_ipl_steps) * step,
			x2 = CAP(i1+1,g_ipl_steps) * step;


	return	(
				pow(pow(x1 - x2,2)*y0 - pow(x0 - x2,2)*y1,2) -
				2*pow(x0 - x1,2)*(pow(x1 - x2,2)*y0 + pow(x0 - x2,2)*y1)*y2 + pow(x0 - x1,4)*pow(y2,2)
			)/(
				4.*(x0 - x1)*(x0 - x2)*(x1 - x2)*(x2*(y0 - y1) + x0*(y1 - y2) + x1*(y2 - y0))
			);

}
