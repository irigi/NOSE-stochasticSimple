/*
 * NOSE-rates-stochastic.h
 *
 *  Created on: Apr 11, 2014
 *      Author: Jan Olsina
 */

#ifndef NOSE_RATES_STOCHASTIC_H_
#define NOSE_RATES_STOCHASTIC_H_

#include <math.h>
#include <complex>
#include <valarray>

#define KB_CM                     0.6950387
#define ENERGY_INTERNAL_TO_CM  5308.8415534785
#define KB    				(KB_CM/ENERGY_INTERNAL_TO_CM)
#define CAP(X,Y) (X >= Y ? Y-1 : X)

typedef std::complex<double> dpc;
typedef std::valarray<dpc> dpcArray;

dpc gfunctMatsubara(double t, double lambda, double tauC, double temp);
void init();
void process();

class site {
public:
	site(int steps, int sites);
	~site();

	dpcArray m_gf;
	double m_en, *m_J, m_lambda, m_tauC;
	double **m_rates; // rates to all others
	double ***m_rates3;
};



#endif /* NOSE_RATES_STOCHASTIC_H_ */
