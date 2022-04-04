#ifndef RFISH_SIMULATOR_H
#define RFISH_SIMULATOR_H

#include "population.h"
#include <vector>
#include <Rcpp.h>
#include "tensor.h"

class Simulator{
	private:
	Population noFishingPop;
	
	public:

	public:
	Simulator(Fish f);

	void setNaturalPopulation(Population & pop); 

	Tensor<double> simulate_multi(Population &pop, std::vector<double> hvec, int nyears, bool re_init);
	std::vector<double> max_avg_utils(std::vector<int> dims, std::vector<double> data);
	std::vector<double> stakeholder_satisfaction(std::vector<int> dims, std::vector<double> data);
	
//	double calcK(Population &pop, double lmin, double h);
//	std::vector<double> calcK_2d(Population &pop, std::vector<double> lminvec, std::vector<double> hvec);

	std::vector<double> simulate_multi_2d(Population &pop, std::vector<double> lminvec, std::vector<double> hvec, int nyears, bool re_init);
	std::vector<double> max_avg_utils_2d(std::vector<int> dims, std::vector<double> data);
	std::vector<double> stakeholder_satisfaction_2d(std::vector<int> dims, std::vector<double> data);
	
	Rcpp::DataFrame simulate_r(Population &pop, double lf, double h, int nyears, bool re_init);

	Rcpp::NumericVector simulate_multi_r(Population &pop, std::vector<double> hvec, int nyears, bool re_init);

};

#endif
