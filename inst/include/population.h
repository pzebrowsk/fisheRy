#ifndef RFISH_POPULATION_H_
#define RFISH_POPULATION_H_

#include "fish.h"
#include <vector>

#include <Rcpp.h>

class PopulationParams {	
	public:
	// reproduction
	double r0 = 21.77072;		// recruitment rate per kg SSB 
	double rmax = 1e20;
	double Bhalf = 365426284;	// Half saturation constant of recruitment

	// management / fishing selectivity
	double sf = 0.1222;	// steepness of selectivity curve
	double lf50 = 61.4806;  // threshold fish length
	
	double b = 0.75;		// density dependence

	// environmental stochasticity
	double sigmaf = 0.4858775;

	// effort dynamics
	double q = 2.83e-6;		// scaling parameter relating to catchability and density
	double dsea = 0.054;	// Required Person-hours per vessel day
	double dmax = 30000;	// max available person-hours

	double h = 0;

	double n = 1e6;	// superfish size

	// ***
	// calculated variables
	double mort_fishing_mature = 0; 
	double mort_fishing_immature = 0; 
	int a_thresh;	// threshold age over which fishing selectivity is > 0.5

	//// function to init
	//PopulationParams(double _h){
	//    h = _h;
	//    mort_fishing_mature = -log(1-h);
	//    mort_fishing_immature = -log(1-h);
	//}
};


// TODO: Enable systematic storing and operating population distributions
class PopulationSummary{
	std::vector<double> vage, vfreq, vlen, vmat;
};

class Population{
	private:
	std::vector<double> vage, vfreq, vlen, vmat;
	std::vector<double> carrying_capacity;

	public:
	double current_year = 1;

	Fish proto_fish;	//  a copy of this fish is always used to initialize new fish in population
	std::vector<Fish> fishes;
	
	PopulationParams par;	
	
	public:
	Population(Fish f);
	
	void set_harvestProp(double _h);
	void init(int n);	// initialize population with n individuals

	std::vector<double> calcK();	

	double calcSSB();	
	double selectivity(double len);

	double calcRealizedFishingMortality();
	double effort(double Nr, double F);

	double update();

	int nfish();
	void summarize();
	void print_summary();
	Rcpp::DataFrame get_state();

};


#endif
