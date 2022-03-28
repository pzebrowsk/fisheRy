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
	double Bhalf = 365426284.0; ///5000;	// Half saturation constant of recruitment

//	double s0 = 0.1126797; //0.11;          // Egg survival propbability

	// management / fishing selectivity
	double sf = 0.1222;	// steepness of selectivity curve
	double lf50 = 45; //61.4806;  // threshold fish length

	// environmental stochasticity
	double sigmaf = 0.4858775;

	// effort dynamics and employment
	double q = 2.83e-6;		// scaling parameter relating to catchability and density
	double dsea = 0.054;	// Required Person-years per vessel day
	double dmax = 30000e20;	// max available person-years
	double dshr = 0.000004;	// FTE/kg
	double b = 0.75;		// density dependence

	// revenue and profit 
	double price_sea = 13.13;		// landing price NOK/kg
	double price_shore = 17.0;		// selling price NOK/kg

	double salary_sea = 1078000;			// employment cost sea NOK/FTE
	double salary_shore = 348000;			// employment cost shore NOK/FTE
	double fixed_costs_sea = 351123000;	// fixed costs sea NOK (= average per unit * #units)
	double fixed_costs_shore = 1032468000;	// fixed costs shore NOK
	double variable_costs_sea = 65000; 		// variable costs NOK/vessel day
	double scale_catch = 0.356; //0.53; 		// percentage of total codfish catch that is cod
	
	double h = 0;

	double n = 1e5;	// superfish size

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
	double K_fishableBiomass = 0;	
	
	public:
	double current_year = 1;

	Fish proto_fish;	//  a copy of this fish is always used to initialize new fish in population
	std::vector<Fish> fishes;
	
	PopulationParams par;	
	
	public:
	Population(Fish f);
	
	void set_superFishSize(double _n);

//	void calc_athresh();
	void set_harvestProp(double _h);
	void set_minSizeLimit(double _lf50);
	void init(int n);	// initialize population with n individuals

	std::vector<double> noFishingEquilibriate();	

	double calcSSB();	
	double selectivity(double len);
	double fishableBiomass();

	double calcRealizedFishingMortality();
	double effort(double Nr, double F);

	std::vector<double> update();

	int nfish();
	void summarize();
	void print_summary();
	Rcpp::DataFrame get_state();

};


#endif
