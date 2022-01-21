#ifndef RFISH_FISH_H_
#define RFISH_FISH_H_

#include <vector>
#include "functions.h"

struct FishParams {
	// growth 
	double gamma1 = 0.33333;
	double gamma2 = 3.519072971;
	double alpha1 = 3.464594;
	double alpha2 = 0.00156;
	double gsi = 0.464097;

	// maturation
	double pmrn_width = 4.241325;
	double pmrn_envelope = 0.25;
	double pmrn_slope = 1.604122;
	double pmrn_intercept = 18.399575;
	double steepness;  // calculated by constructor

	// reproducttion
	double delta = 1820;
	
	// mortality
	double gamma3 = -1.20565;
	double alpha3 = 0.57792;
	double lref   = 18.25037;
	
	// biology
	double theta = 8.10e-6; // kg/cm 
	double zeta  = 3.01; 
	
	double amax = 15;

	std::vector<double> ma = {0,  0.00000000, 0.00000000, 0.00020000, 0.00270054, 0.05375589, 0.28351881, 0.46168639, 0.57433361, 0.65978050, 0.69259962, 0.65432099, 0.35714286, 1.00000000,
								  1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000, 1.00000000,
								  1.00000000, 1.00000000, 1.00000000, 1.00000000};   // maturation probability with age
	std::vector<double> mai = {0, 1.48936, 0.61272, 0.36018, 0.26108, 0.22300, 0.20447, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000,
								  0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000}; // mortality rate of immature indivuduals
	std::vector<double> mam = {0, 1.48936, 0.61272, 0.36018, 0.26108, 0.22300, 0.20447, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000,
								  0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000, 0.20000}; // mortality rate of mature indivuduals

	double l8 = 196;
	double kappa = 0.0555;
	double a0 = -0.937;

	// debug
	int flag = 0;
	
	FishParams(){
		steepness = fish::maturation_steepness(pmrn_width, pmrn_envelope);
	}
};


class Fish{
	//private:
	//static int nfish = 0;

	public:
	FishParams par;
	
	// state variables
	int age = 1;		// age in years
	double length;  // length in cm

	// physiological variables
	double weight;	// weight in kg

	bool isMature = false;
	bool isAlive = true;

	double t_birth;

	Fish(double tb = 0);
	//Fish(double xb, double tb);

	void set_age(int _a);
	void set_length(double s);
	
	void grow();

	//double maturationProbability();   // return the probability of maturation
	bool matureNow();		// check if fish should mature now, based on maturation probability
	void updateMaturity();	
	
	double naturalMortalityRate();
	//double survivalProbability(double external_mortality_rate, double interval);     // return the probability that this individual survives during the given time interval

	void print();
	void print_line();
	void print_header();
	
	std::vector<double> get_state();
};

#endif
