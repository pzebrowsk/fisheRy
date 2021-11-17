#ifndef RFISH_FISH_H_
#define RFISH_FISH_H_

#include <vector>

struct FishParams {
	// biology
	double theta = 8.10e-6; // kg/cm 
	double zeta  = 3.01; 
	
	double amax = 30;

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
	
	//double maturationProbability();   // return the probability of maturation
	bool matureNow();		// check if fish should mature now, based on maturation probability
	void updateMaturity();	
	
	double naturalMortalityRate();
	//double survivalProbability(double external_mortality_rate, double interval);     // return the probability that this individual survives during the given time interval

	void print();
	void print_line();
	static void print_header();
};

#endif
