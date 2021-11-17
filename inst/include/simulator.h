#ifndef RFISH_SIMULATOR_H
#define RFISH_SIMULATOR_H

#include "fish.h"
#include "population.h"
#include <vector>

class Simulator{
	private:
	FishParams fish_par;
	PopulationParams pop_par;
	
	public:
	std::vector<double> K;	// age-specific carrying capacity 


	public:
	Simulator(FishParams fish_par, PopulationParams pop_par);
	~Simulator();

	void init(int amax, int nfish);
};

#endif
