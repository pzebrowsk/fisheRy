#ifndef RFISH_SIMULATOR_H
#define RFISH_SIMULATOR_H

#include "population.h"
#include <vector>

class Simulator{
	private:
	
	public:

	public:
	Simulator();

	std::vector<double> simulate(Population &pop, double h, int nyears, bool re_init);
};

#endif
