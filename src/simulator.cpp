#include "simulator.h"
using namespace std;

Simulator::Simulator(){
}


std::vector<double> Simulator::simulate(Population &pop, double h, int nyears, bool re_init){
	pop.set_harvestProp(h);
	if (re_init) pop.init(1000);
	pop.print_summary();

	std::vector<double> ssb;
	ssb.reserve(nyears);

	for (int i=0; i<nyears; ++i){
		double ssb_now = pop.update();
		ssb.push_back(ssb_now);
	}
	
	return ssb;
}


