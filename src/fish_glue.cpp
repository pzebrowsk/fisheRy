#include <Rcpp.h>
using namespace Rcpp;

#include "fish.h"

RCPP_EXPOSED_CLASS(Fish);
RCPP_EXPOSED_CLASS(FishParams);

RCPP_MODULE(fish_module) {
	class_ <FishParams>("FishParams")
		.constructor()
		.field("flag", &FishParams::flag)
	;

	class_<Fish>("Fish")
		.constructor()
		//.constructor<double,double>()
		
		.field("age", &Fish::age)
		.field("length", &Fish::length)
		.field("weight", &Fish::weight)
		.field_readonly("t_birth", &Fish::t_birth)
		.field("par", &Fish::par)

		.method("print", &Fish::print)
		.method("set_age", &Fish::set_age)
		//.method("set_length", &Fish::set_length)

		.method("matureNow", &Fish::matureNow)
		.method("updateMaturity", &Fish::updateMaturity)
		
		.method("naturalMortalityRate", &Fish::naturalMortalityRate)
		//.method("survivalProbability", &Fish::survivalProbability)
	;
}	
	
#include "population.h"

RCPP_EXPOSED_CLASS(PopulationParams);


////RCPP_EXPOSED_AS(Population);
RCPP_MODULE(population_module){
	class_ <PopulationParams>("PopulationParams")
		.constructor()
		.field("n", &PopulationParams::n)
		.field("h", &PopulationParams::h)
		.field("mort_fishing_mature", &PopulationParams::mort_fishing_mature) 
		.field("mort_fishing_immature", &PopulationParams::mort_fishing_immature) 
	;
	
	class_ <Population>("Population")
		.constructor<Fish>()
		.method("set_harvestProp", &Population::set_harvestProp) 

		.method("calcK", &Population::calcK) 
		.method("selectivity", &Population::selectivity) 
		.method("init", &Population::init) 
		.method("update", &Population::update)
		.method("calcSSB", &Population::calcSSB)

		.method("get_state", &Population::get_state)
		.method("summarize", &Population::summarize)
		.method("print_summary", &Population::print_summary)
		.method("nfish", &Population::nfish)
	;
}


#include "simulator.h"
RCPP_EXPOSED_CLASS(Population);

RCPP_MODULE(simulator_module){
	class_ <Simulator>("Simulator")
		.constructor()
    .method("simulate", &Simulator::simulate_r)
    .method("simulate_multi", &Simulator::simulate_multi)
    .method("max_avg_utils", &Simulator::max_avg_utils)
    .method("stakeholder_satisfaction", &Simulator::stakeholder_satisfaction)
	;
}




