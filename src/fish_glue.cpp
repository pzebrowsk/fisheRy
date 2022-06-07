#include <Rcpp.h>
using namespace Rcpp;

#include "fish.h"
#include "functions.h"

RCPP_EXPOSED_CLASS(Fish);
RCPP_EXPOSED_CLASS(FishParams);

RCPP_EXPOSED_ENUM_NODECL(Model)

RCPP_MODULE(fish_module) {
	
//	function("init_length", &fish::init_length);
	function("length_juvenile", &fish::length_juvenile);
	function("length_adult", &fish::length_adult);
	function("maturation_steepness", &fish::maturation_steepness);
	function("maturation_probability", &fish::maturation_probability);
	function("weight_fish", &fish::weight_fish);
	function("fecundity", &fish::fecundity);
	function("gsi", &fish::gsi);
	function("natural_mortality", &fish::natural_mortality);
	function("survival_probability", &fish::survival_probability);
	
	class_ <FishParams>("FishParams")
		.constructor()
		.field("flag", &FishParams::flag)
//		.field("Bhalf_growth", &FishParams::Bhalf_growth)
		.field("c", &FishParams::c)
		.field("beta1", &FishParams::beta1)
		.field("beta2", &FishParams::beta2)
		.field("s0", &FishParams::s0)
		.field("growth_model", &FishParams::growth_model)
		.field("use_old_model_mat", &FishParams::use_old_model_mat)
		.field("use_old_model_mor", &FishParams::use_old_model_mor)
		.field("use_old_model_fec", &FishParams::use_old_model_fec)
	;

	class_ <Fish>("Fish")
		.constructor()
//		.constructor<double>()
		
		.field("age", &Fish::age)
		.field("length", &Fish::length)
		.field("weight", &Fish::weight)
		.field_readonly("t_birth", &Fish::t_birth)
		.field("par", &Fish::par)

		.method("print", &Fish::print)
		.method("print_line", &Fish::print_line)
		.method("print_header", &Fish::print_header)

		.method("set_age", &Fish::set_age)
		//.method("set_length", &Fish::set_length)

		.method("init", &Fish::init)
		.method("matureNow", &Fish::matureNow)
		.method("updateMaturity", &Fish::updateMaturity)
		.method("grow", &Fish::grow)
		
		.method("naturalMortalityRate", &Fish::naturalMortalityRate)
		//.method("survivalProbability", &Fish::survivalProbability)

		.method("get_state", &Fish::get_state)

	;
}	
	
#include "population.h"

RCPP_EXPOSED_CLASS(PopulationParams);


////RCPP_EXPOSED_AS(Population);
RCPP_MODULE(population_module){
	class_ <PopulationParams>("PopulationParams")
		.constructor()
		.field("n", &PopulationParams::n)
		.field("Bhalf", &PopulationParams::Bhalf)
		.field_readonly("h", &PopulationParams::h)
		.field_readonly("lf50", &PopulationParams::lf50)
//		.field("mort_fishing_mature", &PopulationParams::mort_fishing_mature) 
//		.field("mort_fishing_immature", &PopulationParams::mort_fishing_immature) 
		.field("dsea", &PopulationParams::dsea)
	;
	
	class_ <Population>("Population")
		.constructor<Fish>()
		.field("par", &Population::par)
		.field("verbose", &Population::verbose)
		.field("K", &Population::K_fishableBiomass)
		
		.method("set_superFishSize", &Population::set_superFishSize) 
		
		.method("set_harvestProp", &Population::set_harvestProp) 
		.method("set_minSizeLimit", &Population::set_minSizeLimit) 

		.method("selectivity", &Population::selectivity) 
		.method("init", &Population::init) 
		.method("update", &Population::update)
		.method("calcSSB", &Population::calcSSB)
		.method("fishableBiomass", &Population::fishableBiomass)

		.method("noFishingEquilibriate", &Population::noFishingEquilibriate)

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
	.constructor<Fish>()

	//.field("noFishingPop", &Simulator::noFishingPop)

	.method("setNaturalPopulation", &Simulator::setNaturalPopulation)
	.method("equilibriateNaturalPopulation", &Simulator::equilibriateNaturalPopulation)
	
	.method("simulate", &Simulator::simulate_r)
    
	.method("simulate_multi", &Simulator::simulate_multi_r)
    .method("max_avg_utils", &Simulator::max_avg_utils)
    .method("stakeholder_satisfaction", &Simulator::stakeholder_satisfaction)
	
	.method("simulate_multi_2d", &Simulator::simulate_multi_2d_r)
    .method("max_avg_utils_2d", &Simulator::max_avg_utils_2d)
    .method("stakeholder_satisfaction_2d", &Simulator::stakeholder_satisfaction_2d)
	;
}




