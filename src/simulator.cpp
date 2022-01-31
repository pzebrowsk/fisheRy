#include "simulator.h"
#include <chrono> // microseconds
#include <thread> // sleep_for
using namespace std;

Simulator::Simulator(Fish f) : noFishingPop(f) {
}

void Simulator::setNaturalPopulation(Population &pop){
	noFishingPop = pop;
}


Rcpp::DataFrame Simulator::simulate_r(Population &pop, double lf, double h, int nyears, bool re_init){
	noFishingPop.set_harvestProp(h);
	noFishingPop.set_minSizeLimit(lf);
	double K = noFishingPop.fishableBiomass();

	pop.K_fishableBiomass = K;
	pop.set_harvestProp(h);
	pop.set_minSizeLimit(lf);
	if (re_init) pop.init(1000);
	pop.print_summary();

	std::vector<double> ssb;
	ssb.reserve(nyears);
	
	std::vector<double> yield;
	yield.reserve(nyears);
	
	std::vector<double> emp_sea;
	emp_sea.reserve(nyears);
	std::vector<double> emp_shore;
	emp_shore.reserve(nyears);
	
	std::vector<double> psea;
	psea.reserve(nyears);
	std::vector<double> pshr;
	pshr.reserve(nyears);
	
	Rcpp::DataFrame df = Rcpp::DataFrame::create();

	for (int i=0; i<nyears; ++i){
		std::vector<double> state_now = pop.update();
		ssb.push_back(state_now[0]);
		yield.push_back(state_now[1]);
		emp_sea.push_back(state_now[2]);
		emp_shore.push_back(state_now[3]);
		psea.push_back(state_now[4]);
		pshr.push_back(state_now[5]);
	}
	
	df.push_back(ssb, "ssb");
	df.push_back(yield, "yield");
	df.push_back(emp_sea, "employment.sea");
	df.push_back(emp_shore, "employment.shore");
	df.push_back(psea, "profit.sea");
	df.push_back(pshr, "profit.shore");

	return df;
}

//double Simulator::calcK(Population &pop, double lmin, double h){
//	pop.set_harvestProp(h);
//	pop.set_minSizeLimit(lmin);
//	double K = 0;
//	for (auto f : pop.fishes) K += f.weight * pop.par.n * pop.selectivity(f.length);
//	return K;
//}


vector<double> Simulator::simulate_multi(Population &pop, vector<double> hvec, int nyears, bool re_init){
	int niters = 1;
	Tensor<double> res({niters, 4, hvec.size(), nyears});
	Population pop_ref = pop;

	for (int iter = 0; iter < niters; ++iter){
		for (int ih=0; ih<hvec.size(); ++ih){ // loop over control parameter
			pop = pop_ref;
			pop.set_harvestProp(hvec[ih]);
			
			noFishingPop.set_harvestProp(hvec[ih]);
			noFishingPop.set_minSizeLimit(pop.par.lf50);
			double K = noFishingPop.fishableBiomass();
			pop.K_fishableBiomass = K;
	
			if (re_init) pop.init(1000);
			pop.print_summary();
		
			for (int t=0; t<nyears; ++t){
				std::vector<double> state_now = pop.update();
				
				res({iter, 0, ih, t}) = state_now[0];				// ssb
				res({iter, 1, ih, t}) = state_now[1];				// yield
				res({iter, 2, ih, t}) = state_now[2]+state_now[3];				// employment
				res({iter, 3, ih, t}) = state_now[4]+state_now[5];	// total profit
			}
		}
	}
	res.print();
	return res.avg_dim(3).vec;	// average over iterations

}


vector<double> Simulator::max_avg_utils(vector<int> dims, vector<double> data){
	Tensor<double> res(dims);
	res.vec = data;		// res is {u, c1, t}

	Tensor<double> res2 = res.avg_dim(0).max_dim(0);	// avg over t, then max over c1
	res.transform(2, std::divides<double>(), res2.vec); // divide u dimension by res2

	return res.avg_dim(0).vec;
}


vector<double> Simulator::stakeholder_satisfaction(vector<int> dims, vector<double> data){
	Tensor<double> res(dims);
	res.vec = data;		// res is {u, c1, t}

	Tensor<double> res2 = res.avg_dim(0).max_dim(0);	// avg over t, then max over c1
	res.transform(2, std::divides<double>(), res2.vec); // divide u dimension by res2

	Tensor<double> sp({5,4});	// spvec is {s, u}
	//        ssb yield emp  profit  
	sp.vec = {0.0, 0.3, 0.0, 0.7,	// industrial
			  0.3, 0.5, 0.1, 0.1,	// artisanal
			  0.3, 0.2, 0.5, 0.0,	// employment-maximizing policymakers
			  0.2, 0.2, 0.0, 0.6,	// profit-maximizing policymakers
			  0.5, 0.1, 0.2, 0.2	// conservationists
			 };

	sp.print();

	Tensor<double> Ssucy = sp.repeat_inner(res.dim[1]).repeat_inner(res.dim[2]) * res.repeat_outer(sp.dim[0]);
	//                        ^ {s, u, c}              ^ {s, u, c, y}			      ^ {s, u, c, y}

	Tensor<double> Sscy = Ssucy.accumulate(0.0, 2, std::plus<double>());	// aggregate along u dim to get {s, c, y}

	Tensor<double> Ssc = Sscy.avg_dim(0);
	//                        ^ {s, c}
	Ssc.transform(1, std::divides<double>(), Ssc.max_dim(0).vec);
	//			  ^ s                            ^ {s}
	
	return Ssc.vec;

}



vector<double> Simulator::simulate_multi_2d(Population &pop, vector<double> lminvec, vector<double> hvec, int nyears, bool re_init){
	int niters = 1;
	Tensor<double> res({niters, 4, lminvec.size(), hvec.size(), nyears});
	Population pop_ref = pop;

	for (int iter = 0; iter < niters; ++iter){
		for (int il=0; il<lminvec.size(); ++il){ // loop over control parameter 2
			for (int ih=0; ih<hvec.size(); ++ih){ // loop over control parameter 1
				pop = pop_ref;
				pop.set_harvestProp(hvec[ih]);
				pop.set_minSizeLimit(lminvec[il]);
				
				noFishingPop.set_harvestProp(hvec[ih]);
				noFishingPop.set_minSizeLimit(lminvec[il]);
				double K = noFishingPop.fishableBiomass();
				pop.K_fishableBiomass = K;

				if (re_init) pop.init(1000);
				pop.print_summary();
			
				for (int t=0; t<nyears; ++t){
					std::vector<double> state_now = pop.update();
					
					res({iter, 0, il, ih, t}) = state_now[0];               // ssb
					res({iter, 1, il, ih, t}) = state_now[1];               // yield
					res({iter, 2, il, ih, t}) = state_now[2]+state_now[3];  // employment
					res({iter, 3, il, ih, t}) = state_now[4]+state_now[5];  // total profit
				}
			}
			std::this_thread::sleep_for(std::chrono::microseconds(100));
		}
	}
	//res.print();
	
	return res.avg_dim(4).vec;	// average over iterations

}


vector<double> Simulator::max_avg_utils_2d(vector<int> dims, vector<double> data){
	Tensor<double> res(dims);
	res.vec = data;		// res is {u, c2, c1, t}

	Tensor<double> res2 = res.avg_dim(0).max_dim(0).max_dim(0);	// avg over t, then max over c1, then max over c2
	res.transform(3, std::divides<double>(), res2.vec); // divide u dimension by res2

	return res.avg_dim(0).vec;
}


vector<double> Simulator::stakeholder_satisfaction_2d(vector<int> dims, vector<double> data){
	Tensor<double> res(dims);
	res.vec = data;		// res is {u, c2, c1, t}

	Tensor<double> res2 = res.avg_dim(0).max_dim(0).max_dim(0);  // avg over t, then max over c1, then max over c2
	res.transform(3, std::divides<double>(), res2.vec); // divide u dimension by res2

	Tensor<double> sp({5,4});	// spvec is {s, u}
	//        ssb yield emp  profit  
	sp.vec = {0.0, 0.3, 0.0, 0.7,	// industrial
			  0.3, 0.5, 0.1, 0.1,	// artisanal
			  0.3, 0.2, 0.5, 0.0,	// employment-maximizing policymakers
			  0.2, 0.2, 0.0, 0.6,	// profit-maximizing policymakers
			  0.5, 0.1, 0.2, 0.2	// conservationists
			 };

	sp.print();

	Tensor<double> Ssucy = sp.repeat_inner(res.dim[1]).repeat_inner(res.dim[2]).repeat_inner(res.dim[3]) * res.repeat_outer(sp.dim[0]);
	//                        ^ {s, u, c2}             ^ {s, u, c2, c1}         ^ {s, u, c2, c1, y}            ^ {s, u, c2, c1, y}

	Tensor<double> Sscy = Ssucy.accumulate(0.0, 3, std::plus<double>());	// aggregate along u dim to get {s, c2, c1, y}

	Tensor<double> Ssc = Sscy.avg_dim(0);
	//                        ^ {s, c2, c1}
	Ssc.transform(2, std::divides<double>(), Ssc.max_dim(0).max_dim(0).vec);
	//			  ^ s                            ^ {s, c2}  ^ {s}
	
	return Ssc.vec;

}


