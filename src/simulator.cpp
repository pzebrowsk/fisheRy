#include "simulator.h"
using namespace std;

Simulator::Simulator(){
}


Rcpp::DataFrame Simulator::simulate_r(Population &pop, double h, int nyears, bool re_init){
	pop.set_harvestProp(h);
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

vector<double> Simulator::simulate_multi(Population &pop, vector<double> hvec, int nyears, bool re_init){
	Tensor<double> res({10, 4, hvec.size(), nyears});
	Population pop_ref = pop;

	for (int iter = 0; iter < 10; ++iter){
		for (int ih=0; ih<hvec.size(); ++ih){ // loop over control parameter
			pop = pop_ref;
			pop.set_harvestProp(hvec[ih]);
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

