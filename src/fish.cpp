#include "fish.h"
#include <iostream>
#include <cmath>
#include <stdexcept>
using namespace std;


Fish::Fish(double tb){
	t_birth = tb;
}


//Fish::Fish(double xb, double tb){
	//set_age(
	//set_length(xb);
	//t_birth = tb;
	//age = 0;
//}


void Fish::init(double tsb, double temp){
	par.init();
	
	set_age(1);
	
	if (par.growth_model == Model::Joshi23){
		// calc length at age 1
		double tsb_ano = tsb - par.tsbmean;
		double temp_ano = temp - par.Tmean;
		double dl = fish::dl_power(tsb_ano, temp_ano, par.gamma1, par.gamma2, par.alpha1, par.alpha2, par.beta1, par.beta2);
		double l1 = fish::length_juvenile(par.L0, dl, par.gamma1, par.gamma2);

		set_length(l1);
	}
}

void Fish::set_age(int _a){
	age = _a;

	// In Dankel model, length is a direct function of age
	if (par.growth_model == Model::Dankel22){
		length = par.l8*(1-exp(-par.kappa*(age-par.a0)));
		set_length(length);
	}
}


void Fish::set_length(double s){
	length = s;
	
	if (par.growth_model == Model::Dankel22){ 
		weight = par.theta*pow(length, par.zeta);  // Eq. 2
	}
	else{
		weight = fish::weight_fish(length, par.alpha2, par.gamma2)/1000;  
	}
	
}


//double Fish::maturationProbability(){
	//double prob = (age > par.ma.size())? 1 : par.ma[int(age)];
	//return prob;
//}



double Fish::naturalMortalityRate(){
	double rate;
	if (age > par.amax) return 1e20; // FIXME: use inf
	else {
		if (par.use_old_model_mor){
		   	if (isMature) rate = par.mam[age]; 
			else          rate = par.mai[age]; 
			return rate;
		}
		else {
			return fish::natural_mortality(length, par.M0, par.gamma3, par.alpha3, par.Lref);
		}
	}
}


//double Fish::survivalProbability(double mfi, double mfm, double interval){
	//double external_mortality_rate = (isMature)? mfm : mfi;
	//double mortality_rate = naturalMortalityRate() + external_mortality_rate * selectivity();
	//return exp(-mortality_rate*interval);
//}

bool Fish::matureNow(){
	double runif = rand() / double(RAND_MAX);
	
	double maturation_prob;
	if (par.use_old_model_mat) maturation_prob = (age > par.amax)? 1 : par.ma[age];
	else                       maturation_prob = fish::maturation_probability(age, length, par.steepness, par.pmrn_slope, par.pmrn_intercept);
	//cout << "matureNow(): " << age << " " << length << " " << runif << " " << maturation_prob << "\n";
	
	return runif <= maturation_prob;
}


void Fish::updateMaturity(){
	isMature = isMature || matureNow();	// allow short-circuit evalutaion
}


void Fish::grow(double tsb, double temp){
	if (par.growth_model == Model::Dankel22){
		// do nothing. age is incremented by population update
		gsi_effective = par.gsi; // required if new fecundity model is used in combination with old growth model
	}
	else if (par.growth_model == Model::Joshi23){
		double tsb_ano = tsb - par.tsbmean;
		double temp_ano = temp - par.Tmean;
		double dl = fish::dl_power(tsb_ano, temp_ano, par.gamma1, par.gamma2, par.alpha1, par.alpha2, par.beta1, par.beta2);
		double lnew;
		if (isMature){
			lnew = fish::length_adult(length, dl, par.gamma1, par.gamma2, par.gsi);
		}
		else{
			lnew = fish::length_juvenile(length, dl, par.gamma1, par.gamma2);
		}
		gsi_effective = fish::gsi(lnew, length, par.gamma1, par.gamma2, par.alpha1, par.alpha2);
		set_length(lnew);
	}
	else{
		throw std::runtime_error("Invalid growth model specified");
	}
}


double Fish::produceRecruits(){
	if (par.use_old_model_fec) return par.r0 * weight;
	else                       return (par.delta * gsi_effective * weight) * par.s0; // fish::fecundity(length, par.gamma2, par.alpha2, par.delta, gsi_effective) * par.s0; //
}


void Fish::print(){
	cout << "Fish: \n";
	cout << "  age = " << age << "\n";
	cout << "  length = " << length << "\n";
	cout << "  weight = " << weight << "\n";
	cout << "  t_birth = " << t_birth << "\n";
	cout << "  isMature = " << isMature << "\n";
	cout << "----\n";
}

void Fish::print_line(){
	cout << t_birth << "\t" << age << "\t" << isMature << "\t" << isAlive << "\t" << length << "\t" << weight << "\t";
}	

void Fish::print_header(){
	cout << "t_birth" << "\t" << "age" << "\t" << "isMature" << "\t" << "isAlive" << "\t" << "length" << "\t" << "weight" << "\t";
}	

std::vector<double> Fish::get_state(){
	return {t_birth, age, isMature, isAlive, length, weight};
}




