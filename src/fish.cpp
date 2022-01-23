#include "fish.h"
#include <iostream>
#include <cmath>
using namespace std;

Fish::Fish(double tb){

	t_birth = tb;
	
	set_age(1);
	
	double l0 = fish::init_length(par.gamma1, par.gamma2, par.alpha1, par.alpha2);
	set_length(l0);
	
	//cout << "Creating fish: l0 = " << l0 << ", steepness = " << par.steepness << "\n"; 
}

//Fish::Fish(double xb, double tb){
	//set_age(
	//set_length(xb);
	//t_birth = tb;
	//age = 0;
//}


void Fish::set_age(int _a){
	age = _a;

//	// in an explicity age-length model, this will go.
//	length = par.l8*(1-exp(-par.kappa*(age-par.a0)));
//	set_length(length);
}


void Fish::set_length(double s){
	length = s;
	weight = fish::weight_fish(length, par.alpha2, par.gamma2);  
}


//double Fish::maturationProbability(){
	//double prob = (age > par.ma.size())? 1 : par.ma[int(age)];
	//return prob;
//}



double Fish::naturalMortalityRate(){
	double rate;
	if (age > par.amax) return 1e20;
	else return fish::natural_mortality(length, par.gamma3, par.alpha3, par.Lref);
	  
//   	if (isMature) rate = (age > par.amax)? 1e20 : par.mam[age]; // FIXME: use inf
//	else          rate = (age > par.amax)? 1e20 : par.mai[age]; // FIXME: use inf
//	return rate;
}


//double Fish::survivalProbability(double mfi, double mfm, double interval){
	//double external_mortality_rate = (isMature)? mfm : mfi;
	//double mortality_rate = naturalMortalityRate() + external_mortality_rate * selectivity();
	//return exp(-mortality_rate*interval);
//}

bool Fish::matureNow(){
	double runif = rand() / double(RAND_MAX);
//	double maturation_prob = (age > par.amax)? 1 : par.ma[age];
	double maturation_prob = fish::maturation_probability(age, length, par.steepness, par.pmrn_slope, par.pmrn_intercept);
	assert(par.steepness > 0);
	//cout << "matureNow(): " << age << " " << length << " " << runif << " " << maturation_prob << "\n";
	return runif <= maturation_prob;
}

void Fish::updateMaturity(){
	isMature = isMature || matureNow();	// allow short-circuit evalutaion
}

void Fish::grow(){
	double lnew;
	if (isMature){
		lnew = fish::length_adult(length, par.gamma1, par.gamma2, par.alpha1, par.alpha2, par.gsi);
	}
	else{
		lnew = fish::length_juvenile(length, par.gamma1, par.gamma2, par.alpha1, par.alpha2);
	}
	gsi_effective = fish::gsi(lnew, length, par.gamma1, par.gamma2, par.alpha1, par.alpha2);
	set_length(lnew);
}

double Fish::produceEggs(){
	return par.delta * gsi_effective * weight;
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




