#include "fish.h"
#include <iostream>
#include <cmath>
using namespace std;

Fish::Fish(double tb){
	t_birth = tb;
	set_age(age);
}

//Fish::Fish(double xb, double tb){
	//set_age(
	//set_length(xb);
	//t_birth = tb;
	//age = 0;
//}


void Fish::set_age(int _a){
	age = _a;

	// in an explicity age-length model, this will go.
	length = par.l8*(1-exp(-par.kappa*(age-par.a0)));
	set_length(length);
}


void Fish::set_length(double s){
	length = s;
	weight = par.theta*pow(length, par.zeta);  // Eq. 2
}


//double Fish::maturationProbability(){
	//double prob = (age > par.ma.size())? 1 : par.ma[int(age)];
	//return prob;
//}



double Fish::naturalMortalityRate(){
	double rate;
   	if (isMature) rate = (age > par.amax)? 1e20 : par.mam[age]; // FIXME: use inf
	else          rate = (age > par.amax)? 1e20 : par.mai[age]; // FIXME: use inf
	return rate;
}


//double Fish::survivalProbability(double mfi, double mfm, double interval){
	//double external_mortality_rate = (isMature)? mfm : mfi;
	//double mortality_rate = naturalMortalityRate() + external_mortality_rate * selectivity();
	//return exp(-mortality_rate*interval);
//}

bool Fish::matureNow(){
	double runif = rand() / double(RAND_MAX); 
	double maturation_prob = (age > par.amax)? 1 : par.ma[age];
	return runif <= maturation_prob;
}

void Fish::updateMaturity(){
	isMature = isMature || matureNow();	// allow short-circuit evalutaion
}



void Fish::print(){
	cout << "Fish: \n";
	cout << "  age = " << age << "\n";
	cout << "  length = " << length << "\n";
	cout << "  weight = " << weight << "\n";
	cout << "  t_birth = " << t_birth << "\n";
	cout << "----\n";
}

void Fish::print_line(){
	cout << t_birth << "\t" << age << "\t" << isMature << "\t" << isAlive << "\t" << length << "\t" << weight << "\t";
}	

void Fish::print_header(){
	cout << "t_birth" << "\t" << "age" << "\t" << "isMature" << "\t" << "isAlive" << "\t" << "length" << "\t" << "weight" << "\t";
}	


