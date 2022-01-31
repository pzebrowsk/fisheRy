#ifndef RFISH_FUNCTIONS_H
#define RFISH_FUNCTIONS_H

#include <cmath>
#include <algorithm>
#include <iostream>

#include <Rcpp.h>

namespace fish {

//// [[Rcpp::export]]
//inline double init_length(double gamma_1, double gamma_2, double alpha_1, double alpha_2){
////	double s = gamma_1 * alpha_1 / alpha_2;
////	return pow(s, 1 / gamma_2);
//	return 19.18;
//}


// [[Rcpp::export]]
inline double length_juvenile(double init_length, double gamma_1, double gamma_2, double alpha_1, double alpha_2){
	double t1 = pow(init_length, gamma_1 * gamma_2);
	double t2 = gamma_1 * alpha_1 * pow(alpha_2, -gamma_1);
	return pow(t1 + t2, 1 / (gamma_1 * gamma_2));
}


// [[Rcpp::export]]
inline double length_adult(double init_length, double gamma_1, double gamma_2, double alpha_1, double alpha_2, double gsi){
	double l = length_juvenile(init_length, gamma_1, gamma_2, alpha_1, alpha_2) / pow(1 + gamma_1 * gsi, 1 / (gamma_1 * gamma_2));
	return std::max(l, init_length); 

}


// [[Rcpp::export]]
inline double maturation_steepness(double pmrn_width, double pmrn_envelope){
	return pmrn_width / (log((1- pmrn_envelope) / pmrn_envelope) - log(pmrn_envelope / (1 - pmrn_envelope)));
}


// [[Rcpp::export]]
inline double maturation_probability(double age, double body_length, double steepness, double pmrn_slope, double pmrn_intercept){
	double pmrn_midpoint = pmrn_slope * age + pmrn_intercept;
	double p = 1.0 / (1.0 + exp(- (body_length - pmrn_midpoint) / steepness));
	//std::cout << "p_mat: age = " << age << ", length = " << body_length << ", slope = " << pmrn_slope << ", int = " << pmrn_intercept << ", mid = " << pmrn_midpoint << ", steepness = " << steepness << ", p = " << p << "\n";
	return p;  
}


// [[Rcpp::export]]
inline double weight_fish(double body_length, double alpha_2, double gamma_2){
	return alpha_2 * pow(body_length, gamma_2);
}


// [[Rcpp::export]]
inline double fecundity(double body_length, double gamma_2, double alpha_2, double delta, double gsi) {
	return delta * weight_fish(body_length, alpha_2, gamma_2) * gsi;
}


// effective GDI given initial and final length
// [[Rcpp::export]]
inline double gsi(double body_length, double init_body_length, double gamma_1, double gamma_2, double alpha_1, double alpha_2){
	double t1 = pow(init_body_length, gamma_1 * gamma_2);
	double t2 = pow(body_length, gamma_1 * gamma_2);
	double t3 = gamma_1 * alpha_1 * pow(alpha_2, -gamma_1);
	double d  = gamma_1 * pow(body_length, gamma_1 * gamma_2);
	
	return (t1 - t2 + t3)/d;
}


// [[Rcpp::export]]
inline double natural_mortality(double body_length, double gamma_3, double alpha_3, double body_length_ref){
	return alpha_3 * pow(body_length / body_length_ref, gamma_3);
}


// [[Rcpp::export]]
inline double survival_probability(double body_length, double gamma_3, double alpha_3, double body_length_ref, double fishing_mortality, double dt = 1){
	double total_mortality = natural_mortality(body_length, gamma_3, alpha_3, body_length_ref) + fishing_mortality;
	return exp(-total_mortality * dt);
}


} // namespace fish

#endif





