#include <iostream>
#include <cassert>
#include <vector>
#include <functional>
#include <algorithm>

#include <numeric>

#include <tensor.h>

using namespace std;

bool equals(double x1, double x2, double tol=1e-6){
	return (abs(x1-x2) < tol);
}

template <class T>
bool equals(vector<T> x1, vector<T> x2, double tol=1e-6){
	if (x1.size()!= x2.size()) return false;
	for (int i=0; i<x1.size(); ++i) if (abs(x1[i]-x2[i]) > tol) return false;
	return true;
}


int main(){
	


}

