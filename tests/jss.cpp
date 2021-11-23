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

	Tensor<double> u({2,3,5});
	u.fill_sequence();
	u.print();

	for (int i=0; i<30; ++i){
		auto v = u.index(i);
		cout << i << ": ";
		for (auto x : v) cout << x << " ";
		cout << "\n";
	}


	double d;
	d = u.accumulate_dim(0, u.location({0,0,0}), 0, std::plus<double>());
	cout << "accum0 (0,0,0) = " << d << "\n"; // expected 10
	if (!equals(d,10)) return 1;

	d= u.accumulate_dim(0, u.location(1,1,0),   0, std::plus<double>());
	cout << "accum0 (1,1,0) = " << d << "\n"; // expected 110	// template version
	if (!equals(d,110)) return 1;
	
	d= u.accumulate_dim(0, u.location(1,0,3),   1, std::plus<double>());
	cout << "accum1 (1,0,3) = " << d << "\n"; // expected 69 	// template version
	if (!equals(d,69)) return 1;
	
	d= u.accumulate_dim(0, u.location({0,1,4}), 2, std::plus<double>());
	cout << "accum1 (0,1,4) = " << d << "\n"; // expected 33
	if (!equals(d,33)) return 1;


	vector<int> x = u.plane(0);
	vector<int> expected;
	expected = {0,5,10, 15,20,25};
	cout << "starts dim 0: "; for (auto xx : x) cout << xx << " "; cout << "\n";  
	if (!equals(x, expected)) return 1;

	x = u.plane(1);
	expected = {0,1,2,3,4,15,16,17,18,19};
	cout << "starts dim 1: "; for (auto xx : x) cout << xx << " "; cout << "\n";	
	if (!equals(x, expected)) return 1;

	x = u.plane(1,1);
	expected = {5,6,7,8,9,20,21,22,23,24};
	cout << "starts dim 1, off 1: "; for (auto xx : x) cout << xx << " "; cout << "\n";	
	if (!equals(x, expected)) return 1;

	x = u.plane(2);
	expected = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14};
	cout << "starts dim 2: "; for (auto xx : x) cout << xx << " "; cout << "\n";	
	if (!equals(x, expected)) return 1;

	x = u.plane(2,1);
	expected = {15,16,17,18,19,20,21,22,23,24,25,26,27,28,29};
	cout << "starts dim 2, off 1: "; for (auto xx : x) cout << xx << " "; cout << "\n";	
	if (!equals(x, expected)) return 1;

	
	vector<double> expected1;
	Tensor<double> v = u.accumulate(0, 0, std::plus<double>());
	expected1 = {10,35,60,  
	            85,110,135};
	v.print();
	if (!equals(v.vec, expected1)) return 1;
	
	Tensor<double> v1 = u.accumulate(0, 1, std::plus<double>());
	expected1 = {15, 18, 21, 24, 27, 
                60, 63, 66, 69, 72};
	v1.print();
	if (!equals(v1.vec, expected1)) return 1;

	Tensor<double> v2 = u.accumulate(0, 2, std::plus<double>());
	expected1 = {15, 17, 19, 21, 23, 
                25, 27, 29, 31, 33, 
                35, 37, 39, 41, 43};
	v2.print();
	if (!equals(v2.vec, expected1)) return 1;

	u += 0.1;
	u.print();

	u = u - 0.1;
	u.print();

	Tensor<double> w({2,3,5});
	w.fill_sequence();
	w *= 0.01;
	w.print();

	Tensor<double> z = u+w;
	z.print();
	
	z /= 100;
	z.print();

	Tensor<double> u1 = u.max_dim(0);
	u1.print();

	Tensor<double> u2 = u.max_dim(1);
	u2.print();

	Tensor<double> u3 = u.max_dim(2);
	u3.print();

	Tensor<double> u11 = u.avg_dim(0);
	u11.print();

	Tensor<double> u12 = u.avg_dim(1);
	u12.print();

	Tensor<double> u21 = u.accumulate(0, 1, std::plus<double>(), {3,2,1});
	u21.print();

	Tensor<double> u31({2,3,5});
	u31.fill_sequence();
	u31.transform(0, std::plus<double>(), {4,3,2,1,0});
	u31.print();

	Tensor<double> u32({2,3,5});
	u32.fill_sequence();
	u32.transform(1, std::plus<double>(), {10,5,0});
	u32.print();

	return 0;

}

