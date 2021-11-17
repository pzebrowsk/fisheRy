#ifndef RFISH_JSS_H_
#define RFISH_JSS_H_

#include <iostream>
#include <cassert>
#include <vector>
#include <functional>
#include <algorithm>

#include <numeric>

// Multidimensional vector with
// dimensions [..., N2, N1, N0 ] 
// indexed as [..., i2, i1, i0 ]
// i0 is the lowest dimension (elements along i0 are stored consequtively in memory.
// 
// Example, an arry with dimensions [2,3,5] is stored as:
// loc: index
//	 0: 0 0 0
//	 1: 0 0 1 
//	 2: 0 0 2 
//	 3: 0 0 3 
//	 4: 0 0 4 
//	 5: 0 1 0 
//	 6: 0 1 1 
//	 7: 0 1 2 
//	 8: 0 1 3 
//	 9: 0 1 4 
//	10: 0 2 0 
//	11: 0 2 1 
//	12: 0 2 2 
//	13: 0 2 3 
//	14: 0 2 4 
//	15: 1 0 0 
//	16: 1 0 1 
//	17: 1 0 2 
//	18: 1 0 3 
//	19: 1 0 4 
//	20: 1 1 0 
//	21: 1 1 1 
//	22: 1 1 2 
//	23: 1 1 3 
//	24: 1 1 4 
//	25: 1 2 0 
//	26: 1 2 1 
//	27: 1 2 2 
//	28: 1 2 3 
//	29: 1 2 4  

class MultiVector{
	private:
	std::vector<int> dim;
	std::vector<int> offsets;
	int nelem;
	
	public:
	std::vector<double> vec;

	MultiVector(std::vector<int> _dim){
		dim = _dim;
		nelem = std::accumulate(dim.begin(), dim.end(), 1, std::multiplies<int>());
		vec.resize(nelem);

		std::cout << "nelem = " << nelem << "\n";

		int ndim = dim.size();
		offsets.resize(ndim,0);
		int p = 1;
		for (int i=ndim-1; i>=0; --i){
			offsets[i] = p;
			p *= dim[i];
		}
		
		std::cout << "offsets: "; for (auto o : offsets) std::cout << o << " "; std::cout << "\n";

	}

	double& operator() (std::vector<int> ix){
		int loc = 0;
		int ndim = dim.size();
		for (int i=ndim-1; i>=0; --i){
			loc += offsets[i]*ix[i];
		}
		return vec[loc];
	}

	std::vector<int> index(int loc){
		int ndim = dim.size();
		std::vector<int> id(ndim);
		for (int i=ndim-1; i>=0; --i){
			int ix = loc % dim[i];
			loc = (loc-ix)/dim[i];
			id[i]=ix;
		}
		return id;
	}

	// for debug only
	void fill_sequence(){
		for(int i=0; i<vec.size(); ++i) vec[i]=i;
	}

	template <class BinOp>
	Multivector accumulate(int dim, BinOp binary_op){
			
	}

};

using namespace std;

int main(){

	MultiVector u({2,3,5});
	u.fill_sequence();

	for (int i2 = 0; i2<2; ++i2){
		for (int i1=0; i1<3; ++i1){
			for (int i0=0; i0<5; ++i0){
				cout << u({i2,i1,i0}) << " ";
			}
			cout << "\n";
		}
		cout << "\n";
	}
	cout << "\n";



	for (int i=0; i<30; ++i){
		auto v = u.index(i);
		cout << i << ": ";
		for (auto x : v) cout << x << " ";
		cout << "\n";
	}

	

}

#endif
