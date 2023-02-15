#include <iostream>
#include <fish.h>
using namespace std;

int main(){

	string params_file = "params/cod_params.ini";
	Fish f(params_file);
	f.init(1.93e3, 5.61);

	cout << "Length = " << f.length << endl;
	
}
