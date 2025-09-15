#include "common/seed.hh"
#include <iostream>

using namespace std;

int main(void){
	cout << "Max seeds: " << seed_raw_data::NUM_SEEDS << endl;
	cout << "These are..."<<endl;
	for (unsigned i=0; i<seed_raw_data::NUM_SEEDS; i++){
		cout << i<<" | "<< random_seed(i)<<"\n";
	}
	return 0;
}
