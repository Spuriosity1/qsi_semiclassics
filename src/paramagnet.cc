#include <qsi.hh>
#include <spin.hh>
#include <fourier.hh>
#include <sstream>
#include <vec3.hh>
// Just generates a random spin configuration

// usage paramagnet N seed
int main(int argc, char** argv){

    unsigned n = atol(argv[1]);
    fourier f(n, false); // TODO use it for quadratic estimates?
    qsi simulate(n,&f,NULL,NULL,0);
    vec3* spins = simulate.state();

    std::normal_distribution<double> dist{0,1};
    std::mt19937 rand;
    for (int i=0; i < simulate.n_spin(); i++) {
        vec3 s(dist(rand),dist(rand),dist(rand));
        spins[i] = s.normalise();
    }

    std::stringstream output;
    output << "paramagnet%N=" << n <<"%seed=" <<argv[2]<< ".csv";

    simulate.save_spins(output.str());
    
}
    