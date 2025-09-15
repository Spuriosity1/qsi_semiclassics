// Executable for introducing a random spin configuration and iannealing it to find the gorund state.
// Usage: anneal /path/to/input.txt /path/to/output.csv seed
// Infile format:
// - Initial temperature
// - Annelaing steps



#include <basic_parser.hh>
#include <qsi.hh>
#include <fourier.hh>
#include <vec3.hh>
#include <plaq.hh>
#include <iostream>
#include <cstring>

int main(int argc, char** argv){
    if (argc<4){
        std::cerr << "Invalid number of arguments passed\n";
        std::cerr << "Usage: anneal infile outfile seed" <<std::endl;
        exit(EXIT_FAILURE);
    }


    basic_parser<int, unsigned, double> p;
    unsigned n, seed, sweep;
    double T, gphase;
    std::string prefix;
    p.declare("system_size", &n);
    p.declare("n_sweep", &sweep);
    p.declare("temperature", &T);
    p.declare("gphase", &gphase);
    p.declare_optional("prefix",&prefix,"untitled");

    p.from_file(argv[1]);
    seed = atoi(argv[3]);

    

    std::string output_path(argv[2]);

    fourier f(n, false); 
    qsi simulate(n,&f,NULL,NULL,seed);

    simulate.set_uniform_g( std::polar(1., gphase));

    const unsigned nspin = simulate.n_spin();
    const unsigned block_size = nspin*sizeof(vec3);
    vec3* out_data = new vec3[block_size*sweep];
    for (unsigned i=0;i<sweep;i++){
        simulate.MC(T, 1);
        vec3* spinblock = simulate.state();
        memcpy(out_data + i*block_size, spinblock, block_size);
    }

    
    
    std::string filestem(output_path);
    filestem.append("/");
    filestem.append(prefix);
    filestem.append("+");
    filestem.append(std::to_string(seed));

    std::string s(filestem);
    s.append("cooling.spins");

    FILE* fp = fopen(s.c_str(),"wb");

    fwrite(out_data, sizeof(vec3), nspin*sweep, fp);
    delete[] out_data;
    fclose(fp);

    // save the input params in the same spot
    s = filestem;
    s.append("_input.inp");
    p.into_file(s.c_str());
}