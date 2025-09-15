/* Implementation of functions in class qsi
 *
 * Created on 13/06/2018
 * Copyright (C) 2018, 2019 Attila Szabó <as2372@cam.ac.uk>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the Creative Commons Attribution License (CC-BY),
 * version 4.0 of the License.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *
 * Please find a copy of the Creative Commons CC-BY License on
 * <https://creativecommons.org/licenses/by/4.0/>.
 */

#include "qsi.hh"
#include <vec3.hh>
#include <misc.hh>
#include <random.hh>
#include <seed.hh>
#include <direction.hh>
#include <util.hh>
#include "spin.hh"
#include "tetra.hh"
#include "plaq.hh"
#include "ptetra.hh"
#include "fourier.hh"
#include "static_correlator.hh"
#include "dynamic_correlator.hh"
#include "corrtypes.hh"

#include <algorithm>
#include <queue>
#include <cmath>
#include <random>
#include <complex>
#include <cstring>
#include <list>
#include <fstream>
#include <set>
#include <assert.h>
#include <array>
#include <vector>
#include "debug.hh"
#include <gsl/gsl_errno.h>


void mc_sampler::ban_plaquette(size_t plaquette_idx){
    plaquettes.erase(std::remove(plaquettes.begin(), plaquettes.end(), plaquette_idx), plaquettes.end());
}
void mc_sampler::ban_spin(size_t spin_idx){
    spins.erase(std::remove(spins.begin(), spins.end(), spin_idx), spins.end()); 
}


//-------------- ASSEMBLING STRUCTURE -----------------------------------------
// NEVER CHANGE INPUT TYPE TO CONST POINTER AS IT IS BEING CHANGED INSIDE

/* Returns the spin sitting at some point or throws an error if there is
 * nothing there.
 * Order in array m_spin: cubic[n]+fcc_Dy[i]+e[j] is at m_spin[16*n+4*i+j]
 * Cubic lattice points listed in row major order */
spin* qsi::spin_at (const vec3_int& pos) const {
    return m_spin[spin_idx_at(pos)];
}

size_t qsi::spin_idx_at (vec3_int pos) const {
    using namespace direction;
    int sl = -1;
    for (int i = 0; i < 4; i++) 
        if (alldiv(pos - pyro[i], 4)) {
            sl = i;
            break;
        }
    if (sl<0) throw std::runtime_error("Invalid position of spin");
    pos -= pyro[sl];

    int inc = -1;
    for (int i = 0; i < 4; i++) 
        if (alldiv(pos-fcc_Dy[i], 8)) {
            inc = i;
            break;
        }
    if (inc<0) throw std::runtime_error("Invalid position of fcc site");
    pos -= fcc_Dy[inc];
    
    int x = mod(pos[0]/8, na);
    int y = mod(pos[1]/8, nb);
    int z = mod(pos[2]/8, nc);
    int cell = (x*nb+y)*nc+z;

    return 16*cell+4*inc+sl;
}




/* Returns the plaquette sitting at some point or throws an error if there is
 * nothing there.
 * Order in array m_plaq: cubic[n]+fcc_Ti[i]+e[j] is at m_plaq[16*n+4*i+j] */
plaq* qsi::plaq_at (const vec3_int& pos) const {
    return m_plaq[plaq_idx_at(pos)];
}

size_t qsi::plaq_idx_at (vec3_int pos) const {
    using namespace direction;
    int sl = -1;
    for (int i = 0; i < 4; i++) 
        if (alldiv(pos - pyro[i], 4)) {
            sl = i;
            break;
        }
    if (sl<0) throw std::runtime_error("Invalid position of plaquette");
    pos -= pyro[sl];

    int inc = -1;
    for (int i = 0; i < 4; i++) 
        if (alldiv(pos-fcc_Ti[i], 8)) {
            inc = i;
            break;
        }
    if (inc<0) throw std::runtime_error("Invalid position of fcc site");
    pos -= fcc_Ti[inc];
    
    int x = mod(pos[0]/8, na);
    int y = mod(pos[1]/8, nb);
    int z = mod(pos[2]/8, nc);
    int cell = (x*nb+y)*nc+z;

    return 16*cell+4*inc+sl;
}


/* Returns the tetrahedron centre at some point or throws an error if there is
 * nothing there.
 * Order in array m_tetra: cubic[n]+fcc_Dy[i]+t[j] is at m_tetra[8*n+2*i+j] */
tetra* qsi::tetra_at (const vec3_int& pos) const {
    return m_tetra[tetra_idx_at(pos)];
}


size_t qsi::tetra_idx_at ( vec3_int pos) const {
    using namespace direction;
    int sl = -1;
    for (int i = 0; i < 2; i++) 
        if (alldiv(pos - diamond[i], 4)) {
            sl = i;
            break;
        }
    if (sl<0) throw std::runtime_error("Invalid position of spin");
    pos -= diamond[sl];

    int inc = -1;
    for (int i = 0; i < 4; i++) 
        if (alldiv(pos-fcc_Dy[i], 8)) {
            inc = i;
            break;
        }
    if (inc<0) throw std::runtime_error("Invalid position of fcc site");
    pos -= fcc_Dy[inc];
    
    int x = mod(pos[0]/8, na);
    int y = mod(pos[1]/8, nb);
    int z = mod(pos[2]/8, nc);
    int cell = (x*nb+y)*nc+z;

    return 8*cell+2*inc+sl;
}


/* Returns the dual tetrahedron centre at some point or throws an error if
 * there is nothing there.
 * Order in array m_ptetra: cubic[n]+fcc_Ti[i]+t[j] is at m_ptetra[8*n+2*i+j] */
ptetra* qsi::ptetra_at (const vec3_int& pos) const {
    return m_ptetra[ptetra_idx_at(pos)];
}

size_t qsi::ptetra_idx_at (vec3_int pos) const {
    using namespace direction;
    int sl = -1;
    for (int i = 0; i < 2; i++) 
        if (alldiv(pos - diamond[i], 4)) {
            sl = i;
            break;
        }
    if (sl<0) throw std::range_error("Invalid position of tetra");
    pos -= diamond[sl];

    int inc = -1;
    for (int i = 0; i < 4; i++) 
        if (alldiv(pos-fcc_Ti[i], 8)) {
            inc = i;
            break;
        }
    if (inc<0) throw std::range_error("Invalid position of fcc site");
    pos -= fcc_Ti[inc];
    
    int x = mod(pos[0]/8, na);
    int y = mod(pos[1]/8, nb);
    int z = mod(pos[2]/8, nc);
    int cell = (x*nb+y)*nc+z;

    return 8*cell+2*inc+sl;
}

// quick and dirty, not elegant at all
// also broken???
vec3_int qsi::nearest_ptetra(const vec3_int& pos, int dmnd_sl, const vec3_int& tiebreak) const {
    using namespace direction;
    vec3_int cell_estimate = pos / 8;
    std::array<vec3_int,4> best;
    for (auto& b : best){
        b = cell_estimate + fcc_Ti[0];
    } 

 
    int best_d2 = std::numeric_limits<int>::max();

    int num_best = 0;
    const vec3_int L(this->na*8, this->nb*8, this->nc*8);
    
    for (int IX = -1; IX<=1; IX++){
        for (int IY = -1; IY<=1; IY++){
            for (int IZ = -1; IZ<=1; IZ++){
                for (int i = 0; i < 4; i++) {
                    vec3_int X = (cell_estimate + vec3_int(IX,IY,IZ))*8 + fcc_Ti[i] + diamond[dmnd_sl];
                
                    if ( t3_metric2(X, pos, L) < best_d2) {
                        best_d2 = (X-pos).len2();
                        best[0] = X;
                        num_best=1;
                    } else if (t3_metric2(X, pos, L) == best_d2){
                        if (num_best < 6){
                            best[num_best] = X;
                        }
                        num_best++;
                    }
                }
            }
        }
    }

    if (num_best > 6){
        throw std::runtime_error("Site is equidistant from more than 6 dual lattice sites- that should be impossible!");
    }

    if (best_d2 == std::numeric_limits<int>::max()) {
        throw std::runtime_error("Failed to find closest ptetra.");
    }
    
    // TIEBREAKER ROUND
    best_d2 = std::numeric_limits<int>::max();
    vec3_int best_one;
    for (int i=0; i<num_best; i++){
        int d2 = (best[i] - pos - tiebreak).len2();
        if (d2 < best_d2){
            best_one = best[i];
            best_d2 =d2;
        }
    }
    if (best_d2 == std::numeric_limits<int>::max()) {
        throw std::runtime_error("Failed to tiebreak closest ptetra.");
    }

    assert( this->ptetra_at(best_one)->pos() == best_one );

    return best_one;
}
    
   /* 
vec3_int qsi::nearest_tetra(const vec3_int& pos) const {
    throw std::runtime_error("Not Implemented");
}
vec3_int qsi::nearest_spin(const vec3_int& pos) const {
    throw std::runtime_error("Not Implemented");
}
vec3_int qsi::nearest_plaq(const vec3_int& pos) const {
    throw std::runtime_error("Not Implemented");
}
*/



//-------------- CONSTRUCTOR, DESTRUCTOR --------------------------------------

qsi::qsi (int n, fourier* f, static_correlator* s, dynamic_correlator* d,
          size_t seed):
    qsi(n,n,n,f,s,d,seed)
{}

qsi::qsi (int na, int nb, int nc, fourier* f, static_correlator* s,
          dynamic_correlator* d, size_t seed):
    na(na),
    nb(nb),
    nc(nc),
    N(16ul*na*nb*nc),
    random(random_seed(seed)),
    uniform(0.0,1.0),
    normal(0.0,1.0),
    site(0,N-1),
    m_fourier(*f),
    m_static(s),
    m_dynamic(d)
{
    // Complain if no correlators would be stored
    //if ((!s) && (!d))
    //throw ("Specify either static or dynamic correlator objects");
    
    // Arrays that will keep pointers to dynamically constructed objects
    m_spin = new spin*[N];
    m_plaq = new plaq*[N];
    m_tetra = new tetra*[N/2];
    m_ptetra = new ptetra*[N/2];

    // State vector of vectors
    m_state = new vec3[N];

    // CONSTRUCT OBJECTS
    // Running indices
    size_t
        i_spin = 0,
        i_plaq = 0,
        i_tetra = 0,
        i_ptetra = 0;
    for (int x = 0; x < na; ++x) {
        for (int y = 0; y < nb; ++y) {
            for (int z = 0; z < nc; ++z) {
                using namespace direction;
                vec3_int cubic(8*x,8*y,8*z);
                
                // Spins and tetrahedra
                for (int fcc = 0; fcc < 4; ++fcc) {
                    for (int ssl = 0; ssl < 4; ++ssl) {
                        m_spin[i_spin] = new spin(cubic+fcc_Dy[fcc]+pyro[ssl],
                                                  (spin::sublattice)ssl,
                                                  i_spin, m_state);
                        ++i_spin;
                    }
                    m_tetra[i_tetra] = new tetra(cubic+fcc_Dy[fcc]+diamond[0],
                                                 tetra::sublattice::A);
                    ++i_tetra;
                    m_tetra[i_tetra] = new tetra(cubic+fcc_Dy[fcc]+diamond[1],
                                                 tetra::sublattice::B);
                    ++i_tetra;
                }

                // Plaquettes and dual tetrahedra
                for (int fcc = 0; fcc < 4; ++fcc) {
                    for (int psl = 0; psl < 4; ++psl) {
                        m_plaq[i_plaq] = new plaq(cubic+fcc_Ti[fcc]+pyro[psl],
                                                  (plaq::sublattice)psl);
                        ++i_plaq;
                    }
                    m_ptetra[i_ptetra]=new ptetra(cubic+fcc_Ti[fcc]+diamond[0],
                                                  ptetra::sublattice::A,
                                                  i_ptetra);
                    ++i_ptetra;
                    m_ptetra[i_ptetra]=new ptetra(cubic+fcc_Ti[fcc]+diamond[1],
                                                  ptetra::sublattice::B,
                                                  i_ptetra);
                    ++i_ptetra;
                }
            }
        }
    }

    // REGISTER OBJECTS WITH ONE ANOTHER
    // spins <-> tetrahedra
    for (size_t i = 0; i < N/2; ++i) {
        tetra* t = m_tetra[i];
        for (int j = 0; j < 4; ++j) {
            /* t[i]->sublat() evaluates to +/- 1,
               multiplying SL A based e[j] gives distance from actual centre */
            vec3_int r = t->pos() + static_cast<int>(t->sublat_enum()) * direction::pyro[j];
            spin* s = spin_at(r);
            s->register_plaquette(t);
            t->reg(s);
        }
    }

    // plaquettes <-> dual tetrahedra
    for (size_t i = 0; i < N/2; ++i) {
        ptetra* t = m_ptetra[i];
        for (int j = 0; j < 4; ++j) {
            /* t[i]->sublat() evaluates to +/- 1,
               multiplying SL A based e[j] gives distance from actual centre */
            vec3_int r = t->pos() + static_cast<int>(t->sublat_enum()) * direction::pyro[j];
            plaq* p = plaq_at(r);
            p->reg(t);
            t->reg(p);
        }
    }

    // spins <-> plaquettes
    for (size_t i = 0; i < N; ++i) {
        plaq* p = m_plaq[i];
        for (int j = 0; j < 6; ++j) {
            vec3_int r = p->pos() + direction::plaqt[p->sublat()][j];
            spin* s = spin_at(r);
            s->register_plaquette(p,j);
            p->reg(s,j);
        }
    }
}


qsi::~qsi() {
    delete[] m_state;
    
    for (size_t i = 0; i < N; ++i) {
        delete m_spin[i];
        delete m_plaq[i];
    }

    for (size_t i = 0; i < N/2; ++i) {
        delete m_tetra[i];
        delete m_ptetra[i];
    }

    delete m_spin;
    delete m_plaq;
    delete m_tetra;
    delete m_ptetra;
}


//-------------- CLASSICAL DYNAMICS ------------------------------------
// RHS of classical equation of motion
int qsi::evol(double t, const double* y, double* dydt) {
    const vec3* yv = (const vec3*) y;
    vec3* dv = (vec3*) dydt;

    for (size_t i = 0; i < N; ++i)
	dv[i] = m_spin[i] -> diff(yv);
    
    return GSL_SUCCESS;
}


//-------------- MONTE CARLO -------------------------------------------
// Gauge invariant XY angle shuffling
void qsi::MC_gauge() {
    for (size_t i = 0; i < N/2; ++i) {
        m_tetra[i]->rotate( 2*M_PI*uniform(random) );
    }
}

// XY angle update for all spins using von-Mises distribution
void qsi::MC_angle(double T, const mc_sampler& sampling) {
    const std::vector<size_t>& allowed_spins = sampling.get_spins();
    switch (sampling.mode) {
    case mc_sampler::sequential:
        // for (size_t i = 0; i < N; ++i){
        for (auto& i : allowed_spins){
            m_spin[i]->MC_angle(T, random);
        }
        break;
    case mc_sampler::random:
    // NB this could be sped up by drawing from 0... len(spins) and using allowed_spins[] as a lookup table
    // but this is almost certainly overkill
        
        for (unsigned j=0; j< allowed_spins.size(); j++){
            size_t i = site(random);
            if (std::binary_search(allowed_spins.begin(), allowed_spins.end(), i)){
                m_spin[i]->MC_angle(T, random);
            }
            // reject immediately if i is banned
        }
        break;
    default:
        throw std::runtime_error("Invalid randomisation method");
    }
}

// Propose a random Ising tilting: Gaussian with spread dependent on T
double qsi::MC_propose_angle(double T) {
    double stdev = 0.5*sqrt(T);
    return stdev * normal(random);
}

// Monte Carlo descent algorithm
size_t qsi::MC_angle_small(double T, const mc_sampler& sampling) {
    const std::vector<size_t>& allowed_spins = sampling.get_spins();
    size_t success = 0;
    spin* s;
    size_t i;
    for (auto& j : allowed_spins) {
        switch(sampling.mode) {
        case mc_sampler::sequential:
            s = m_spin[j];
            break;
        case mc_sampler::random:
            
            i = site(random);
            if ( !std::binary_search(allowed_spins.begin(), allowed_spins.end(), i)){
                continue;
            }
            s = m_spin[i];
            break;
        default:
            throw std::runtime_error("Invalid randomisation method");
        }
        double rot = MC_propose_angle(T);
        // The distribution sampled is ~exp(beta|h||S| cos(arg(S)-arg(h)))
        double dE = -std::real( std::conj(s->field_cplx()) * s->xy() *
                                (std::polar(1.0,rot) - 1.0) );
        if ((dE < 0) || (exp(-dE/T) > uniform(random)) ) {
            s->xy() *= std::polar(1.0,rot);
            ++success;
        }
        // NB rotational moves leave Ising intact, so ignore them
    }
    return success;
}

// Monte Carlo descent with effective zero temperature
size_t qsi::MC_angle_freeze(double T, const mc_sampler& sampling) {
    const std::vector<size_t>& allowed_spins = sampling.get_spins();
    size_t success = 0;
    for (auto& j : allowed_spins) {
        spin* s;
        size_t i;
        switch(sampling.mode) {
        case mc_sampler::sequential:
            s = m_spin[j];
            break;
        case mc_sampler::random:
            i = site(random);
                if ( !std::binary_search(allowed_spins.begin(), allowed_spins.end(), i)){
                    continue;
                }
            s = m_spin[i];
            break;
        default:
            throw std::runtime_error("Invalid randomisation method");
        }
        double rot = MC_propose_angle(T);
        // The distribution sampled is ~exp(beta|h||S| cos(arg(S)-arg(h)))
        double dE = -std::real( std::conj(s->field_cplx()) * s->xy() *
                                (std::polar(1.0,rot) - 1.0) );
        if ( dE < 0 ) {
            s->xy() *= std::polar(1.0,rot);
            ++success;
        }
    }
    return success;
}
    
// Propose a random Ising tilting: Gaussian with spread dependent on T
double qsi::MC_propose_Ising(double T) {
    double stdev = 0.25*sqrt(T);
    // The Ising component can't be larger than 1 even if T is really large...
    if (stdev > 0.5) 
        stdev = 0.5;
    return stdev * normal(random);
}

size_t qsi::MC_Ising(double T, const mc_sampler& sampling) {
    // Auxilliary array to store new configuration in each step
    vec3* trial = new vec3[N];
    memcpy(trial,m_state,N*sizeof(vec3));

    // Count successful tilts
    size_t n_tilt = 0;

    const std::vector<size_t>& allowed_plaqs = sampling.get_plaquettes();

    plaq* pp;
    size_t i;
    // Try to tilt each plaquette
    for (auto& j : allowed_plaqs) {
        switch(sampling.mode) {
        case mc_sampler::sequential:
            pp = m_plaq[j];
            break;
        case mc_sampler::random:
            i = site(random);
            if ( !std::binary_search(allowed_plaqs.begin(), allowed_plaqs.end(), i)){
                continue;
            }
            pp = m_plaq[i];
            break;
        default:
            throw std::runtime_error("Invalid randomisation method");
        }
        plaq& p = *pp;
        // Propose a tilting amount
        double tilt = MC_propose_Ising(T);
            double dE;

        for (int j = 0; j < 6; ++j) {
            // Perform tilting
            double z = (p[j].ising(trial) += tilt);
            // Reject immediately if z-component goes over 1
            if (fabs(z) > 1.0)
            goto reject_MC_Ising;
            // Ising component fine, normalise xy component
            double xy = sqrt(1-z*z);
            p[j].xy(trial) = std::polar(xy, std::arg(p[j].xy()));
                // Next spin must be tilted with opposite sign
                tilt = -tilt;
        }

        // Compare energies
        dE = p.energy_MC(trial) - p.energy_MC();
        /* Metropolis decision
        * true = rejection, false = acceptance
        * At the end of this step, trial == state */
        if ((dE>0) && (exp(-dE/T) < uniform(random)) ) {
            reject_MC_Ising:
            for (int j = 0; j < 6; ++j)
            p[j].heis(trial) = p[j].heis();
            continue;
        } else {
            for (int j = 0; j < 6; ++j)
            p[j].heis() = p[j].heis(trial);
            ++n_tilt;
        }
    }

    // Destroy trial; state should contain the updated state by now
    delete[] trial;

    return n_tilt;
}

size_t qsi::MC_Ising_freeze(double T, const mc_sampler& sampling) {
    // Auxilliary array to store new configuration in each step
    vec3* trial = new vec3[N];
    memcpy(trial,m_state,N*sizeof(vec3));

    // Count successful tilts
    size_t n_tilt = 0;

    const std::vector<size_t>& allowed_plaqs = sampling.get_plaquettes();
    plaq* pp;
    size_t i;
    // Try to tilt each plaquette
    for ( auto& j : allowed_plaqs ) {
        switch(sampling.mode) {
        case mc_sampler::sequential:
            pp = m_plaq[j];
            break;
        case mc_sampler::random:
            i = site(random);
            if ( !std::binary_search(allowed_plaqs.begin(), allowed_plaqs.end(), i)){
                continue;
            }
            pp = m_plaq[i];
            break;
        default:
            throw std::runtime_error("Invalid randomisation method");
        }
	plaq& p = *pp;
	// Propose a tilting amount
	double tilt = MC_propose_Ising(T);
        double dE;

	for (int j = 0; j < 6; ++j) {
	    // Perform tilting
	    double z = (p[j].ising(trial) += tilt);
	    // Reject immediately if z-component goes over 1
	    if (fabs(z) > 1.0)
		goto reject_MC_Ising_freeze;
	    // Ising component fine, normalise xy component
	    double xy = sqrt(1-z*z);
	    p[j].xy(trial) = std::polar(xy, std::arg(p[j].xy()));
            // Next spin must be tilted with opposite sign
            tilt = -tilt;
	}

	// Compare energies
	dE = p.energy_MC(trial) - p.energy_MC();
	/* Metropolis decision (T = 0)
	 * true = rejection, false = acceptance
	 * At the end of this step, trial == state */
	if ( dE>0 ) {
	reject_MC_Ising_freeze:
	    for (int j = 0; j < 6; ++j)
		p[j].heis(trial) = p[j].heis();
	    continue;
	} else {
	    for (int j = 0; j < 6; ++j)
		p[j].heis() = p[j].heis(trial);
	    ++n_tilt;
	}
    }

    // Destroy trial; state should contain the updated state by now
    delete[] trial;

    return n_tilt;
}

inline double sqr(double x) {return x*x;}

void qsi::MC(double T, unsigned sweep, const mc_sampler& sampling,
            mc_sz_mode method_Sz, mc_angle_mode method_angle
            )
             {
    // Randomly shuffle phases (preserves loop fluxes)
    // MC_gauge();
    //double prev = energy(), E;
    for (unsigned i = 0; i < sweep; ++i) {
        size_t success = 0;
        // Angular update: Rotate spins about z axis according to local field
        switch (method_angle) {
        case mc_angle_mode::none: break;
        case mc_angle_mode::finite_t_von_mises:
            MC_angle(T, sampling);
            break;
        case mc_angle_mode::finite_t_gradient:
            success += MC_angle_small(T, sampling);
            break;
        case mc_angle_mode::zero_t:
            success += MC_angle_freeze(T, sampling);
            db_printerr("[MC Freeze] success %.4f %%",100.*(double)success/N);
            break;
        default:
            throw std::runtime_error("Invalid angular MC method");
        }
        // Ising update: Rotate Z axis into the XY plane in a ringflip pattern
        switch (method_Sz) {
        case mc_sz_mode::none: break;
        case mc_sz_mode::finite_t:
            success = MC_Ising(T, sampling);
            db_printerr("[MC Metropolis] success %.4f %%\n",100.*(double)success/N);
            break;
        case mc_sz_mode::zero_t:
            success = MC_Ising_freeze(T, sampling);
            db_printerr("[MC Ising Freeze] success %.4f%%\n",100.*(double)success/N);
            break;
        default:
            throw std::runtime_error("Invalid Ising MC method");
        }
        fflush(stderr);
        //E = energy();
        db_printerr("[MC SWEEP] completed %u / %u\n", i+1, sweep); fflush(stderr);
        //prev = E;
    }
}



//-------------- CHARGE MANAGEMENT --------------------------------------------


void qsi::add_spinon_pair(const vec3_int& pos, const vec3_int& neg, double q) {
    tetra* plus = tetra_at(pos);
    tetra* minus = tetra_at(neg);

    // Trivial but bug-prone case
    if (plus == minus)
        return;
    
    for (size_t i = 0; i < N/2; ++i)
        m_tetra[i]->prev = NULL;

    // breadth first search from minus to plus
    std::list<tetra*> bfs;
    bfs.push_back(minus);
    while(true) { // break this loop once plus is encountered
        tetra* cur = bfs.front();
        bfs.pop_front();
        
        for (int i = 0; i < 4; ++i) {
            tetra* other = cur -> neighbour(i) -> t_other(cur);
            if ((other != minus) && (other->prev == NULL)) {
                other->prev = cur->neighbour(i);
                bfs.push_back(other);
                if (other == plus)
                    // can only break from nested loop using goto
                    goto after_bfs_spinon;
            }
        }
    }
    
 after_bfs_spinon:
    // walk back from plus to minus and perform vison flips along the way
    tetra* cur = plus;
    while (cur != minus) {
        cur->prev->spinon_hop(q,cur);
        cur = cur->prev->t_other(cur);
    }
}


void qsi::add_vison_pair(const vec3_int& pos, const vec3_int& neg, double charge, bool global) {
    // BFS to find shortest path, then manually create a solenoid
    ptetra* plus;
    ptetra* minus;
    try {
        plus = ptetra_at(pos);
    } catch (const std::range_error& e) {
        std::stringstream ss;
        ss << "Position " << pos << " is not a ptetra, nearest are (A) "<<
            this->nearest_ptetra(pos, 1) << ", (B) "<<
            this->nearest_ptetra(pos, 1)<<"\n";
        throw std::range_error(ss.str());
    }
    try {
        minus = ptetra_at(neg);
    } catch (const std::range_error& e) {
        std::stringstream ss;
        ss << "Position " << neg << " is not a ptetra, nearest are (A) "<< 
            this->nearest_ptetra(neg, 0) << ", (B) "<<
            this->nearest_ptetra(neg, 1)<<"\n";
        throw std::range_error(ss.str());
    }

    // Trivial but bug-prone case
    if (plus == minus){
        return;
    }
    
    for (size_t i = 0; i < N/2; ++i)
        m_ptetra[i]->prev = NULL;

    // breadth first search from minus to plus
    std::list<ptetra*> bfs;
    bfs.push_back(minus);
    while(true) { // break this loop once plus is encountered
        ptetra* cur = bfs.front();
        bfs.pop_front();
        
        for (int i = 0; i < 4; ++i) {
            ptetra* other = cur -> neighbour(i) -> t_other(cur);
            if ((other != minus) && (other->prev == NULL)) {
                other->prev = cur->neighbour(i);
                bfs.push_back(other);
                if (other == plus)
                    // can only break from nested loop using goto
                    goto after_bfs_vison;
            }
        }
    }
    
    after_bfs_vison:
    // walk back from plus to minus and perform vison flips along the way
    ptetra* cur = plus;
    while (cur != minus) {
        if (global) {
            add_vison_dipole(cur, cur->prev->t_other(cur));
        } else {
            cur->prev->vison_hop(cur, charge);
        }
        cur = cur->prev->t_other(cur);
    }
}

void qsi::add_vison_dipole(ptetra* pos_site, ptetra* neg_site){
    // Adds in a Dirac dipole solution connected by a straight line.
    vec3 b = pos_site->pos() - neg_site->pos();
    vec3 x = 0.5*(pos_site->pos() + neg_site->pos());

    
    for (unsigned j=0; j<this->n_spin(); j++){
        spin* s = this->m_spin[j];
        double phi = Dirac_dipole(x, b, s->tA()->pos(),s->tB()->pos());
        s->xy() *= std::polar(1.0, phi);
    } 
}

// New approach: explicitly solve the lattice PDE using superpositions
// of the Dirac monopole solution.



int qsi::vison(const vec3_int& pos) const {
    return ptetra_at(pos) -> vison();
}

// ptetra::sublattice qsi::ptetra_sl(const vec3_int& pos) const {
//     return ptetra_at(pos)->sublat_enum();
// }

inline unsigned sqr(int x) {return x*x;}

// Total number of visons (both positive and negative) in the system
// should always be even (if it's not, God only knows what has happened)
size_t qsi::num_visons() const {
    size_t retval = 0;
    for (size_t i = 0; i < N/2; ++i)
        retval += sqr(m_ptetra[i] -> vison());
    // if (retval%2 == 1){
    //     throw std::logic_error("Magnetic monopoles have been 
    //     found in the vacuum, call Dirac");
    // }
    return retval;
}

/**
 * @brief Counts number of visons and number of vison dipoles
 * 
 * @param v number of visons
 * @param p number of vison dipoles
 */
void qsi::vison_pair(size_t& v, size_t& p) const {
    v = 0;
    for (size_t i = 0; i < N/2; ++i)
        v += sqr(m_ptetra[i] -> vison());
    p = 0;
    for (size_t i = 0; i < N; ++i)
        /* Since vison charge is virtually always 1,0, or -1,
         * there only is a dipole if the product of the two charges is -1
         * As vison charges were all updated above, can use memoised values */
        if (((m_plaq[i]->tA()->vison(true)) *
             (m_plaq[i]->tB()->vison(true))  ) == -1)
            ++p;
}

double qsi::energy() const {
    double retval = 0.0;
    for (size_t i = 0; i < N; ++i){
        // ring-exchange contribution
        retval += m_plaq[i]->ring_energy();
        // ising contribution (double counts, oh well)
        // retval += m_spin[i]->zfield_ising_3nn()*m_spin[i]->ising()/2;
    }
    
    return retval;
}


double qsi::spin_aiao_order() const {
    double retval = 0;
    for (size_t i = 0; i < N; ++i)
        retval += m_spin[i]->ising();
    return retval;
}

vec3 qsi::heat_current() const {
    vec3 retval;
    for (size_t i = 0; i < N; ++i)
        retval += m_plaq[i]->heat_current();
    return retval;
}

double qsi::vison_aiao_order() const {
    double order = 0;
    for (size_t i=0; i<n_tetra(); i++) {
        const ptetra* p = ptetra_no(i);
        order += p->vison()*static_cast<int>(p->sublat_enum());
    }
    order /= 1.*n_tetra();
    return order;
}



std::complex<double> qsi::avg_cisb() const {
    std::complex<double> order = 0;
    for (size_t j=0; j<this->n_spin(); j++) {
        order += plaq_no(j)->ring();
    }
    order /= this->n_spin();
    return order;
}


//-------------- SAVING CORRELATORS -------------------------------------------

void qsi::save(bool stat, bool dyn, size_t t) {
    // Don't attempt saving into NULL pointers
    stat = stat && m_static;
    dyn  = dyn && m_dynamic;

    if (!stat && !dyn)
        return;
    
    using namespace corrtypes;
    // Only generate correlators actually asked for
    unsigned what = 0;
    if (stat) what |= m_static -> what;
    if (dyn) {
        what |= m_dynamic -> what;
        if ((m_dynamic -> what) & NOISE) what |= SZZ;
    }

    // Sublattice index
    for (int i = 0; i < 4; ++i) {
        
        if (what & SZZ) { // S^z
            // Clean Fourier transform array
            m_fourier.clear_input();
            for (size_t s = i; s < N; s+=4) {
                /* Fill input of m_static with spins of sublattice i
                 * Move them to fcc lattice sites of "A" tetrahedra */
                m_fourier( m_spin[s]->fcc() / 4) = m_spin[s]->get_g_factor() *
                    m_spin[s]->ising();
            }

            m_fourier.transform();
            // If static correlators are collected, store results there
            if (stat)
                m_static -> set(i, SZZ);
            // If dynamic correlators are collected, store results there
            if (dyn)
                m_dynamic -> set(t, i, SZZ);
        }

        if (what & SXX) { // S^x
            // Clean Fourier transform array
            m_fourier.clear_input();
            for (size_t s = i; s < N; s+=4) {
                /* Fill input of m_static with spins of sublattice i
                 * Move them to fcc lattice sites of "A" tetrahedra */
                m_fourier( m_spin[s]->fcc() / 4 ) = m_spin[s]->get_g_factor() *
                    m_spin[s]->heis()[0];
            }
            m_fourier.transform();
            // If static correlators are collected, store results there
            if (stat)
                m_static -> set(i, SXX);
            // If dynamic correlators are collected, store results there
            if (dyn)
                m_dynamic -> set(t, i, SXX);
        }

        if (what & SYY) { // S^y
            // Clean Fourier transform array
            m_fourier.clear_input();
            for (size_t s = i; s < N; s+=4) {
                /* Fill input of m_static with spins of sublattice i
                 * Move them to fcc lattice sites of "A" tetrahedra */
                m_fourier( m_spin[s]->fcc() / 4 ) = m_spin[s]->get_g_factor() *
                    m_spin[s]->heis()[1];
            }
            m_fourier.transform();
            // If static correlators are collected, store results there
            if (stat)
                m_static -> set(i, SYY);
            // If dynamic correlators are collected, store results there
            if (dyn)
                m_dynamic -> set(t, i, SYY);
        }

        if (what & SPM) { // S^+S^-
            m_fourier.clear_input();
            for (size_t s = i; s < N; s+=4) {
                m_fourier( m_spin[s]->fcc() / 4 ) = m_spin[s]->get_g_factor() *
                    m_spin[s]->xy();
            }
            m_fourier.transform();
            if (stat)
                m_static -> set(i, SPM);
            if (dyn)
                m_dynamic -> set(t, i, SPM);
        }

        if (what & IMAG) { // b = Im(S^+S^-...)
            m_fourier.clear_input();
            for (size_t s = i; s < N; s+=4) {
                m_fourier( m_plaq[s]->fcc() / 4) = m_spin[s]->get_g_factor() *
                    m_plaq[s]->ring().imag();
            }
            m_fourier.transform();
            if (stat)
                m_static -> set(i, IMAG);
            if (dyn)
                m_dynamic -> set(t, i, IMAG);
        }

        if (what & ARG) { // B = arg(S^+S^-...)
            m_fourier.clear_input();
            for (size_t s = i; s < N; s+=4){
                m_fourier( m_plaq[s]->fcc() / 4) = m_plaq[s]->B();
            }
            m_fourier.transform();
            if (stat)
                m_static -> set(i, ARG);
            if (dyn)
                m_dynamic -> set(t, i, ARG);
        }

        if (what & SIN) { // sin B, essentially the same as IMAG
            m_fourier.clear_input();
            for (size_t s = i; s < N; s+=4){
                m_fourier( m_plaq[s]->fcc() / 4) = m_plaq[s]->ring().imag();
            }
            m_fourier.transform();
            if (stat)
                m_static -> set(i, SIN);
            if (dyn)
                m_dynamic -> set(t, i, SIN);
        }

        if (what & COS) { // cos B
            m_fourier.clear_input();
            for (size_t s = i; s < N; s+=4)
                m_fourier( m_plaq[s]->fcc() / 4) = m_plaq[s]->ring().real();
            m_fourier.transform();
            if (stat)
                m_static -> set(i, COS);
            if (dyn)
                m_dynamic -> set(t, i, COS);
        }
    }

    // If static correlators are stored, store them now
    if (stat)
        m_static -> add_correlator();
}

// Trivial instantiations of save()
void qsi::save_static() {
    save(true,false,0);
}

void qsi::save_dynamic(size_t t) {
    save(false,true,t);
}

// void initialise_random(){
//     throw "Not Implemented!";
// }

// -------------- Setters -----------------
void qsi::set_uniform_g(std::complex<double> g){
    for (size_t i=0; i<this->n_spin(); i++){
        this->m_plaq[i]->g_constant = g;
    }
}

// void qsi::set_g(const unsigned idx, std::complex<double> g){
//     this->m_plaq[idx]->g_constant = g;
// }

// void qsi::set_g(const vec3_int& x, std::complex<double> g){
//     set_g(x, g);
// }

void qsi::set_uniform_RK(double RK){
    for (size_t i=0; i<this->n_spin(); i++){
        this->m_plaq[i]->RK_potential = RK;
    }
}


// ---------- File IO -----------

void qsi::save_B( const std::string& bfile){
    FILE* B_out = fopen(bfile.c_str(), "w");
    // Print plaquette data in format X Y Z B Energy
    fprintf(B_out, "#   X   Y   Z       B       Energy\n");
    for (size_t i = 0; i < this->n_spin(); i++) {
        const plaq* p = this->plaq_no(i);
        fprintf(B_out, "%4d %4d %4d %.16e %.16e \n", p->pos()[0], p->pos()[1], p->pos()[2], p->B(), p->ring_energy()+1);
    }
    fclose(B_out);
}

/*
void qsi::tetra_spanning_tree(std::vector<const spin*>& tree) const {
    // BFS 
    assert( this->n_spin() > 0 );
    bool visited[ this->n_spin() ]; // set of vertices in the tree, tagged with 'false' if not there
    for (int j=0; j<this->n_spin(); j++){
        visited[j] = false;
    }
    visited[0] = 1;

    // BFS
    const tetra* curr;
    const spin* s;
    const tetra* t_neighbour;
    // curr->prev = nullptr;

    std::queue<const tetra*> to_visit;
    to_visit.push(this->m_tetra[0]);

    while(!to_visit.empty()){
        curr = to_visit.front();
        to_visit.pop();
        for (unsigned mu=0; mu<4; mu++){
            s = curr->neighbour(mu);
            t_neighbour = s->t_other(curr);
            unsigned t_neighbour_idx = this->tetra_idx_at(t_neighbour->pos()); //O(1) operation

            if (!visited[t_neighbour_idx]){
                to_visit.push(t_neighbour);
                visited[t_neighbour_idx] = true;
                tree.push_back(s);
                // t_neighbour->prev = s;
            }
        }
    }
}
*/

void qsi::save_A(const std::string& afile){
    FILE* spin_out = fopen(afile.c_str(), "w");
    
    fprintf(spin_out, "X   Y   Z pyro_sl A \n");
    for (size_t i = 0; i < n_spin(); i++) {
        const spin* s = spin_no(i);
        vec3_int pos = s->pos();
        double vi = arg(s->xy());
        fprintf(spin_out, "%4d %4d %4d %1d %.16e \n", 
            pos[0], pos[1], pos[2], s->sublat(), vi);
    }
    fclose(spin_out);
}



void qsi::save_visons( const std::string& vfile){
    FILE* vison_out = fopen(vfile.c_str(), "w");
    fprintf(vison_out, "X Y Z sl divB Energy\n");

    for (size_t i=0; i<this->n_tetra(); i++) {
        const ptetra* p = this->ptetra_no(i);
        fprintf(vison_out, "%4d %4d %4d %+1d %+d %.16e \n",
            p->pos()[0], p->pos()[1], p->pos()[2], 
            static_cast<int>(p->sublat_enum()), p->vison(), p->energy()+1);
    }
    fclose(vison_out);
}

void qsi::save_spins(const std::string& filename){
    FILE* spin_out = fopen(filename.c_str(), "w");
    
    fprintf(spin_out, "#   X   Y   Z   sX   sY   sZ   g_factor\n");
    for (size_t i = 0; i < n_spin(); i++) {
        const spin* s = spin_no(i);
        vec3_int pos = s->pos();
        vec3 vi = s->heis();
        fprintf(spin_out, "%4d %4d %4d %.16e %.16e %.16e %.16e \n", 
                pos[0], pos[1], pos[2], 
                vi[0], vi[1], vi[2], 
                s->get_g_factor());
    }
    fclose(spin_out);
}

void qsi::load_spins(const std::string& filename, bool strict){
    std::string line;
    std::ifstream myfile;
    myfile.open(filename);

    if(!myfile.is_open()) {
        perror("Could not open spinfile");
        exit(EXIT_FAILURE);
    }
    
    vec3_int pos;
    vec3 v;
    double g_factor;

    std::set<const spin*> uninitialised;
    if (strict){
        for (size_t j=0; j<this->n_spin(); j++){
            uninitialised.insert(this->spin_no(j));
        }
    }


    while(getline(myfile, line)) {
        if (line[0] =='#') continue;
        if (line.size() == 0) continue;
            
        std::istringstream iss(line);
        iss >> pos[0] >> pos[1] >> pos[2] >> v[0] >> v[1] >> v[2];
        
        spin* s = spin_at(pos);

        if (strict){
            size_t n = uninitialised.erase(s);
            if (n == 0){
                // attempting double-specification, the user probably doesn't want to do this!
                std::cerr<<"Spin Loading Failed: Attempting to re-initialise spin at [ "<<
                pos[0] << ", " << pos[1] << ", " << pos[2] << " ] " << std::endl;
                throw std::runtime_error("Bad Spin Intialiser");
            }
        }

        s->heis() = v;

        if ( !iss.eof() ) {
            iss >> g_factor;
            s->set_g_factor(g_factor);
        }
    }

    if (strict && !uninitialised.empty()) {
        std::cerr << "Spin Loading Failed: Spins at"<<std::endl;
        for (auto&p : uninitialised){
            const vec3_int where = p->pos();
            std::cerr << " [ "
                << where[0] << ", " << where[1] << ", " << where[2]
                << " ]" << std::endl;
        }
        std::cout << "Were not initialised. " << std::endl;
        throw std::runtime_error("Bad Spin Intialiser");
    }

}


void qsi::load_couplings( const std::string& couplingfile, bool require_specify_all ){
    std::ifstream ifs;
    ifs.open(couplingfile,std::ios::in);
    if (!ifs.is_open()){
        std::cerr<<"Failed to open coupling file"<<couplingfile<<std::endl;
        throw std::runtime_error("Failed to open couplingfile");
    }
    std::string line;
    vec3_int R;

    std::set<const plaq*> uninitialised;
    if (require_specify_all){
        for (size_t j=0; j<this->n_spin(); j++){
            auto res = uninitialised.insert(this->plaq_no(j));
            if (res.second == false){
                std::cerr << "No insert!" << this->plaq_no(j)->pos();
                throw std::runtime_error("plaq_no returned a duplicate pointer");
            }
        }
    }
    

    while ( std::getline(ifs, line) ){
        if ( line[0] == '#' ) continue;
        if (line.size() == 0) continue;

        std::istringstream ss(line);
        ss >> R[0] >> R[1] >> R[2];
        plaq* plaquette;
        try {
            plaquette = this->plaq_at(R);
        } catch (const std::runtime_error& e){
            std::cerr<<"Plaquette could not be found " << e.what() << std::endl;
            throw std::runtime_error("Plaquette could not be found");
        }
        std::cout << R[0] << " " << R[1] << " " << R[2] << " " << plaquette->pos() << "\n";

        if (require_specify_all){
            size_t n = uninitialised.erase(plaquette);
            if (n == 0){
                // attempting double-specification, the user probably doesn't want to do this!
                std::cerr<<"Plaquette Loading failed: Attempting to re-initialise plaquette at [ "<<
                R[0] << ", " << R[1] << ", " << R[2] << " ] " << std::endl;
                throw std::runtime_error("Bad Plaquette Intialiser");
            }
        }
        
        ss >> plaquette->g_constant >> plaquette->RK_potential;   
    }

    if (require_specify_all && !uninitialised.empty()) {
        std::cerr << "Plaquette Loading Failed. Plaqs at"<<std::endl;
        for (auto&p : uninitialised){
            const vec3_int where = p->pos();
            std::cerr << " [ "
                << where[0] << ", " << where[1] << ", " << where[2]
                << " ]" << std::endl;
        }
        std::cout << "Were not initialised. " << std::endl;
        throw std::runtime_error("Bad Plaquette Intialiser");
    }
}

void qsi::save_couplings( const std::string& couplingfile ){
    FILE* plaq_out = fopen(couplingfile.c_str(),"w");
    
    fprintf(plaq_out, "# PLAQUETTES\n");
    fprintf(plaq_out, "# x y z g RK_potential\n");
    for (size_t j=0; j<this->n_spin(); j++){
        const plaq* pl = this->plaq_no(j);
        vec3_int pos = pl->pos();
        // ofs << v[0] << " " << v[1] << " " <<v[2] << " ";
        // ofs <<  << " " << pl->RK_potential << "\n";
        fprintf(plaq_out, "%4d %4d %4d (%.16lf,%.16lf) %.16e\n", pos[0], pos[1], pos[2], 
            pl->g_constant.real(), pl->g_constant.imag(), pl->RK_potential);
    }
    fclose(plaq_out);
}
