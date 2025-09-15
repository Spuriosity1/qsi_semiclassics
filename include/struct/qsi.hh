/* Class qsi
 * Implements a cuboidal block of QSI classical dynamics
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

#ifndef qsi_hh
#define qsi_hh

#include <vec3.hh>
#include <random>
#include <set>

class spin;
class plaq;
class tetra;
class ptetra;
class fourier;
class static_correlator;
class dynamic_correlator;



 struct mc_sampler {
    const static unsigned sequential = 0u;
    const static unsigned random = 1u;
    mc_sampler() = delete;
    mc_sampler(unsigned mode, unsigned num_spins) : 
        plaquettes(num_spins), spins(num_spins)
    { 
        this->mode = mode;
        // num_plaquettes is always equal to the number of spins
        for (size_t j=0; j<num_spins; j++){
            plaquettes[j] = j;
            spins[j] = j;
        }
    }
    
    
    bool operator==(const mc_sampler& b) { return mode == b.mode; }
/*
    mc_sampler operator=(const mc_sampler& other) { 
        this->mode = other.mode;
        return *this;
    }
*/
    void ban_plaquette(size_t plaquette_idx);
    void ban_spin(size_t spin_idx);
        
    
    unsigned mode;

    const std::vector<size_t>& get_plaquettes() const { return plaquettes; };
    const std::vector<size_t>& get_spins() const { return spins; };
private:
    std::vector<size_t> plaquettes;
    std::vector<size_t> spins;
};

enum class mc_angle_mode{
    none, finite_t_von_mises, finite_t_gradient, zero_t
};

enum class mc_sz_mode {
    none, finite_t, zero_t
};


class qsi {
    // Number of unit cells along different directions
    const int na, nb, nc;

    // Number of spins/plaquettes
    const size_t N;

    // Spins
    spin** m_spin;

    // Plaquettes
    plaq** m_plaq;

    // Tetrahedra
    tetra** m_tetra;
    ptetra** m_ptetra;

    // State of spins
    vec3* m_state;

    // Random number engine
    std::mt19937 random;
    // Uniform sampling on [0,1[
    std::uniform_real_distribution<double> uniform;
    // Standard normal sampling
    std::normal_distribution<double> normal;
    // Uniform sampling of integers 0..(N-1)
    std::uniform_int_distribution<size_t> site;

    // Correlator management
    fourier& m_fourier;
    static_correlator* m_static;
    dynamic_correlator* m_dynamic;

    //---------- ROUTINES TO ASSEMBLE STRUCTURE ------------------------------
    // NB basic unit of distance is a_fcc/8
    
    /* Returns spin pointer corresponding to given position,
     * provided there is a spin there */
    spin* spin_at (const vec3_int& pos) const;
    
    /* Returns plaquette pointer corresponding to given position,
     * provided there is a plaquette centre there */
    plaq* plaq_at (const vec3_int& pos) const;

    /* Returns tetrahedron pointer at given position,
     * provided there is a tetrahedron centre there */
    tetra* tetra_at (const vec3_int& pos) const;

    /* Returns dual tetrahedron pointer at given position,
     * provided there is a dual tetrahedron centre there */
    ptetra* ptetra_at (const vec3_int& pos) const;


    //---------- PRIVATE IMPLEMENTATION OF MONTE CARLO ------------------------
    // Reshuffles XY angles on all tetrahedra (gauge invariant)
    void MC_gauge();

    // Does a MC step on all XY angles (von Mises distribution)
    void MC_angle(double T, const mc_sampler& sampling);

    // Propose a random XY rotation depending on temperature
    double MC_propose_angle (double T);

    // Does a MC step on all XY angles (small angle Metropolis by spin)
    // Return value: successful moves
    size_t MC_angle_small (double T, const mc_sampler& sampling);

    // Does a zero-temperature small angle Metropolis step on all XY angles
    // T is a handle on step sizes to ensure continuity
    // Return values: successful moves
    size_t MC_angle_freeze (double T, const mc_sampler& sampling);

    // Propose a random Ising tilting depending on temperature
    double MC_propose_Ising (double T);

    // Does a MC step on all Ising components (Metropolis by plaquette)
    // Return value: successful plaquette moves
    size_t MC_Ising (double T, const mc_sampler& sampling);
    
    // Does a zero-T MC step on all Ising components (Metropolis by plaquette)
    // T is a handle on step sizes to ensure continuity
    // Return value: successful plaquette moves
    size_t MC_Ising_freeze (double T, const mc_sampler& sampling);


    //---------- PUBLIC INTERFACE ---------------------------------------------
public:
    vec3_int shape(){
        return vec3_int(na, nb ,nc);
    }
    
    /**
     * @brief Construct a new qsi object
     * 
     * @param n Linear dimension of cube
     * @param f Fourier transformer object - this gets used and reused (not sure why it has to be an argument???)
     * @param s static_correlator pointer - can be null if no correlators are stored
     * @param d dynamic_correlator pointer - cann be null if no correlators are stored
     * @param seed Used to seed the RNG
     */
    qsi(int n, fourier* f, static_correlator* s = NULL,
        dynamic_correlator* d = NULL, size_t seed = 0);
    /**
     * @brief Construct a new qsi object
     * 
     * @param na # unit cells X dimension of cube
     * @param nb # unit cells Y dimension of cube
     * @param nc # unit cells Z dimension of cube
     * @param f Fourier transformer object - this gets used and reused (not sure why it has to be an argument???)
     * @param s static_correlator pointer - can be null if no correlators are stored
     * @param d dynamic_correlator pointer - cann be null if no correlators are stored
     * @param seed Used to seed the RNG
     */
    qsi(int na, int nb, int nc, fourier* f, static_correlator* s = NULL,
        dynamic_correlator* d = NULL, size_t seed = 0);

    // Destructor: clear up dynamically constructed objects
    ~qsi();

    // Number of spins/plaquettes
    size_t n_spin() const {return N;}
    // Number of tetrahedra
    size_t n_tetra() const {return N/2;}

    // ---------- Constant getters ----------

    const spin* spin_no(size_t n) const {return m_spin[n];}
    const plaq* plaq_no(size_t n) const {return m_plaq[n];}
    const  tetra*  tetra_no(size_t n) const {return  m_tetra[n];}
    const ptetra* ptetra_no(size_t n) const {return m_ptetra[n];}

    const spin* c_spin_at(const vec3_int& v) const {return spin_at(v);}
    const plaq* c_plaq_at(const vec3_int& v) const {return plaq_at(v);}
    const  tetra*  c_tetra_at(const vec3_int& v) const {return  tetra_at(v);}
    const ptetra* c_ptetra_at(const vec3_int& v) const {return ptetra_at(v);}


    /**
     * @brief Returns coordinats of the closest dual diamond lattice 
     * 
     * @param pos Position to hunt around
     * @param dmnd_sl On which diamond sublattice should we check?
     * @param tiebreak If equidistant from more than one site, choose the site that is closest to pos + tiebreak
     * @return vec3_int 
     */
    vec3_int nearest_ptetra(const vec3_int& pos, int dmnd_sl, const vec3_int& tiebreak=vec3_int(1,1,1)) const;
    vec3_int nearest_tetra(const vec3_int& pos) const;
    vec3_int nearest_spin(const vec3_int& pos) const;
    vec3_int nearest_plaq(const vec3_int& pos) const;


    size_t spin_idx_at  (vec3_int pos) const;
    size_t plaq_idx_at  (vec3_int pos) const;
    size_t tetra_idx_at (vec3_int pos) const;
    size_t ptetra_idx_at(vec3_int pos) const;

    // -------------- Setters -----------------
    void set_uniform_g(std::complex<double> g);
    void set_uniform_RK(double RK);

    // ---------- Classical dynamics ----------

    // Returns the full state of the QSI object as a non-const pointer
    vec3* state() {return m_state;}    
    
    /* Time derivative compatible with GSL ODE solvers */
    int evol(double t, const double* y, double* dydt);

    // ---------- Monte Carlo ----------

    static const unsigned
        MC_NONE = 0u, // No update
        MC_LARGE = 1u, // Large angle updates (von Mises)
        MC_SMALL = 2u, // Small angle/S^z updates
        MC_FREEZE = 3u, // Freezing
        
        MC_SEQ = 0u, // Traverse spins/plaquettes in order
        MC_RANDOM = 1u; // Pick spins/plaquettes at random

    /* Full Monte Carlo step
     * "sweep" rounds of XY angle and Ising MC; 1 round of  gauge shuffling
     * Type of MC algorithm set by argument method_Sz, method_angle
     * Sampling of spins/plaquettes set by argument sampling*/



    /**
     * @brief Runs a full Monte Carlo sweep
     * 
     * Performs 1 round of gauge shuffling, then n Ising MC steps
     * 
     * @param T Temperature
     * @param sweep Number of Ising MC steps
     * @param method_Sz Bit flag for Ising flip strategy. may be
     *                  - MC_NONE for no flips
     *                  - MC_SMALL for Metropolis by plaquette
     *                  - MC_FREEZE to flip only if energy is strictly smaller
     * @param method_angle Bit flag for angular relaxation strategy.
     *                  - MC_NONE for no relaxation
     *                  - MC_LARGE for von Mises sampled gradient descent
     *                  - MC_SMALL for uniform MC gradient descent
     *                  - MC_FREEZE for MC at zero temperature
     * @param sampling  Bit flag for order to do MC steps in
     *                  - MC_SEQ to go through all sites in order
     *                  - MC_RANDOM to visit random sites
     */
    void MC(double T, unsigned sweep, const mc_sampler& sampling,
            mc_sz_mode method_Sz = mc_sz_mode::finite_t,
            mc_angle_mode method_angle = mc_angle_mode::finite_t_von_mises
            );

    // ---------- Handling charges ----------
    
    // Introduce a positive and a negative spinon at given positions
    void add_spinon_pair(const vec3_int& pos, const vec3_int& neg, double q);
    
    /**
     * @brief Adds a vison pair between ptetras 'pos' and 'neg' using a dual string operator.
     * 
     * @param pos the location of the +ve vison
     * @param neg the location of the -ve vison
     * @param charge vison charge, living in [-2, 2] (mod 4)
     */
    void add_vison_pair(const vec3_int& pos, const vec3_int& neg, double charge=1, bool global=false);
    
    /**
     * @brief Peroforms a global singular gauge transform correesponding to a Dirac monopole solution between nearest neighbours pos and neg.
     * 
     * @param pos The location of the positive dipole.
     * @param neg The location of hte negative dipole.
     */
    void add_vison_dipole(const vec3_int& pos, const vec3_int& neg){
        add_vison_dipole(this->ptetra_at(pos), this->ptetra_at(neg));
    }
    void add_vison_dipole(ptetra* pos_site, ptetra* neg_site);
    

    // Vison charge at given position
    int vison(const vec3_int& pos) const;

    // geometric convenience functions

    // ptetra::sublattice ptetra_sl(const vec3_int& pos) const;

    
    // Total number of visons
    size_t num_visons() const;

    // Number of visons and vison pairs
    void vison_pair(size_t& vison, size_t& pair) const;

    // Total energy
    double energy() const;

    // Total instantaneous heat current
    vec3 heat_current() const;

    // Measures the Q_t inspired order parameter, a measure of vison asymmetry
    double vison_aiao_order() const;

    // Measures the Q_t inspired order parameter, a measure of vison asymmetry
    double spin_aiao_order() const;

    // Measure the chart-independent order parameter
    std::complex<double> avg_cisb() const;

    // void tetra_spanning_tree(std::vector<const spin*>& tree) const;
    
    // ---------- Saving correlators ----------

    // Stores static spin correlators in static_correlator object
    void save_static();
    // Stores time slice in dynamic_correlator object
    void save_dynamic(size_t t);
    // Stores either or both types of correlators
    void save(bool stat, bool dyn, size_t t = 0);

    // Saves the B values associated with all of the plaquettes
    void save_B( const std::string& bfile );

    
    // Saves the raw vector potential associated with all pyro sites
    void save_A( const std::string& afile );
    // Saves a file indicating vison occupancy of all tetrahedra
    void save_visons( const std::string& vfile );

    // ------------- STATE IO ------------
    // Save and load Heisenberg components of spins (with respect to local axes)
    void save_spins( const std::string& spinfile );
    void load_spins( const std::string& spinfile, bool strict=true );

    // coupling save/load
    void save_couplings( const std::string& couplingfile);
    void load_couplings( const std::string& couplingfile, bool require_specify_all=true );
};

/* Time derivative as an external function, for compatibility with GSL */
inline int qsi_evol(double t, const double* y, double* dydt, void* q) {
    return ((qsi*)q)->evol(t,y,dydt);
}

#endif
