#include "defect_strategy.hh"
#include <algorithm>
#include <cstdint>
#include <random>
#include <sys/types.h>



defect_strategy::defect_strategy(const qsi& sim, 
        const std::string& choice,
        double _defect_probability,
        std::random_device& rd,
        uint8_t _sublat_mask
        ):
    qsi_sim(sim), defect_probability(_defect_probability),sublat_mask(_sublat_mask), rng(rd())
{
    if (choice == "single_impurity") {
        this->_choice = choice_t::single_impurity;
    } else if (choice == "random") {
        this->_choice = choice_t::random;
    } else if (choice == "fixed_number") {
        this->_choice = choice_t::fixed_number;
    } else /*
              if (choice == "avoid_collision") {
              this->_choice = choice_t::avoid_collision;
              } else */
    {
        throw std::runtime_error("invalid defect strategy");
    }
}

std::vector<size_t> defect_strategy::choose_spin_indices(){
    std::vector<size_t> sp;
    std::uniform_real_distribution<> unif(0, 1);
    std::uniform_int_distribution<>  randint(0, 1);
    
    
    int n_defects = qsi_sim.n_spin()*defect_probability;
    std::vector<int> possible_indices;

    for (size_t idx=0; idx<qsi_sim.n_spin(); idx++){
            possible_indices.push_back(idx);
    }
    std::shuffle(possible_indices.begin(), possible_indices.end(), rng); 
    // suboptimal but straightforward approach

    switch (this->_choice) {
        case choice_t::single_impurity: 
            sp.push_back(0);
            break;
        case choice_t::random:
            for (size_t idx=0; idx<qsi_sim.n_spin(); idx++){
                if (unif(rng) < defect_probability){    
                    selected_spin_idx.push_back(idx);
                    sp.push_back(idx);
                }
            }
            break;
        case choice_t::fixed_number:
            for (size_t j=0; j<n_defects; j++) {
                sp.push_back(possible_indices[j]);
            }
            break;
        case choice_t::avoid_collision:
        default:
            throw std::runtime_error("Not Implemented");
    }
    return sp;
}
