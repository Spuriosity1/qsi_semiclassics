#pragma once 
#include <qsi.hh>
#include <random>

class defect_strategy {
    public:
        enum class choice_t {
            single_impurity,
            random,
            fixed_number,
            avoid_collision
        };

        defect_strategy(const qsi& sim, 
                const std::string& choice,
                double _defect_probability,
                std::random_device& rd,
                uint8_t subtlat_mask = 0x7 
                );
        std::vector<size_t> choose_spin_indices();
    private:
        choice_t _choice;
        const qsi& qsi_sim;
        double defect_probability;
        uint8_t sublat_mask;
        std::mt19937 rng;

        // stores spins that have already been chosen
        std::vector<unsigned int> selected_spin_idx;
};

