#include "vec3.hh"
#include <fstream>
#include <nlohmann/json.hpp>

using namespace nlohmann;

// Function to convert JSON data to plaquette coupling file
void convert_json_to_plaquette_couplings(const std::string& json_filename, const std::string& output_filename) {
    // Read JSON data
    std::ifstream json_file(json_filename);
    if (!json_file.is_open()) {
        std::cerr << "Failed to open JSON file: " << json_filename << std::endl;
        return;
    }
    
    json data;
    try {
        json_file >> data;
    } catch (const std::exception& e) {
        std::cerr << "Error parsing JSON: " << e.what() << std::endl;
        return;
    }
    
    // Extract plaquette positions from JSON data
    std::vector<vec3_int> plaquette_positions;
    if (data.contains("plaqs") && data["plaqs"].is_array()) {
        for (const auto& plaq : data["plaqs"]) {
            if (plaq.contains("pos") && plaq["pos"].is_array() && plaq["pos"].size() == 3) {
                vec3_int pos = {
                    plaq["pos"][0].get<int>(),
                    plaq["pos"][1].get<int>(),
                    plaq["pos"][2].get<int>()
                };
                plaquette_positions.push_back(pos);
            }
        }
    }
    
    std::cout << "Found " << plaquette_positions.size() << " plaquettes in JSON data" << std::endl;
    
    // Create output file
    std::ofstream output_file(output_filename);
    if (!output_file.is_open()) {
        std::cerr << "Failed to open output file: " << output_filename << std::endl;
        return;
    }
    
    // Write header
    output_file << "# Plaquette couplings generated from JSON data\n";
    output_file << "# Format: x y z g_constant RK_potential\n";
    
    // Write plaquette data
    for (const auto& pos : plaquette_positions) {
        output_file << pos[0] << " " << pos[1] << " " << pos[2] << " 1 0\n";
    }
    
    std::cout << "Successfully wrote plaquette couplings to " << output_filename << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cout << "Usage: " << argv[0] << " <input_json_file> <output_coupling_file>" << std::endl;
        return 1;
    }
    
    convert_json_to_plaquette_couplings(argv[1], argv[2]);

    return 0;
}

