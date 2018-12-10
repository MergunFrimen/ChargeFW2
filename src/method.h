#include <utility>

#include <utility>

//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <vector>

#include "structures/molecule.h"
#include "parameters.h"

class Method {
protected:
    const std::string name_{};
    const std::vector<std::string> common_parameters_{};
    const std::vector<std::string> atom_parameters_{};
    const std::vector<std::string> bond_parameters_{};

    const Parameters *parameters_{nullptr};

public:
    Method(std::string name, std::vector<std::string> common, std::vector<std::string> atom, std::vector<std::string> bond) :
            name_{std::move(name)},
            common_parameters_{std::move(common)},
            atom_parameters_{std::move(atom)},
            bond_parameters_{std::move(bond)} {}

    void set_parameters(const Parameters *parameters);

    bool has_parameters() {
        return (common_parameters_.size() + atom_parameters_.size() + bond_parameters_.size()) != 0;
    }

    virtual std::vector<double> calculate_charges(const Molecule &molecule) = 0;

    std::string name() { return name_; }
};
