//
// Created by krab1k on 6.11.18.
//

#include<iostream>

#include "Method.h"
#include "Parameters.h"
#include "utility/Utility.h"


void Method::set_parameters(const Parameters *parameters) {
    if (common_parameters_.size() + atom_parameters_.size() + bond_parameters_.size() == 0 and parameters == nullptr) {
        return;
    }
    if (parameters->common() != nullptr and parameters->common()->names() != common_parameters_) {
        std::cerr << "Invalid common parameters provided" << std::endl;
        std::cerr << "Expected: " << common_parameters_ << std::endl;
        std::cerr << "Got: " << parameters->common()->names() << std::endl;
        exit(EXIT_FAILURE);
    }

    if (parameters->atom() != nullptr and parameters->atom()->names() != atom_parameters_) {
        std::cerr << "Invalid atom parameters provided" << std::endl;
        std::cerr << "Expected: " << atom_parameters_ << std::endl;
        std::cerr << "Got: " << parameters->atom()->names() << std::endl;
        exit(EXIT_FAILURE);
    }

    if (parameters->bond() != nullptr and parameters->bond()->names() != bond_parameters_) {
        std::cerr << "Invalid bond parameters provided" << std::endl;
        std::cerr << "Expected: " << bond_parameters_ << std::endl;
        std::cerr << "Got: " << parameters->bond()->names() << std::endl;
        exit(EXIT_FAILURE);
    }
    parameters_ = parameters;
}