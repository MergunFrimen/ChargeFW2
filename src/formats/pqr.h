//
// Created by krab1k on 31.1.19.
//

#pragma once

#include <string>

#include "writer.h"
#include "../structures/molecule_set.h"
#include "../charges.h"


class PQR : public Writer {
public:
    void save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) override;
    void append_charges_to_file(const MoleculeSet &ms, const Charges &charges, const std::string &filename);
};
