//
// Created by danny305 on May 3rd, 2021.
//

#include <iostream>
#include <string>
#include <vector>
#include <stdexcept>
#include <fstream>
#include <filesystem>

#include <gemmi/align.hpp>  // assign_label_seq_id
#include <gemmi/cif.hpp>
#include <gemmi/mmcif.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/to_cif.hpp>
#include <gemmi/modify.hpp>

// suppresses warnings from gemmi/to_mmcif.hpp
// unfortunately only works with clang
#define GEMMI_WRITE_IMPLEMENTATION
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wcpp"
#include <gemmi/to_mmcif.hpp>  // make_mmcif_block
#pragma GCC diagnostic pop
#undef GEMMI_WRITE_IMPLEMENTATION

#include "chargefw2.h"
#include "cif.h"
#include "../structures/molecule_set.h"
#include "../charges.h"
#include "../config.h"
#include "../utility/strings.h"

static std::string convert_bond_order_to_mmcif_value_order_string(int order) {
    // See link for bond order values for PDBx/mmCIF category _chem_comp_bond.value_order
    // https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Items/_chem_comp_bond.value_order.html#papwtenum
    switch (order) {
        case 1:
            return "SING";
        case 2:
            return "DOUB";
        case 3:
            return "TRIP";
        case 4:
            return "AROM";
        default:
            return ".";  // unknown
    }
}

static void append_charges_to_block(const Molecule &molecule, const Charges &charges, gemmi::cif::Block &block) {
    const std::string partial_atomic_charges_meta_prefix = "_partial_atomic_charges_meta.";
    const std::string partial_atomic_charges_prefix = "_partial_atomic_charges.";
    
    const std::vector<std::string> partial_atomic_charges_meta_attributes = {
        "id",
        "type",
        "method",
    };

    const std::vector<std::string> partial_atomic_charges_attributes = {
        "type_id",
        "atom_id",
        "charge",
    };

    const auto& atom_charges = charges[molecule.name()];

    // _partial_atomic_charges_meta
    auto& metadata_loop = block.init_loop(partial_atomic_charges_meta_prefix, partial_atomic_charges_meta_attributes);
    const auto id = "1";
    const auto type = "empirical";
    const auto method = fmt::format("'{}/{}'", charges.method_name(), charges.parameters_name());
    metadata_loop.add_row({
        id,
        type,
        method,
    });

    // _partial_atomic_charges  
    auto& charges_loop = block.init_loop(partial_atomic_charges_prefix, partial_atomic_charges_attributes);
    for (size_t i = 0; i < molecule.atoms().size(); ++i) {
        const auto &atom = molecule.atoms()[i];
        const auto id = "1";
        const auto atomId = fmt::format("{}", atom.index() + 1);
        const auto charge = fmt::format("{: 1.4f}", atom_charges[i]);
        charges_loop.add_row({
            id,
            atomId,
            charge,
        });
    }
}

static void filter_out_altloc_atoms(gemmi::cif::Block &block) {
    auto table = block.find_mmcif_category("_atom_site.");
    auto &loop = *table.get_loop();

    const int label_alt_id_index_col = table.find_column_position("_atom_site.label_alt_id");
    const int id_col = table.find_column_position("_atom_site.id");
    
    std::set<int> rows_to_remove;
    const auto& label_alt_id = table.bloc.find_loop(loop.tags[label_alt_id_index_col]);
    for (unsigned i = 0; i < loop.length(); ++i) {
        if (label_alt_id.at(i) != "." and label_alt_id.at(i) != "A") {
            rows_to_remove.insert(i);
        }
    }

    const size_t new_length = loop.length() - rows_to_remove.size();
    std::vector<std::vector<std::string>> new_columns(loop.tags.size());

    for (unsigned j = 0; j != loop.tags.size(); ++j) {
        auto column = table.bloc.find_loop(loop.tags[j]);
        for (unsigned i = 0; i != loop.length(); ++i) {
            if (rows_to_remove.find(i) == rows_to_remove.end()) {
                new_columns[j].push_back(column[i]);
            }
        }
    }

    // reset ids
    for (unsigned i = 0; i < new_length; ++i) {
        new_columns[id_col][i] = fmt::format("{}", i + 1);
    }

    loop.set_all_values(new_columns);
}

static void generate_mmcif_from_block(gemmi::cif::Block &block, const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    std::filesystem::path out_dir{config::chg_out_dir};
    std::string out_filename = std::filesystem::path(filename).filename().string() + ".charges.cif";
    std::string out_file{(out_dir / out_filename).string()};
    std::ofstream out_stream{out_file};

    filter_out_altloc_atoms(block);
    append_charges_to_block(ms.molecules()[0], charges, block);

    gemmi::cif::write_cif_block_to_stream(out_stream, block);
}

static void generate_mmcif_from_cif_file(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto document = gemmi::cif::read_file(filename);
    auto& block = document.sole_block();
    
    generate_mmcif_from_block(block, ms, charges, filename);
}

static void generate_mmcif_from_pdb_file(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto structure = gemmi::read_pdb_file(filename);

    if (structure.models.empty() || structure.models[0].chains.empty()) {
        fmt::print(stderr, "ERROR: No models or no chains in PDB file.\n");
        exit(EXIT_FILE_ERROR);
    }

    gemmi::setup_entities(structure);
    gemmi::assign_label_seq_id(structure, false);

    auto block = gemmi::make_mmcif_block(structure);

    // remove pesky _chem_comp category >:(
    block.find_mmcif_category("_chem_comp.").erase();

    generate_mmcif_from_block(block, ms, charges, filename);
}

static void generate_mmcif_from_atom_and_bond_data(const MoleculeSet &ms, const Charges &charges) {
    const std::string atom_site_prefix = "_atom_site.";
    const std::string chem_comp_prefix = "_chem_comp.";    
    const std::string chem_comp_bond_prefix = "_chem_comp_bond.";

    const std::vector<std::string> atom_site_attributes = {
        "group_PDB",
        "id",
        "type_symbol",
        "label_atom_id",
        "label_comp_id",
        "label_seq_id",
        // "label_alt_id",
        "label_asym_id",
        "label_entity_id",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        // "auth_asym_id",
    };

    const std::vector<std::string> chem_comp_attributes = {
        "id",
    };

    const std::vector<std::string> chem_comp_bond_attributes = {
        "comp_id",
        "atom_id_1",
        "atom_id_2",
        "value_order",
    };


    for (const auto& molecule : ms.molecules()) {
        std::string molecule_name = to_lowercase(molecule.name());    
        std::filesystem::path out_dir{config::chg_out_dir};
        std::string out_filename = molecule_name + ".charges.cif";
        std::string out_file{(out_dir / out_filename).string()};
        std::ofstream out_stream{out_file};

        auto document = gemmi::cif::Document{};
        auto& block = document.add_new_block(molecule_name);

        // _atom_site
        auto& atom_site_loop = block.init_loop(atom_site_prefix, atom_site_attributes);
        for (size_t i = 0; i < molecule.atoms().size(); i++) {
            const auto &atom = molecule.atoms()[i];
            const std::string group_PDB = atom.hetatm() ? "HETATM" : "ATOM";
            const std::string id = fmt::format("{}", atom.index() + 1);
            const std::string type_symbol = atom.element().symbol();
            const std::string label_atom_id = id;
            const std::string label_comp_id = atom.residue();
            const std::string label_seq_id = fmt::format("{}", atom.residue_id());
            // const std::string label_alt_id = ".";
            const std::string label_asym_id = atom.chain_id() == "" ? "." : atom.chain_id();
            const std::string label_entity_id = "1";
            const std::string cartn_x = fmt::format("{:.3f}", atom.pos()[0]);
            const std::string cartn_y = fmt::format("{:.3f}", atom.pos()[1]);
            const std::string cartn_z = fmt::format("{:.3f}", atom.pos()[2]);
            // const std::string auth_asym_id = ".";
            atom_site_loop.add_row({
                group_PDB,
                id,
                type_symbol,
                label_atom_id,
                label_comp_id,
                label_seq_id,
                // label_alt_id,
                label_asym_id,
                label_entity_id,
                cartn_x,
                cartn_y,
                cartn_z,
                // auth_asym_id,
            });
        }

        // _chem_comp
        auto& chem_comp_loop = block.init_loop(chem_comp_prefix, chem_comp_attributes);
        const std::string comp_id = "UNL";
        chem_comp_loop.add_row({
            comp_id,
        });

        // _chem_comp_bond
        auto& chem_comp_bond_loop = block.init_loop(chem_comp_bond_prefix, chem_comp_bond_attributes);
        for (const auto& bond: molecule.bonds()) {
            const std::string comp_id = "UNL";
            const std::string atom_id_1 = fmt::format("{}", bond.first().index() + 1);
            const std::string atom_id_2 = fmt::format("{}", bond.second().index() + 1);
            const std::string value_order = convert_bond_order_to_mmcif_value_order_string(bond.order());
            chem_comp_bond_loop.add_row({
                comp_id,
                atom_id_1,
                atom_id_2,
                value_order,
            });
        }

        append_charges_to_block(molecule, charges, block);

        gemmi::cif::write_cif_block_to_stream(out_stream, block);
    }
}

void CIF::save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto ext = std::filesystem::path(filename).extension().string();
    
    if (ext == ".cif") {
        generate_mmcif_from_cif_file(ms, charges, config::input_file);
    } else if (ext == ".pdb" or ext == ".ent") {
        generate_mmcif_from_pdb_file(ms, charges, config::input_file);
    } else if (ext == ".mol2" or ext == ".sdf") {
        generate_mmcif_from_atom_and_bond_data(ms, charges);
    }
}
