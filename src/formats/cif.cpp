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

// class MCRA {
//     const int _model;
//     const std::string _chain;
//     const std::string _res_num;
//     const std::string _residue;
//     const std::string _atom;

// public:
//     MCRA(const int model,
//          const std::string &chain,
//          const std::string &res_num,
//          const std::string &residue,
//          const std::string &atom)
//         : _model(model), _chain(chain), _res_num(res_num), _residue(residue), _atom(atom)
//         {}

//     // finds _atom_site.id value corresponding to this MCRA
//     int find_atomsite_id(gemmi::cif::Table &table) const {
//         const int id_idx = table.find_column_position("_atom_site.id");

//         for (const auto row: table) {
//             if (is_row(table, row)) {
//                 return std::stoi(row[id_idx]);
//             }
//         }

//         return -1;
//     }

//     friend bool operator==(const MCRA &lhs, const MCRA &rhs) {
//         return lhs._model == rhs._model &&
//                lhs._chain == rhs._chain &&
//                lhs._res_num == rhs._res_num &&
//                lhs._residue == rhs._residue &&
//                lhs._atom == rhs._atom;
//     }

// private:

//     std::string remove_surrounding_quotes(const std::string &string) const {
//         if (string[0] != '"') {
//             return string;
//         }
//         return string.substr(1, string.size() - 2);
//     }

//     // returns true if the row is the same as this MCRA
//     bool is_row(const gemmi::cif::Table &table, const gemmi::cif::Table::Row &row) const {
//         const int pdbx_PDB_model_num_idx = table.find_column_position("_atom_site.pdbx_PDB_model_num");
//         const int auth_asym_id_idx = table.find_column_position("_atom_site.auth_asym_id");
//         const int auth_seq_id_idx = table.find_column_position("_atom_site.auth_seq_id");
//         const int label_comp_id_idx = table.find_column_position("_atom_site.label_comp_id");
//         const int label_atom_id_idx = table.find_column_position("_atom_site.label_atom_id");

//         const int model = static_cast<int>(std::stoul(row[pdbx_PDB_model_num_idx]));
//         const std::string &chain = row[auth_asym_id_idx];
//         const std::string &res_num = row[auth_seq_id_idx];
//         const std::string &residue = row[label_comp_id_idx];
//         const std::string &atom = remove_surrounding_quotes(row[label_atom_id_idx]);

//         // fmt::print(stderr, "{} {} | {} {} | {} {} | {} {} | {} {}\n",
//         //         model, _model,
//         //         chain, _chain,
//         //         res_num, _res_num,
//         //         residue, _residue,
//         //         atom, _atom);

//         return *this == MCRA(model, chain, res_num, residue, atom);
//     }
// };


// static std::map<int, double> get_id_to_charge_mapping(const Molecule &molecule, const Charges &charges, gemmi::cif::Block &block) {
//     auto table = block.find_mmcif_category("_atom_site.");
//     if (!table.ok()) {
//         fmt::print(stderr, "_atom_site category is empty\n");
//         exit(EXIT_FILE_ERROR);
//     }
//     const auto &atom_charges = charges[molecule.name()];
//     const int model = 1;
//     std::map<int, double> p_charge = {};
//     // loop over atoms in molecule and find corresponding row in table
//     // then update the partial charge and vdw radius
//     for (size_t i = 0; i < molecule.atoms().size(); ++i) {
//         const auto &atom = molecule.atoms()[i];

//         // check that the atom index is correct
//         if (atom.index() != i) {
//             fmt::print(stderr, "Atom index mismatch: {} != {}\n", atom.index(), i);
//             exit(EXIT_FILE_ERROR);
//         }

//         MCRA mcra{
//             model, 
//             atom.chain_id(), 
//             std::to_string(atom.residue_id()), 
//             atom.residue(), 
//             atom.name()
//         };

//         const int id = mcra.find_atomsite_id(table);
//         if (id == -1){
//                 fmt::print(stderr, "Failed to find Atom(1, {}, {}, {}, {})\n", 
//                 atom.chain_id(), atom.residue_id(), atom.residue(), atom.name());
//             exit(EXIT_FILE_ERROR);
//         }
//         assert(id > 0);
//         p_charge[id] = atom_charges[atom.index()];
//     }

//     return p_charge;
// }

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
    auto structure = gemmi::make_structure_from_block(block);
    if (structure.models.empty()) {
        throw std::runtime_error("Not enough information to create a structure");
    }

    // filter out the altloc atoms
    // retains the _atom_site.id order
    auto &model = structure.models[0];
    for (gemmi::Chain &chain: model.chains) {
        for (gemmi::Residue &residue: chain.residues) {
            bool hetatm = residue.het_flag == 'H';
            for (auto it = residue.atoms.begin(); it != residue.atoms.end(); ) {
                const auto& atom = *it;
                if (not atom.has_altloc() or atom.altloc == 'A') {
                    if ((not hetatm) or
                        (config::read_hetatm and residue.name != "HOH") or
                        (config::read_hetatm and not config::ignore_water)) {
                        ++it;
                    }
                } else {
                    residue.atoms.erase(it);
                }
            }
        }
    }

    gemmi::update_mmcif_block(structure, block);
}

static void generate_mmcif_from_block(gemmi::cif::Block &block, const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    const Molecule &molecule = ms.molecules()[0];
    
    filter_out_altloc_atoms(block);
    append_charges_to_block(molecule, charges, block);
    
    std::filesystem::path out_dir{config::chg_out_dir};
    std::string molecule_name = to_lowercase(molecule.name());    
    std::string out_filename = molecule_name + ".fw2.cif";
    std::string out_file{(out_dir / out_filename).string()};
    std::ofstream out_stream{out_file};

    // remove pesky _chem_comp category >:(
    block.find_mmcif_category("_chem_comp.").erase();

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
        std::filesystem::path out_dir{config::chg_out_dir};
        std::string molecule_name = to_lowercase(molecule.name());    
        std::string out_filename = molecule_name + ".fw2.cif";
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
