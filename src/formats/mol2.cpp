//
// Created by krab1k on 24.1.19.
//

#include <memory>
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <iterator>
#include <stdexcept>
#include <set>
#include <fmt/format.h>
#include <filesystem>
#include <gemmi/to_cif.hpp>

#include "common.h"
#include "chargefw2.h"
#include "cif.h"
#include "../config.h"
#include "../charges.h"
#include "../periodic_table.h"
#include "mol2.h"
#include "../utility/strings.h"


void Mol2::read_record(std::ifstream &file, std::string &line, std::unique_ptr<std::vector<Atom>> &atoms,
                       std::unique_ptr<std::vector<Bond>> &bonds) {

    std::getline(file, line);

    std::stringstream ss(line);
    size_t n_atoms, n_bonds;
    ss >> n_atoms >> n_bonds;

    /* Skip the rest until atom records */
    do {
        std::getline(file, line);
    } while (line != "@<TRIPOS>ATOM" and not file.eof());

    atoms->reserve(n_atoms);
    for (size_t i = 0; i < n_atoms; i++) {
        std::getline(file, line);
        std::stringstream as(line);
        size_t idx;
        std::string atom_name, atom_type;
        double x, y, z;

        size_t residue_id = 0;
        std::string residue = "UNL";
        as >> idx >> atom_name >> x >> y >> z >> atom_type;
        as >> residue_id >> residue;

        std::string element_symbol;
        auto it = atom_type.find('.');
        if (it != std::string::npos) {
            element_symbol = atom_type.substr(0, it);
        } else {
            element_symbol = atom_type;
        }

        auto element = PeriodicTable::pte().get_element_by_symbol(element_symbol);

        atoms->emplace_back(i, element, x, y, z, atom_name, residue_id, residue, "", false);
    }

    /* Read @<TRIPOS>BOND */
    std::getline(file, line);
    if (line != "@<TRIPOS>BOND") {
        throw std::runtime_error("No BOND record");
    }

    bonds->reserve(n_bonds);
    for (size_t i = 0; i < n_bonds; i++) {
        std::getline(file, line);
        std::stringstream bs(line);

        size_t idx, first, second;
        std::string type;
        bs >> idx >> first >> second >> type;

        int order;
        try {
            order = std::stoi(type);
        } catch (std::invalid_argument &) {
            /* Set bond order to 1 for aromatic and other bond types */
            order = 1;
        }
        bonds->emplace_back(&((*atoms)[first - 1]), &((*atoms)[second - 1]), order);
    }
}


void Mol2::read_until_end_of_record(std::ifstream &file) {
    std::string line;
    while (std::getline(file, line)) {
        if (line == "@<TRIPOS>MOLECULE") {
            break;
        }
    }
}


MoleculeSet Mol2::read_file(const std::string &filename) {
    std::ifstream file(filename);
    if (!file) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    auto molecules = std::make_unique<std::vector<Molecule>>();

    std::string line;
    std::string name;
    std::set<std::string> molecule_names;

    /* Skip comments or empty lines*/
    do {
        std::getline(file, line);
    } while (line != "@<TRIPOS>MOLECULE" and not file.eof());

    while (std::getline(file, line)) {
        try {
            name = sanitize_name(line);
            name = get_unique_name(name, molecule_names);
            molecule_names.insert(name);

            auto atoms = std::make_unique<std::vector<Atom>>();
            auto bonds = std::make_unique<std::vector<Bond>>();
            read_record(file, line, atoms, bonds);

            if (atoms->empty()) {
                throw std::runtime_error("No atoms were loaded");
            }

            molecules->emplace_back(name, std::move(atoms), std::move(bonds));
        }
        catch (std::exception &e) {
            fmt::print(stderr, "Error when reading {}: {}\n", name, e.what());
        }
        read_until_end_of_record(file);
    }

    return MoleculeSet(std::move(molecules));
}


void Mol2::save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto file = std::fopen(filename.c_str(), "w");
    if (!file) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    for (const auto &molecule: ms.molecules()) {

        try {
            auto chg = charges[molecule.name()];
            fmt::print(file, "@<TRIPOS>MOLECULE\n");
            fmt::print(file, "{}\n", molecule.name());
            fmt::print(file, "{} {}\n", molecule.atoms().size(), molecule.bonds().size());

            /* Try to guess if the molecule is protein or not */
            if (molecule.atoms()[0].chain_id().empty()) {
                fmt::print(file, "SMALL\n");
            } else {
                fmt::print(file, "PROTEIN\n");
            }
            fmt::print(file, "USER_CHARGES\n");
            fmt::print(file, "****\n");
            fmt::print(file, "Charges calculated by ChargeFW2 {}, method: {}\n", VERSION, charges.method_name());

            fmt::print(file, "@<TRIPOS>ATOM\n");
            for (size_t i = 0; i < molecule.atoms().size(); i++) {
                const auto &atom = molecule.atoms()[i];
                fmt::print(file, "{:>5d} {:<6s} {:>8.3f} {:>8.3f} {:>8.3f} {:s} {:>3d} {:>3s} {:>6.3f}\n",
                           i + 1, atom.name(), atom.pos()[0], atom.pos()[1], atom.pos()[2], atom.element().symbol(),
                           atom.residue_id(), atom.residue(), chg[i]);
            }

            fmt::print(file, "@<TRIPOS>BOND\n");
            for (size_t i = 0; i < molecule.bonds().size(); i++) {
                const auto &bond = molecule.bonds()[i];
                fmt::print(file, "{:>5d} {:>5d} {:>5d} {:>2d}\n",
                           i + 1, bond.first().index() + 1, bond.second().index() + 1, bond.order());
            }
        }
        catch (std::out_of_range &) {
            /* Do nothing */
        }
    }
    fclose(file);
}

static std::string convert_bond_order_to_mmcif_string(int order) {
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

void Mol2::append_charges_to_file(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    const std::string atom_site_prefix = "_atom_site";
    const std::string chem_comp_prefix = "_chem_comp";    
    const std::string chem_comp_bond_prefix = "_chem_comp_bond";

    const std::vector<std::string> atom_site_attributes = {
        ".group_PDB",
        ".id",
        ".type_symbol",
        ".label_atom_id",
        ".label_comp_id",
        ".label_seq_id",
        ".label_alt_id",
        ".label_asym_id",
        ".label_entity_id",
        ".Cartn_x",
        ".Cartn_y",
        ".Cartn_z",
        ".auth_asym_id",
    };

    const std::vector<std::string> chem_comp_attributes = {
        ".id",
    };

    const std::vector<std::string> chem_comp_bond_attributes = {
        ".comp_id",
        ".atom_id_1",
        ".atom_id_2",
        ".value_order",
    };

    const auto& molecule = ms.molecules()[0];

    std::filesystem::path out_dir{config::chg_out_dir};
    std::string out_filename = std::filesystem::path(filename).filename().replace_extension(".charges.cif").string();
    std::string out_file{(out_dir / out_filename).string()};
    std::ofstream out_stream{out_file};

    auto document = gemmi::cif::Document{};
    auto& block = document.add_new_block(to_lowercase(molecule.name()));

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
        const std::string label_alt_id = ".";
        const std::string label_asym_id = atom.chain_id() == "" ? "." : atom.chain_id();
        const std::string label_entity_id = "1";
        const std::string cartn_x = fmt::format("{:.3f}", atom.pos()[0]);
        const std::string cartn_y = fmt::format("{:.3f}", atom.pos()[1]);
        const std::string cartn_z = fmt::format("{:.3f}", atom.pos()[2]);
        const std::string auth_asym_id = ".";
        atom_site_loop.add_row({
            group_PDB,
            id,
            type_symbol,
            label_atom_id,
            label_comp_id,
            label_seq_id,
            label_alt_id,
            label_asym_id,
            label_entity_id,
            cartn_x,
            cartn_y,
            cartn_z,
            auth_asym_id,
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
        const std::string value_order = convert_bond_order_to_mmcif_string(bond.order());
        chem_comp_bond_loop.add_row({
            comp_id,
            atom_id_1,
            atom_id_2,
            value_order,
        });
    }

    gemmi::cif::write_cif_block_to_stream(out_stream, block);
}