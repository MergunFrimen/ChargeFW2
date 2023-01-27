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

namespace fs = std::filesystem;

static const std::vector<std::string> atom_site_columns{
        "group_PDB",
        "auth_asym_id",
        "auth_seq_id",
        "pdbx_PDB_ins_code",
        "auth_comp_id",
        "auth_atom_id",
        "label_alt_id",
        "type_symbol",
        "Cartn_x",
        "Cartn_y",
        "Cartn_z",
        "pdbx_PDB_model_num",
    };


class MCRA {
    const int _model;
    const std::string _chain;
    const std::string _res_num;
    const std::string _residue;
    const std::string _atom;

public:
    MCRA(const int model,
         const std::string &chain,
         const std::string &res_num,
         const std::string &residue,
         const std::string &atom)
        : _model(model), _chain(chain), _res_num(res_num), _residue(residue), _atom(atom)
        {}

    int find_row(gemmi::cif::Table &table, int start_idx = 0) const {
        int idx = start_idx; 
        int total_rows = static_cast<int>(table.length());

        assert(idx < total_rows);

        while (idx < total_rows) {
            if (this->is_row(table[idx])) {
                return idx;
            }
            ++idx;
        }

        fmt::print(stderr, "In MCRA::find_row(), unable to find Atom(1, {}, {}, {}, {})\n"
                    "Looping through entire table to double check...",
                    this->_chain, this->_res_num, this->_residue, this->_atom); 

        idx = 0;
        for (const auto site : table) {
            if (this->is_row(site)) {
                return idx;
            }
            ++idx;
        }
        return -1;
    }

    friend std::ostream &operator<<(std::ostream &os, const MCRA &mcra) {
        os << "Atom(" << mcra._model << " " << mcra._chain << " " << mcra._res_num << " [" << mcra._residue << "] " << mcra._atom << ")";
        return os;
    }

private:
    bool is_row(const gemmi::cif::Table::Row &row) const {
        int model = static_cast<int>(std::stoul(row[11]));
        const std::string &chain = row[1];
        const std::string &residue = row[4];
        const std::string &res_num = row[2];
        const std::string atom = row[5][0] != '"' ? row[5] : row[5].substr(1, row[5].size() - 2);

        return _model == model && _chain == chain && _res_num == res_num && _residue == residue && _atom == atom;
    }
};


void CIF::append_fw2_config(gemmi::cif::Block &block, const std::string &method, const std::string &parameters) {

    std::string config_prefix = "_chargeFW2_config.";

    std::vector<std::string> config_tags{
        "method", "parameter_set", 
        "read_hetam", "ignore_water",
        "permissive_types"
    };

    std::string read_hetatm = config::read_hetatm ? "True": "False";
    std::string ignore_water = config::ignore_water ? "True": "False";
    std::string permissive_types = config::permissive_types ? "True": "False";

    std::vector<std::string> config_data{
        fmt::format("\"{}\"", method), fmt::format("\"{}\"", parameters),
        read_hetatm, ignore_water,
        permissive_types
    };

    for (unsigned i = 0; i != config_tags.size(); ++i) {
            block.set_pair(config_prefix + config_tags[i], config_data[i]);
    }
}


void CIF::replace_fw2_columns(gemmi::cif::Table &table, 
                              const std::vector<std::string> &p_charge,
                              const std::vector<std::string> &vdw_radii,
                              const std::vector<std::string> &fw2_tags) {

    std::vector<std::vector<std::string>> fw2_columns = {p_charge, vdw_radii};

    assert(fw2_columns.size() == fw2_tags.size());

    for (unsigned i = 0; i != fw2_tags.size(); ++i){
        auto column = table.bloc.find_loop(fw2_tags[i]);
        std::copy(fw2_columns[i].begin(), fw2_columns[i].end(), column.begin());
    }
}


void CIF::append_fw2_columns(gemmi::cif::Table &table,
                             const std::vector<std::string> &p_charge,
                             const std::vector<std::string> &vdw_radii,
                             const std::vector<std::string> &fw2_tags) {

    auto &loop = *table.get_loop();

    unsigned long orig_tag_size = loop.tags.size();
    unsigned long new_tag_size = orig_tag_size + 2;

    // Creates a new table full of empty strings with the correct number of dimensions
    // Outside vector size is the # of columns, inside vector size is the # of rows.
    std::vector<std::vector<std::string>> new_columns(new_tag_size, {loop.length(), {"Empty"}});

    // Copies data from original columns to their respecitve column in the new table filled with empty strings.
    // Leaving only the new appended columns as empty strings
    for (unsigned i = 0; i != orig_tag_size; ++i) {
        auto column = table.bloc.find_loop(loop.tags[i]);
        std::copy(column.begin(), column.end(), new_columns[i].begin());
    }

    new_columns[new_tag_size - 2] = p_charge;
    new_columns[new_tag_size - 1] = vdw_radii;

    for (const auto &tag : fw2_tags)
        loop.tags.push_back(tag);

    loop.set_all_values(new_columns);
}                                


void CIF::write_cif_block(std::ostream &out,
                          gemmi::cif::Table &table, 
                          const std::vector<std::string> &p_charge,
                          const std::vector<std::string> &vdw_radii,
                          const std::string &method,
                          const std::string &parameters) {

    append_fw2_config(table.bloc, method, parameters);

    std::vector<std::string> fw2_tags{
        "_atom_site.fw2_charge",
        "_atom_site.fw2_vdw_radius"};

    auto &loop = *table.get_loop();

    if (loop.has_tag(fw2_tags[0]) && loop.has_tag(fw2_tags[1])){
        replace_fw2_columns(table, p_charge, vdw_radii, fw2_tags);
    } else {
        append_fw2_columns(table, p_charge, vdw_radii, fw2_tags);
    }

    gemmi::cif::write_cif_block_to_stream(out, table.bloc);
}


void CIF::save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {

    fs::path out_dir{config::chg_out_dir};
    std::string out_filename = fs::path(filename).filename().replace_extension(".fw2.cif").string();
    std::string out_file{(out_dir / out_filename).string()};
    std::ofstream out_stream{out_file};

    auto doc = gemmi::cif::read_file(filename);
    auto& block = doc.sole_block();
    auto table = block.find("_atom_site.", atom_site_columns);

    if (!table.ok()) {
        fmt::print(stderr, "_atom_site category is empty for file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    auto p_charge  = std::vector<std::string>{table.length(), "?"};
    auto vdw_radii = std::vector<std::string>{table.length(), "?"};

    const auto &molecule = ms.molecules()[0];

    // ChargeFW2 is hardcoded to only read first model.
    const int model = 1;
    int row_num = 0;

    try {
        auto chg = charges[molecule.name()];
        for (size_t i = 0; i < molecule.atoms().size(); i++) {
            const auto &atom = molecule.atoms()[i];

            MCRA mcra{
                model, 
                atom.chain_id(), 
                std::to_string(atom.residue_id()), 
                atom.residue(), 
                atom.name()
            };

            row_num = mcra.find_row(table, row_num);

            if (row_num == -1){
                 fmt::print(stderr, "Failed to find Atom(1, {}, {}, {}, {})\n", 
                    atom.chain_id(), atom.residue_id(), atom.residue(), atom.name());
                exit(EXIT_FILE_ERROR);
            }
            p_charge[row_num]  = fmt::format("{:.3f}", chg[i]);
            vdw_radii[row_num] = fmt::format("{:.3f}", atom.element().vdw_radius());
        }
        write_cif_block(out_stream, table, p_charge, vdw_radii, charges.method_name(), charges.parameters_name());
    }
    catch (std::out_of_range &) {
        /* Do nothing */
    }
}

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

static void append_charges_to_block(const MoleculeSet &ms, const Charges &charges, gemmi::cif::Block &block) {
    const std::string partial_atomic_charges_meta_prefix = "_partial_atomic_charges_meta";
    const std::string partial_atomic_charges_prefix = "_partial_atomic_charges";
    
    const std::vector<std::string> partial_atomic_charges_meta_attributes = {
        ".id",
        ".type",
        ".method",
    };

    const std::vector<std::string> partial_atomic_charges_attributes = {
        ".type_id",
        ".atom_id",
        ".charge",
    };

    const auto& molecule = ms.molecules()[0];
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
        const auto charge = fmt::format("{}", atom_charges[i]);
        charges_loop.add_row({
            id,
            atomId,
            charge,
        });
    }
}

static void generate_from_cif_file(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    std::filesystem::path out_dir{config::chg_out_dir};
    std::string out_filename = std::filesystem::path(filename).filename().string() + ".charges.cif";
    std::string out_file{(out_dir / out_filename).string()};
    std::ofstream out_stream{out_file};

    auto document = gemmi::cif::read_file(filename);
    auto& block = document.sole_block();

    append_charges_to_block(ms, charges, block);

    gemmi::cif::write_cif_block_to_stream(out_stream, block);
}

static void generate_from_pdb_file(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto structure = gemmi::read_pdb_file(filename);

    if (structure.models.empty() || structure.models[0].chains.empty()) {
        fmt::print(stderr, "ERROR: No models or no chains in PDB file.\n");
        exit(EXIT_FILE_ERROR);
    }

    gemmi::setup_entities(structure);
    gemmi::assign_label_seq_id(structure, false);

    auto block = gemmi::make_mmcif_block(structure);
    
    std::filesystem::path out_dir{config::chg_out_dir};
    std::string out_filename = std::filesystem::path(filename).filename().string() + ".charges.cif";
    std::string out_file{(out_dir / out_filename).string()};
    std::ofstream out_stream{out_file};

    append_charges_to_block(ms, charges, block);
    
    gemmi::cif::write_cif_block_to_stream(out_stream, block);
}

static void generate_from_atom_and_bond_data(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
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
        // ".label_alt_id",
        ".label_asym_id",
        ".label_entity_id",
        ".Cartn_x",
        ".Cartn_y",
        ".Cartn_z",
        // ".auth_asym_id",
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
    std::string out_filename = std::filesystem::path(filename).filename().string() + ".charges.cif";
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

    append_charges_to_block(ms, charges, block);

    gemmi::cif::write_cif_block_to_stream(out_stream, block);
}

void CIF::generate_mmcif_file_with_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto ext = std::filesystem::path(filename).extension().string();
    
    if (ext == ".cif") {
        generate_from_cif_file(ms, charges, config::input_file);
    } else if (ext == ".pdb" or ext == ".ent") {
        generate_from_pdb_file(ms, charges, config::input_file);
    } else if (ext == ".mol2" or ext == ".sdf") {
        generate_from_atom_and_bond_data(ms, charges, config::input_file);
    }
}
