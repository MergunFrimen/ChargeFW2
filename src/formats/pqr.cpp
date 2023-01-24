//
// Created by krab1k on 31.1.19.
//

#include <string>
#include <fstream>
#include <filesystem>
#define GEMMI_WRITE_IMPLEMENTATION // TODO: add sprintf.h to usr/include/gemmi/third_party libs
#include <gemmi/to_mmcif.hpp>  // make_mmcif_document
#undef GEMMI_WRITE_IMPLEMENTATION
#include <gemmi/to_cif.hpp>
#include <gemmi/pdb.hpp>
#include <gemmi/align.hpp>  // assign_label_seq_id
#include <fmt/format.h>

#include "chargefw2.h"
#include "pqr.h"
#include "cif.h"
#include "../config.h"
#include "../structures/molecule_set.h"
#include "../charges.h"


void PQR::save_charges(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto file = std::fopen(filename.c_str(), "w");
    if (!file) {
        fmt::print(stderr, "Cannot open file: {}\n", filename);
        exit(EXIT_FILE_ERROR);
    }

    const auto &molecule = ms.molecules()[0];

    try
    {
        auto chg = charges[molecule.name()];
        for (size_t i = 0; i < molecule.atoms().size(); i++) {
            const auto &atom = molecule.atoms()[i];

            /* Match PDB format weirdness */
            std::string name;
            if (atom.element().symbol().length() == 1 and atom.name().length() < 4) {
                name = " ";
            }
            name.append(atom.name());

            fmt::print(file, "{:<6}{:>5d} {:<4s} {:>3s} {:1s} {:>3d}    {:>8.3f}{:>8.3f}{:>8.3f} {:>6.3f} {:>6.3f}\n",
                       atom.hetatm() ? "HETATM": "ATOM", i + 1, name, atom.residue(), atom.chain_id(), atom.residue_id(),
                       atom.pos()[0], atom.pos()[1], atom.pos()[2], chg[i], atom.element().vdw_radius());
        }
    }
    catch (std::out_of_range &) {
        /* Do nothing */
    }

    fclose(file);
}

// converts PDB file to CIF file
// TODO: append charges to converted CIF file
void PQR::append_charges_to_file(const MoleculeSet &ms, const Charges &charges, const std::string &filename) {
    auto structure = gemmi::read_pdb_file(filename);

    if (structure.models.empty() || structure.models[0].chains.empty()) {
        fmt::print(stderr, "ERROR: No models or no chains in PDB file.\n");
        exit(EXIT_FILE_ERROR);
    }

    gemmi::setup_entities(structure);
    gemmi::assign_label_seq_id(structure, false);

    const auto& document = gemmi::make_mmcif_document(structure);
    
    std::filesystem::path out_dir{config::chg_out_dir};
    std::string out_filename = std::filesystem::path(filename).filename().replace_extension(".charges.cif").string();
    std::string out_file{(out_dir / out_filename).string()};
    std::ofstream out_stream{out_file};
    
    write_cif_to_stream(out_stream, document);
}