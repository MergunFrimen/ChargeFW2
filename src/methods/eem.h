//
// Created by krab1k on 31/10/18.
//

#pragma once

#include <Eigen/Core>
#include <vector>
#include <boost/config.hpp>

#include "../structures/molecule.h"
#include "../method.h"


class EEM : public EEMethod {
    enum common {kappa};
    enum atom {A, B};
    [[nodiscard]] Eigen::VectorXd EE_system(const std::vector<const Atom *> &atoms, double total_charge) const;

    [[nodiscard]] std::vector<double> calculate_charges(const Molecule &molecule) const override;

public:
    explicit EEM() : EEMethod("EEM", {"kappa"}, {"A", "B"}, {}, {}) {}

    virtual ~EEM() = default;
};

extern "C" BOOST_SYMBOL_EXPORT EEM method;
EEM method;
