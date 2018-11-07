//
// Created by krab1k on 23/10/18.
//

#include <utility>
#include "Atom.h"
#include "../Element.h"
#include "../PeriodicTable.h"

#include <iostream>

Atom::Atom(int idx, const Element *element, double x, double y, double z) {
    index_ = idx;
    element_ = element;
    pos_[0] = x;
    pos_[1] = y;
    pos_[2] = z;
}

std::ostream &operator<<(std::ostream &str, const Atom &atom) {
    auto pte = PeriodicTable::pte();
    return str << "Atom " << atom.element_->symbol() << " Idx: " << atom.index_;
}
