//
// Created by krab1k on 6.11.18.
//

#pragma once

#include <iostream>
#include <iomanip>
#include <vector>

template<typename T>
std::ostream &operator<<(std::ostream &os, const std::vector<T> &vec) {
    os << std::fixed << std::setprecision(3);
    for (const auto &val: vec) {
        os << val << " ";
    }
    return os;
}