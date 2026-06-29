#pragma once

#include <string>
#include <string_view>
#include <vector>

struct net {
    std::string name;
    int idx;
    std::vector<int> connected_cells;

    net(std::string_view n, int i) : name(n), idx(i) {
        connected_cells.reserve(10);
    }

    net(std::string n, int i) : net(std::string_view(n), i) {}
};

struct cell {
    std::string name;
    int idx;
    int weight;
    std::vector<int> connected_nets;
    std::vector<int> net_pin_idx;

    cell(std::string_view n, int i, int w = 1) : name(n), idx(i), weight(w) {}

    cell(std::string n, int i, int w = 1) : cell(std::string_view(n), i, w) {}
};

struct parsed_input {
    double balance_factor = 0.0;
    std::vector<net> nets;
    std::vector<cell> cells;
};
