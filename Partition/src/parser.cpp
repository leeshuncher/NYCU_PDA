#include "parser.hpp"

#include <algorithm>
#include <cctype>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <stdexcept>
#include <string_view>
#include <unordered_map>
#include <vector>

using namespace std;

static string trim_copy(string s) {
    auto not_space = [](unsigned char ch) { return !std::isspace(ch); };
    s.erase(s.begin(), std::find_if(s.begin(), s.end(), not_space));
    s.erase(std::find_if(s.rbegin(), s.rend(), not_space).base(), s.end());
    return s;
}

static inline void skip_spaces(const char*& p, const char* end) {
    while (p < end && std::isspace(static_cast<unsigned char>(*p)))
        ++p;
}

static int next_stamp(vector<int>& net_seen_stamp, int& stamp) {
    if (stamp == std::numeric_limits<int>::max()) {
        std::fill(net_seen_stamp.begin(), net_seen_stamp.end(), 0);
        stamp = 1;
    }
    return stamp++;
}

struct transparent_string_hash {
    using is_transparent = void;

    size_t operator()(string_view value) const noexcept {
        return hash<string_view>{}(value);
    }

    size_t operator()(const string& value) const noexcept {
        return (*this)(string_view(value));
    }
};

struct transparent_string_equal {
    using is_transparent = void;

    bool operator()(string_view lhs, string_view rhs) const noexcept {
        return lhs == rhs;
    }

    bool operator()(const string& lhs, const string& rhs) const noexcept {
        return lhs == rhs;
    }
};

using cell_id_map = unordered_map<string, int, transparent_string_hash, transparent_string_equal>;

static cell_id_map::iterator find_cell(cell_id_map& cell_id, string_view name) {
#if defined(__cpp_lib_generic_unordered_lookup) && __cpp_lib_generic_unordered_lookup >= 201811L
    return cell_id.find(name);
#else
    return cell_id.find(string(name));
#endif
}

static bool next_token_fast(const char*& p, const char* end, string_view& tok, bool& is_semi,
                            bool& ends_with_semi) {
    skip_spaces(p, end);
    if (p >= end)
        return false;

    if (*p == ';') {
        ++p;
        tok = {};
        is_semi = true;
        ends_with_semi = false;
        return true;
    }

    const char* start = p;
    while (p < end && !std::isspace(static_cast<unsigned char>(*p)) && *p != ';')
        ++p;

    tok = string_view(start, static_cast<size_t>(p - start));
    is_semi = false;
    ends_with_semi = (p < end && *p == ';');
    if (ends_with_semi)
        ++p;
    return true;
}

static int ensure_cell_id(cell_id_map& cell_id, vector<cell>& cells, vector<int>& net_seen_stamp,
                          string_view name) {
    auto it = find_cell(cell_id, name);
    if (it != cell_id.end())
        return it->second;

    const int cid = static_cast<int>(cells.size());
    cell_id.emplace(string(name), cid);
    cells.emplace_back(name, cid);
    net_seen_stamp.push_back(0);
    return cid;
}

static void append_net_cell(net& n, int cid, bool fast_dedup, vector<int>& net_seen_stamp,
                            int cur_stamp) {
    if (!fast_dedup) {
        n.connected_cells.push_back(cid);
        return;
    }

    if (net_seen_stamp[cid] != cur_stamp) {
        net_seen_stamp[cid] = cur_stamp;
        n.connected_cells.push_back(cid);
    }
}

static void finalize_net(parsed_input& out, net& n, bool fast_dedup) {
    if (!fast_dedup) {
        sort(n.connected_cells.begin(), n.connected_cells.end());
        n.connected_cells.erase(unique(n.connected_cells.begin(), n.connected_cells.end()),
                                n.connected_cells.end());
    }

    for (int pin = 0; pin < static_cast<int>(n.connected_cells.size()); pin++) {
        const int cid = n.connected_cells[pin];
        out.cells[cid].connected_nets.push_back(n.idx);
        out.cells[cid].net_pin_idx.push_back(pin);
    }

    out.nets.push_back(std::move(n));
}

parsed_input parse_input(const string& filename) {
    ifstream infile(filename, ios::binary);
    if (!infile.is_open())
        throw runtime_error("Error opening file: " + filename);

    infile.seekg(0, ios::end);
    const streamoff file_size = infile.tellg();
    infile.seekg(0, ios::beg);

    parsed_input out;

    string first_line;
    if (!getline(infile, first_line))
        throw runtime_error("Empty input file: " + filename);
    out.balance_factor = stod(trim_copy(first_line));

    const size_t estimated_cells =
        max<size_t>(1 << 16, static_cast<size_t>(max<streamoff>(file_size / 32, 0)));
    const size_t estimated_nets =
        max<size_t>(1 << 14, static_cast<size_t>(max<streamoff>(file_size / 48, 0)));
    const bool fast_dedup = (file_size >= (8LL << 20));
    const bool fast_parse = (file_size >= (8LL << 20));

    cell_id_map cell_id;
    cell_id.max_load_factor(0.7f);
    cell_id.reserve(estimated_cells);
    out.cells.reserve(estimated_cells);
    out.nets.reserve(estimated_nets);
    vector<int> net_seen_stamp;
    net_seen_stamp.reserve(estimated_cells);
    int stamp = 1;

    if (fast_parse) {
        string buffer(static_cast<size_t>(file_size), '\0');
        infile.seekg(0, ios::beg);
        if (!buffer.empty())
            infile.read(&buffer[0], file_size);
        if (!infile && infile.gcount() != file_size)
            throw runtime_error("Error reading file: " + filename);

        const char* p = buffer.c_str();
        const char* end = p + buffer.size();
        char* balance_end = nullptr;
        out.balance_factor = std::strtod(p, &balance_end);
        if (balance_end == p)
            throw runtime_error("Empty input file: " + filename);
        p = balance_end;

        string_view tok;
        while (true) {
            bool is_semi = false;
            bool ends_with_semi = false;
            if (!next_token_fast(p, end, tok, is_semi, ends_with_semi))
                break;
            if (is_semi || tok != "NET")
                continue;

            string_view net_name;
            if (!next_token_fast(p, end, net_name, is_semi, ends_with_semi) || is_semi)
                throw runtime_error("Malformed NET line: missing net name");

            net n(net_name, static_cast<int>(out.nets.size()));
            const int cur_stamp = fast_dedup ? next_stamp(net_seen_stamp, stamp) : 0;

            while (next_token_fast(p, end, tok, is_semi, ends_with_semi)) {
                if (is_semi)
                    break;

                const int cid = ensure_cell_id(cell_id, out.cells, net_seen_stamp, tok);
                append_net_cell(n, cid, fast_dedup, net_seen_stamp, cur_stamp);

                if (ends_with_semi)
                    break;
            }

            finalize_net(out, n, fast_dedup);
        }

        return out;
    }

    infile.clear();
    infile.seekg(0, ios::beg);
    string tok;
    while (infile >> tok) {
        if (tok != "NET")
            continue;

        string net_name;
        if (!(infile >> net_name))
            throw runtime_error("Malformed NET line: missing net name");

        net n(net_name, static_cast<int>(out.nets.size()));
        const int cur_stamp = fast_dedup ? next_stamp(net_seen_stamp, stamp) : 0;

        while (infile >> tok) {
            if (tok == ";")
                break;

            bool ends_with_semi = (!tok.empty() && tok.back() == ';');
            if (ends_with_semi)
                tok.pop_back();

            if (!tok.empty()) {
                const int cid = ensure_cell_id(cell_id, out.cells, net_seen_stamp, tok);
                append_net_cell(n, cid, fast_dedup, net_seen_stamp, cur_stamp);
            }

            if (ends_with_semi)
                break;
        }

        finalize_net(out, n, fast_dedup);
    }

    return out;
}
