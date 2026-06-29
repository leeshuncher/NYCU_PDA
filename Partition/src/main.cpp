#include <array>
#include <chrono>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "parser.hpp"
#include "partitioner.hpp"

using namespace std;

static void write_group(ofstream& out, const string& label, const vector<cell>& cells,
                        const vector<int>& cell_ids) {
    out << label << ' ' << cell_ids.size() << "\n";
    for (int cid : cell_ids)
        out << cells[cid].name << ' ';
    out << ";\n";
}

int main(int argc, char** argv) {
    ios::sync_with_stdio(false);
    cin.tie(nullptr);

    if (argc != 3) {
        cerr << "Usage: ./Lab1 <input file name> <output file name>\n";
        return 1;
    }

    const string inputFile = argv[1];
    const string outputFile = argv[2];

    try {
        auto parsed = parse_input(inputFile);
        const auto deadline = chrono::steady_clock::now() + chrono::seconds(3000);

        auto result = partition(parsed, deadline);

        array<vector<int>, 3> groups;
        const array<string, 3> labels = {"G1", "G2", "G3"};
        for (auto& group : groups)
            group.reserve(parsed.cells.size());

        for (int cid = 0; cid < static_cast<int>(parsed.cells.size()); cid++) {
            const int gid = static_cast<int>(result.part[cid]);
            if (gid >= 0 && gid < static_cast<int>(groups.size()))
                groups[gid].push_back(cid);
        }

        ofstream out(outputFile);
        if (!out.is_open()) {
            cerr << "Error opening output file: " << outputFile << "\n";
            return 1;
        }

        out << "Cutsize = " << result.cut << "\n";
        for (int gid = 0; gid < static_cast<int>(groups.size()); gid++)
            write_group(out, labels[gid], parsed.cells, groups[gid]);

        const int s1 = static_cast<int>(groups[0].size());
        const int s2 = static_cast<int>(groups[1].size());
        const int s3 = static_cast<int>(groups[2].size());
        auto in_range = [&](int s) { return result.minSize <= s && s <= result.maxSize; };
        const bool balance_ok = in_range(s1) && in_range(s2) && in_range(s3);
        cerr << "Cutsize = " << result.cut << "\n";
        cerr << "Balance " << (balance_ok ? "PASS" : "FAIL")
             << " (constraint_min=" << result.minSize << ", constraint_max=" << result.maxSize
             << ", G1=" << s1 << ", G2=" << s2 << ", G3=" << s3 << ")\n";
    } catch (const exception& e) {
        cerr << e.what() << "\n";
        return 1;
    }

    return 0;
}
