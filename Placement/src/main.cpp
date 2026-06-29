#include "Legalizer.h"

#include <iostream>
#include <string>

int main(int argc, char **argv) {
    if (argc != 5) {
        std::cerr << "Usage: " << argv[0] << " <alpha> <threshold> <input_file> <output_file>\n";
        return 1;
    }

    double alpha = std::stod(argv[1]);
    double threshold = std::stod(argv[2]);
    std::string inputPath = argv[3];
    std::string outputPath = argv[4];

    Legalizer legalizer;
    std::cout << "[Stage] Reading input: " << inputPath << '\n';
    if (!legalizer.readInput(inputPath)) return 1;
    std::cout << "[Stage] Legalizing placement\n";
    if (!legalizer.legalize(alpha, threshold)) return 1;
    std::cout << "[Stage] Writing output: " << outputPath << '\n';
    if (!legalizer.writeOutput(outputPath)) return 1;
    return 0;
}
