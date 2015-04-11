#include "generation.hpp"

#include <exception>
#include <iostream>
#include <string>

int main(int argc, char** argv) {
    bool quick = false;
    std::string out_dir = "output";

    for (int i = 1; i < argc; ++i) {
        const std::string arg(argv[i]);
        if (arg == "--quick") {
            quick = true;
        } else if (arg == "--out" && i + 1 < argc) {
            out_dir = argv[++i];
        } else if (arg == "--help") {
            std::cout << "Usage: qtoa [--quick] [--out DIR]\n";
            return 0;
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return 2;
        }
    }

    try {
        qtoa::generate_all(out_dir, quick);
    } catch (const std::exception& e) {
        std::cerr << "qtoa: " << e.what() << "\n";
        return 1;
    }

    std::cout << "Done. CSV, SVG, and verification files are in " << out_dir
              << "\n";
    return 0;
}
