#include <chrono>
#include <iostream>
#include <string>

#include "boost/program_options.hpp"
#include "qrot/decomposition.h"
#include "qrot/diophantine.h"
#include "qrot/gate.h"
#include "qrot/grid_solver.h"
#include "qrot/matrix.h"
#include "qrot/number.h"
#include "qrot/parser.h"

using namespace qrot;

std::string GridSynth(const AST& ast, const std::uint32_t digits) {
    namespace ch = std::chrono;
    using Ms = ch::duration<double, std::milli>;

    const auto t1 = ch::high_resolution_clock::now();
    const auto theta = ast.Value();
    auto grid_solver =
        TwoDimGridSolver::New(-theta / Float{2}, Float("1e-" + std::to_string(digits)));
    auto diophantine = Diophantine();
    auto decomposer = UnitaryDecomposer();

    const auto t2 = ch::high_resolution_clock::now();
    grid_solver.EnumerateAllSolutions();

    const auto t3 = ch::high_resolution_clock::now();
    auto solutions = std::vector<std::pair<CD2, CD2>>();
    while (solutions.empty()) {
        const auto& grid_solutions = grid_solver.GetSolutions();
        for (const auto& u : grid_solutions) {
            auto t = CD2();
            const auto xi = D2(1) - (u * u.Adj()).Real();
            const auto solved = diophantine.Solve(xi, t);
            if (solved) {
                solutions.emplace_back(u, t);
                break;
            }
        }
        if (!solutions.empty()) {
            break;
        } else {
            grid_solver.EnumerateNextLevelAllSolutions();
        }
    }

    const auto t4 = ch::high_resolution_clock::now();
    auto output = Gate();
    auto t_count = std::numeric_limits<std::size_t>::max();
    for (const auto& [u, t] : solutions) {
        auto mat = MCD2(u, -t.Adj(), t, u.Adj());
        const auto tmp = decomposer.Decompose(mat);
        std::cout << "TCount = " << tmp.CountT() << std::endl;
        if (tmp.CountT() < t_count) {
            t_count = tmp.CountT();
            output = tmp;
        }
    }

    const auto t5 = ch::high_resolution_clock::now();
#ifdef QROT_VERBOSE
    std::cout << "----------------------------------" << std::endl;
    std::cout << "Total Elapsed time = " << ch::duration_cast<Ms>(t5 - t1) << std::endl;
    std::cout << "  Setup            = " << ch::duration_cast<Ms>(t2 - t1) << std::endl;
    std::cout << "  GridProblem      = " << ch::duration_cast<Ms>(t3 - t2) << std::endl;
    std::cout << "  Diophantine      = " << ch::duration_cast<Ms>(t4 - t3) << std::endl;
    std::cout << "  Decompose        = " << ch::duration_cast<Ms>(t5 - t4) << std::endl;
#endif

    return output.ToString();
}

int main(int argc, char** argv) {
    namespace po = boost::program_options;

    // Define description
    // clang-format off
    auto description = po::options_description("Approximate z-rotation for arbitrary precision");
    description.add_options()
        ("help,h", "Display available options")
        ("theta", po::value<std::string>(), "z-rotation angle")
        ("digits,d", po::value<std::uint32_t>()->default_value(10), "Set precision in decimal digits")
    ; // NOLINT
    // clang-format on

    auto p = po::positional_options_description();
    p.add("theta", -1);

    auto vm = po::variables_map();
    po::store(po::command_line_parser(argc, argv).options(description).positional(p).run(), vm);
    po::notify(vm);

    if (vm.count("help") > 0) {
        std::cout << description << std::endl;
        return 0;
    }
    auto ast = AST();
    if (vm.count("theta") > 0) {
        const auto str = vm["theta"].as<std::string>();
        try {
            ast = AST::Parse(str);
        } catch (std::exception& ex) {
            std::cerr << "Failed to parse theta: " << str << std::endl;
            std::cerr << "    Error message: " << ex.what() << std::endl;
            std::cerr << "Examples of z-rotation angle: 1.5*pi, -pi/128, 0.56" << std::endl;
            return 1;
        }
    } else {
        std::cerr << "Z-rotation angle is not set" << std::endl;
        return 1;
    }
    const auto digits = vm["digits"].as<std::uint32_t>();

    const auto output = GridSynth(ast, digits);
    std::cout << output << std::endl;

    return 0;
}
