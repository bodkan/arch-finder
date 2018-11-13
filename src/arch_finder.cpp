#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>



// Split a single VCF line into individual tab-separated elements.
std::vector<std::string> split_line(std::string& line)
{
    std::istringstream ss(line);
    std::vector<std::string> tokens{
        std::istream_iterator<std::string>{ss},
        std::istream_iterator<std::string>{}
    };
    return tokens;
}


/* // Convert genotypes from a VCF format into an EIGENSTRAT format by applying */
/* // a series of regex substitutions on each line as a whole. */
/* int convert_gt(const std::string& s) */
/* { */
/*     std::string converted(s); */

/*     converted = std::regex_replace(converted, std::regex("0\\|0"),     "0"); */
/*     converted = std::regex_replace(converted, std::regex("0/0"),       "0"); */

/*     converted = std::regex_replace(converted, std::regex("0\\|1"),     "1"); */
/*     converted = std::regex_replace(converted, std::regex("1\\|0"),     "1"); */
/*     converted = std::regex_replace(converted, std::regex("0/1"),       "1"); */
/*     converted = std::regex_replace(converted, std::regex("1/0"),       "1"); */

/*     converted = std::regex_replace(converted, std::regex("1\\|1"),     "2"); */
/*     converted = std::regex_replace(converted, std::regex("1/1"),       "2"); */

/*     return std::stoi(converted); */
/* } */

// Convert genotypes from a VCF format into an EIGENSTRAT format by applying
// a series of regex substitutions on each line as a whole.
int convert_gt(const std::string& s)
{
    if      (s == "0|0") { return 2; }
    else if (s == "0/0") { return 2; }

    else if (s == "0|1") { return 1; }
    else if (s == "1|0") { return 1; }
    else if (s == "0/1") { return 1; }
    else if (s == "1/0") { return 1; }

    else if (s == "1|1") { return 0; }
    else                 { return 0; } // (s == "1/1")
}


int main(int argc, char* argv[])
{
    // convert command-line args to vector
    std::vector<std::string> args(argv, argv + argc);

    // find positions of command-line arguments
    int vcf_cli = 0, bed_cli = 0, arch_cli = 0, afr_cli = 0;
    for (int i = 0; i < argc; i++) {
        if (args[i] == "--vcf") {
            vcf_cli = i + 1;
        } else if (args[i] == "--bed") {
            bed_cli = i + 1;
        } else if (args[i] == "--archaics") {
            arch_cli = i + 1;
        } else if (args[i] == "--africans") {
            afr_cli = i + 1;
        } else
            continue;
    }

    // extract path to VCF and list of archaic and African individual names
    std::string vcf_path(args[vcf_cli]);
    std::string bed_path(args[bed_cli]);
    std::vector<std::string> archaics, africans;
    for (int i = arch_cli; i < afr_cli - 1; i++) archaics.push_back(args[i]);
    for (int i = afr_cli; i < argc; i++) africans.push_back(args[i]);

    if (!vcf_cli || !bed_cli || !archaics.size() || !africans.size()) {
        std::cerr << "Usage:\n\t./arch-finder --vcf <path to VCF> --archaics <archaic samples> --africans <African samples>\n\n\tArguments must be provided in the specified order.";
        return 0;
    }

    /* std::cout << vcf_path << "\n"; */
    /* std::cout << bed_path << "\n"; */
    /* for (auto x : archaics) std::cout << x << "\n"; */
    /* for (auto x : africans) std::cout << x << "\n"; */

    std::unique_ptr<std::istream> vcf_file;
    if (!std::ifstream(vcf_path)) {
        std::cerr << "File '" << vcf_path << "' not found.\n";
        return 1;
    }

    // setup boost machinery that will be used for reading gzipped input
    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
    std::ifstream file(vcf_path, std::ios_base::in | std::ios_base::binary);
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);

    vcf_file = std::make_unique<std::istream>(&in);

    std::string line;

    // extract all individual names from the last line of the VCF header
    std::vector<std::string> inds;
    while (std::getline(*vcf_file, line)) {
        if (line.find("##") == 0) continue;

        if (line.find("#") == 0) {
            auto elems = split_line(line);
            for (std::size_t i = 9; i < elems.size(); i++) {
                inds.push_back(elems[i]);
            }
            break;
        }
    }

    // extract positions of archaics and Africans in the GT fields
    std::vector<int> arch_pos, afr_pos;
    for (auto x : archaics) for (std::size_t i = 0; i < inds.size(); i++) if (x == inds[i]) arch_pos.push_back(i + 9);
    for (auto x : africans) for (std::size_t i = 0; i < inds.size(); i++) if (x == inds[i]) afr_pos.push_back(i + 9);

    if (arch_pos.size() != archaics.size()) {
        std::cerr << "Not all Neandertals given on the command-line are present in the VCF.\n";
        return 1;
    }
    if (afr_pos.size() != africans.size()) {
        std::cerr << "Not all Africans given on the command-line are present in the VCF.\n";
        return 1;
    }

    /* for (int i : arch_pos) std::cout << i << "\n"; */
    /* for (int i : afr_pos) std::cout << i << "\n"; */

    // process the genotype part of the VCF
    std::vector<std::pair<std::string, long>> positions;
    for (long i = 0; std::getline(*vcf_file, line); i++) {
        if (i % 10000 == 1) std::cout << i << " sites processed\r" << std::flush;

        auto elems = split_line(line);

        // keep only biallelic SNPs
        if (elems[3].length() != 1 || elems[4].length() != 1 || elems[4] == ".") continue;

        std::vector<int> arch_gt, afr_gt;
        for (int i : arch_pos) if (elems[i] != "./." && elems[i] != ".|.") arch_gt.push_back(convert_gt(elems[i]));
        for (int i : afr_pos) if (elems[i] != "./." && elems[i] != ".|.") afr_gt.push_back(convert_gt(elems[i]));
        if (arch_gt.size() != arch_pos.size() || !afr_gt.size()) continue;

        // calculate allele frequencies
        float arch_freq = std::accumulate(arch_gt.begin(), arch_gt.end(), 0) / float(2 * arch_gt.size());
        float afr_freq = std::accumulate(afr_gt.begin(), afr_gt.end(), 0) / float(2 * afr_gt.size());

        if (arch_freq != 1 || afr_freq != 0) continue;

        positions.push_back(std::make_pair(elems[0], std::stoul(elems[1])));

        /* std::cout << arch_freq << "\t" << afr_freq << "\t#\t"; */
        /* for (auto x : arch_gt) std::cout << x << "\t"; */
        /* for (auto x : afr_gt) std::cout << x << "\n"; */
    }

    // write coordinates
    std::ofstream bed(bed_path);
    for (const auto & pos : positions) {
        bed << pos.first << "\t" << pos.second - 1 << "\t" << pos.second << "\n";
    }
    bed.close();

    return 0;
}
