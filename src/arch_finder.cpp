#include <iostream>
#include <sstream>
#include <fstream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/filter/gzip.hpp>

const int GT_POS = 9;

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



// Convert VCF genotypes into ALT allele counts.
int convert_gt(const std::string& s)
{
    if      (s == "0|0") { return 0; }
    else if (s == "0/0") { return 0; }

    else if (s == "0|1") { return 1; }
    else if (s == "1|0") { return 1; }
    else if (s == "0/1") { return 1; }
    else if (s == "1/0") { return 1; }

    else if (s == "1|1") { return 2; }
    else if (s == "1/1") { return 2; }
    else throw std::runtime_error(std::string("Unsupported genotype: ") + s);
}



bool missing_samples(auto& outgroups, auto& out_pos,
                     auto& archaics, auto& arch_pos,
                     auto& africans, auto& afr_pos) {
    int missing = 0;
    if (out_pos.size() != outgroups.size()) {
        std::cerr << "Not all outgroup samples are present in the VCF.\n";
        missing++;
    }
    if (arch_pos.size() != archaics.size()) {
        std::cerr << "Not all archaics are present in the VCF.\n";
        missing++;
    }
    if (afr_pos.size() != africans.size()) {
        std::cerr << "Not all Africans are present in the VCF.\n";
        missing++;
    }
    return missing > 0;
}



int main(int argc, char* argv[])
{
    // convert command-line args to vector
    std::vector<std::string> args(argv, argv + argc);

    // find positions of command-line arguments
    int vcf_cli = 0, bed_cli = 0, out_cli = 0, arch_cli = 0, afr_cli = 0;
    for (int i = 0; i < argc; i++) {
        if      (args[i] == "--vcf")       vcf_cli = i + 1;
        else if (args[i] == "--bed")       bed_cli = i + 1;
        else if (args[i] == "--outgroups") out_cli = i + 1;
        else if (args[i] == "--archaics")  arch_cli = i + 1;
        else if (args[i] == "--africans")  afr_cli = i + 1;
        else continue;
    }

    // extract path to VCF and lists of outgroup, archaic and African sample names
    std::string vcf_path(args[vcf_cli]), bed_path(args[bed_cli]);
    std::vector<std::string> outgroups, archaics, africans;
    for (int i = out_cli; i < arch_cli - 1; i++) outgroups.push_back(args[i]);
    for (int i = arch_cli; i < afr_cli - 1; i++) archaics.push_back(args[i]);
    for (int i = afr_cli; i < argc; i++)         africans.push_back(args[i]);

    if (!vcf_cli || !bed_cli || !out_cli || !outgroups.size() || !archaics.size() || !africans.size()) {
        std::cerr << "Usage:\n\t./arch-finder --vcf <path to VCF> --bed <path to BED> --outgroups <names> --archaics <names> --africans <names>\n\n\tArguments must be provided in the specified order.\n";
        return 0;
    }

    std::clog << "Input VCF: " << vcf_path << "\n";
    std::clog << "Output BED: " << bed_path << "\n";
    std::clog << "outgroups: "; for (auto x : outgroups) std::clog << x << " ";
    std::clog << "\narchaics: ";  for (auto x : archaics) std::clog << x << " ";
    std::clog << "\nAfricans: ";  for (auto x : africans) std::clog << x << " ";
    std::clog << "\n";

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
            for (std::size_t i = GT_POS; i < elems.size(); i++) {
                inds.push_back(elems[i]);
            }
            break;
        }
    }

    // extract positions of archaics and Africans in the GT fields
    std::vector<int> out_pos, arch_pos, afr_pos;
    for (auto x : outgroups) for (std::size_t i = 0; i < inds.size(); i++) if (x == inds[i]) out_pos.push_back(i + GT_POS);
    for (auto x : archaics) for (std::size_t i = 0; i < inds.size(); i++) if (x == inds[i]) arch_pos.push_back(i + GT_POS);
    for (auto x : africans) for (std::size_t i = 0; i < inds.size(); i++) if (x == inds[i]) afr_pos.push_back(i + GT_POS);

    if (missing_samples(outgroups, out_pos, archaics, arch_pos, africans, afr_pos)) return 1;

    // process the genotype part of the VCF
    std::vector<std::tuple<std::string, long, double, double, double>> sites;
    for (long i = 0; std::getline(*vcf_file, line); i++) {
        auto elems = split_line(line);

        // keep only biallelic SNPs
        if (elems[3].length() != 1 || elems[4].length() != 1 || elems[4] == "." || elems[4] == "-") continue;

        std::vector<int> out_gt, arch_gt, afr_gt;
        for (int i : out_pos) if (elems[i] != "./." && elems[i] != ".|.") out_gt.push_back(convert_gt(elems[i]));
        for (int i : arch_pos) if (elems[i] != "./." && elems[i] != ".|.") arch_gt.push_back(convert_gt(elems[i]));
        for (int i : afr_pos) if (elems[i] != "./." && elems[i] != ".|.") afr_gt.push_back(convert_gt(elems[i]));
        if (out_gt.size() != out_pos.size() || arch_gt.size() != arch_pos.size() || !afr_gt.size()) continue;

        // calculate allele frequencies
        float out_freq = std::accumulate(out_gt.begin(), out_gt.end(), 0) / float(2 * out_gt.size());
        float arch_freq = std::accumulate(arch_gt.begin(), arch_gt.end(), 0) / float(2 * arch_gt.size());
        float afr_freq = std::accumulate(afr_gt.begin(), afr_gt.end(), 0) / float(2 * afr_gt.size());

        // keep only sites at which outgroups == Africans == 0 (reference state) and
        // archaics are fixed for the derived state
        if (out_freq != 0 || arch_freq != 1 || afr_freq != 0) continue;

        sites.push_back(std::make_tuple(elems[0], std::stoul(elems[1]), out_freq, arch_freq, afr_freq));

        /* std::cout << elems[0] << "\t" << std::stoul(elems[1]) - 1 << "\t" << std::stoul(elems[1]) << "\t" << */
        /*              anc << "\t" << der << "\t" << */
        /*              out_freq << "\t" << arch_freq << "\t" << afr_freq << "\n"; */

        /* for (auto x : out_gt) std::cout << x << "\t"; */
        /* for (auto x : arch_gt) std::cout << x << "\t"; */
        /* for (auto x : afr_gt) std::cout << x << "\n"; */

        /* if (sites.size() > 100) break; */
    }

    // write coordinates
    std::ofstream bed(bed_path);
    for (const auto & s : sites) {
        bed << std::get<0>(s) << "\t" <<      // chrom
               std::get<1>(s) - 1 << "\t" <<  // start
               std::get<1>(s) << "\t" <<      // end
               std::get<2>(s) << "\t" <<      // outgroup freq
               std::get<3>(s) << "\t" <<      // archaic freq
               std::get<4>(s) << "\n";        // African freq
    }
    bed.close();

    return 0;
}
