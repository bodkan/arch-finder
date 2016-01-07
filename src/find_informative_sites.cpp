#include <iostream>
#include <cstdlib>

#include "VcfFileReader.h"

//
// Check if Africans share the same genotype at this site.
//
bool
fixed_in_africa(VcfRecord& rec)
{
    int afr_allele = rec.getGT(0, 0);

    // since only African samples were loaded from the VCF file previously,
    // scan through all samples in the record and test if all of them are
    // homozygous for the same allele
    for (int i = 0; i < rec.getNumSamples(); i++) {
        int allele_1 = rec.getGT(i, 0);
        int allele_2 = rec.getGT(i, 1);

        if (! ((allele_1 == afr_allele) &&
               (allele_2 == afr_allele || allele_2 == VcfGenotypeSample::INVALID_GT)))
            return false;
    }

    return true;
}

//
// Does this position carry a simple single-nucleotide allele (as
// opposed to an indel or a structural change)?
//
bool
has_simple_allele(VcfRecord& rec)
{
    return ((std::strlen(rec.getRefStr()) == 1) &&
            (std::strlen(rec.getAltStr()) == 1));
}

//
// Get a vector of positions where all Africans are fixed together with their
// allele states at these positions.
//
std::vector<std::pair<int, char>>
get_fixed_afr_sites(const std::string& hg1k_filename, const char* afr_list)
{
    VcfHeader hg1k_header;
    VcfFileReader hg1k_vcf;
    VcfRecord hg1k_rec;

    hg1k_vcf.open(hg1k_filename.c_str(), hg1k_header, afr_list, NULL, NULL);

    std::vector<std::pair<int, char>> fixed_sites;
    while (hg1k_vcf.readRecord(hg1k_rec)) {
        // consider only alleles which
        //   * are biallelic SNPs
        //   * are fixed in all Africans,
        if (has_simple_allele(hg1k_rec)
                && fixed_in_africa(hg1k_rec)) {
            int pos = hg1k_rec.get1BasedPosition();

            // get index of an African allele
            short afr_allele_index = hg1k_rec.getGT(0, 0);

            // get alleles of an African genotype and an alternative allele
            char afr_allele = *hg1k_rec.getAlleles(afr_allele_index);

            fixed_sites.push_back(std::make_pair(pos, afr_allele));
        }
    }

    return fixed_sites;
}

//
// Find a variant at a given position in a given VCF file.
// Return true/false when it has/has not been found and set the output
// arguments rec and skip_reading accordingly.
//
bool
find_matching_variant(VcfFileReader& vcf,
                      VcfRecord& rec,
                      const int position,
                      bool& skip_reading)
{
    bool variant_found = false;

    while (skip_reading || vcf.readRecord(rec)) {
        // a matching record has been found
        if (rec.get1BasedPosition() == position) {
            skip_reading = false;
            variant_found = true;
            break;
        }

        // a matching record is not present
        if (rec.get1BasedPosition() > position) {
            skip_reading = true;
            variant_found = false;
            break;
        }

        // a matching record has not yet been found
        if (rec.get1BasedPosition() < position) {
            skip_reading = false;
            continue;
        }
    }

    return variant_found;
}

//
// Get a vector of coordinates where a given archaic differs from Africans
// at sites where all Africans are fixed.
//
std::map<int, std::pair<char, char>>
check_archaic_states(const std::string& vcf_filename,
                     const short sample_id,
                     const std::vector<std::pair<int, char>>& afr_fixed_sites)
{
    VcfHeader header;
    VcfFileReader vcf;
    VcfRecord rec;
    vcf.open(vcf_filename.c_str(), header);

    // a vector where all valid positions will be accumulated
    std::map<int, std::pair<char, char>> result;

    bool skip_reading = false;
    // go through all positions which carry alleles fixed in Africans and look
    // for matching records in a given VCF file
    for (auto& afr_site : afr_fixed_sites) {

        int afr_pos = afr_site.first;
        char afr_allele = afr_site.second;

        if (find_matching_variant(vcf, rec, afr_pos, skip_reading)) {
            // to include this variant for further analysis, archaic
            //   * has to carry a valid allele (biallelic SNP),
            //   * has to be homozygous at this site,
            //   * can't have a missing genotype,
            //   * must be different than all Africans
	        if (((rec.getNumAlts() == 0) || (rec.getNumAlts() == 1))
                    && has_simple_allele(rec)
                    && (rec.getGT(sample_id, 0) == rec.getGT(sample_id, 1))
                    && (rec.getGT(sample_id, 0) != VcfGenotypeSample::MISSING_GT)
                    && (*rec.getAlleles(rec.getGT(sample_id, 0)) != afr_allele)) {
                // add this position to the final list of sites
                result.emplace(afr_pos, std::make_pair(
                    afr_allele,
                    *rec.getAlleles(rec.getGT(sample_id, 0)))
                );
            }
        }
    }

    return result;
}

//
// Get a set of sites which at which Altai and Denisovan and Africans
// all differ from each other (i.e. triallelic sites that haven't been
// detected earlier by pairwise comparisons).
//
std::set<int>
get_triallelic_positions(std::map<int, std::pair<char, char>> & altai,
                         std::map<int, std::pair<char, char>> & denisovan)
{
    std::set<int> triallelic_positions;

    // iterate through all informative sites from Altai...
    for (auto & altai_site : altai) {
      int pos = altai_site.first;
      char altai_allele = altai_site.second.second;
      // ... and check if Denisovan (if also carries an informative allele)
      // has the same allele as Altai
      if ((denisovan.count(pos) > 0) && (altai_allele != denisovan[pos].second))
          triallelic_positions.insert(pos);
    }

    return triallelic_positions;
}

//
// Filter out sites that are known not to be truly biallelic.
//
void
filter_out_triallelic(std::map<int, std::pair<char, char>> & sites,
                      std::set<int> & positions_to_remove)
{
    for (int pos : positions_to_remove)
        if (sites.count(pos)) sites.erase(pos);
}

int
main(int argc, char** argv)
{
    if (argc != 5) {
        std::cout << "Usage\n\t./find_informative_sites chr_ID 1000G_VCF Altai_VCF Denisovan_VCF\n";
        return 0;
    }

    std::string chr(argv[1]);

    std::string hg1k_vcf_file(argv[2]);
    std::string altai_vcf_file(argv[3]);
    std::string denisovan_vcf_file(argv[4]);

    std::clog << "[Chromosome " << chr << "] Started scanning the 1000 genomes VCF file.\n";
    std::vector<std::pair<int, char>> fixed_sites = get_fixed_afr_sites(hg1k_vcf_file, "tmp/afr_samples.list");
    std::clog << "[Chromosome " << chr << "] Analysis of the 1000 genomes VCF file DONE (" << fixed_sites.size() << " sites)!\n";

    std::map<int, std::pair<char, char>> altai, denisovan;

    std::clog << "[Chromosome " << chr << "] Started scanning the Altai VCF file.\n";
    altai = check_archaic_states(altai_vcf_file, 0, fixed_sites);
    std::clog << "[Chromosome " << chr << "] Analysis of the Altai VCF file DONE (" << altai.size() << " sites)!\n";

    std::clog << "[Chromosome " << chr << "] Started scanning the Denisovan VCF file.\n";
    denisovan = check_archaic_states(denisovan_vcf_file, 0, fixed_sites);
    std::clog << "[Chromosome " << chr << "] Analysis of the Denisovan VCF file DONE (" << denisovan.size() << " sites)!\n";

    // get a set of sites at which Altai, Denisovan and Africans
    // all differ from each other
    std::set<int> positions_to_remove = get_triallelic_positions(altai, denisovan);

    // remove such triallelic sites
    filter_out_triallelic(altai, positions_to_remove);
    filter_out_triallelic(denisovan, positions_to_remove);

    // initialize the table of final results
    std::map<int, std::tuple<char, char, bool, bool, bool, bool>> table;
    for (auto & site : altai) table.emplace(site.first, std::make_tuple(site.second.first, site.second.second, false, false, false, false));
    for (auto & site : denisovan) table.emplace(site.first, std::make_tuple(site.second.first, site.second.second, false, false, false, false));

    // set boolean-flag at positions where an archaic differs from Africans
    for (auto & site : altai)     std::get<2>(table[site.first]) = true;
    for (auto & site : denisovan) std::get<3>(table[site.first]) = true;

    std::clog << "[Chromosome " << chr << "] Printing out the results.\n";
    // print out the final table in a BED-like format
    for (auto & site : table) {
        std::cout << chr << "\t" << site.first - 1 << "\t" << site.first << "\t"
            << std::get<0>(site.second) << "\t"
            << std::get<1>(site.second) << "\t"
            << std::get<2>(site.second) << "\t"
            << std::get<3>(site.second) << "\n";
    }
    std::clog << "[Chromosome " << chr << "] Printing out the results DONE " <<
        "(" << table.size() << " sites in total)!\n";

    return 0;
}
