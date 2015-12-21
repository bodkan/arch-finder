#include <iostream>
#include <cstdlib>

#include "VcfFileReader.h"

//
// Check that a given variant is not fixed in the global population.
//
bool
polymorphic_outside_africa(VcfRecord& rec)
{
    double amr_freq = std::stod(*(rec.getInfo().getString("AMR_AF")));
    double eas_freq = std::stod(*(rec.getInfo().getString("EAS_AF")));
    double eur_freq = std::stod(*(rec.getInfo().getString("EUR_AF")));
    double sas_freq = std::stod(*(rec.getInfo().getString("SAS_AF")));
    return ((0 < amr_freq && amr_freq < 1.0)
              || (0 < eas_freq && eas_freq < 1.0)
              || (0 < eur_freq && eur_freq < 1.0)
              || (0 < sas_freq && sas_freq < 1.0));
}

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
std::vector<std::tuple<int, char, char>>
get_fixed_afr_sites(const std::string& hg1k_filename, const char* afr_list)
{
    VcfHeader hg1k_header;
    VcfFileReader hg1k_vcf;
    VcfRecord hg1k_rec;

    hg1k_vcf.open(hg1k_filename.c_str(), hg1k_header, afr_list, NULL, NULL);

    std::vector<std::tuple<int, char, char>> fixed_sites;
    while (hg1k_vcf.readRecord(hg1k_rec)) {
        // consider only alleles which
        //   * are biallelic SNPs
        //   * are fixed in all Africans,
        //   * are not fixed outside Africa
        if (has_simple_allele(hg1k_rec)
                && fixed_in_africa(hg1k_rec)
                && polymorphic_outside_africa(hg1k_rec)) {
            int pos = hg1k_rec.get1BasedPosition();

            // get African allele index and index of alternative allele
            // (0 - REF; 1 - ALT)
            short afr_allele_index = hg1k_rec.getGT(0, 0);
            short other_allele_index = afr_allele_index == 0 ? 1 : 0;

            // get alleles of an African genotype and an alternative allele
            char afr_allele = *hg1k_rec.getAlleles(afr_allele_index);
            char other_allele = *hg1k_rec.getAlleles(other_allele_index);

            fixed_sites.push_back(std::make_tuple(pos, afr_allele, other_allele));
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
std::vector<std::tuple<int, char, char>>
check_archaic_states(const std::string& vcf_filename,
                     const short sample_id,
                     const std::vector<std::tuple<int, char, char>>& afr_fixed_sites)
{
    VcfHeader header;
    VcfFileReader vcf;
    VcfRecord rec;
    vcf.open(vcf_filename.c_str(), header);

    // a vector where all valid positions will be accumulated
    std::vector<std::tuple<int, char, char>> result;

    bool skip_reading = false;
    // go through all positions which carry alleles fixed in Africans and look
    // for matching records in a given VCF file
    for (auto& afr_site : afr_fixed_sites) {

        int afr_pos = std::get<0>(afr_site);
        char afr_allele = std::get<1>(afr_site);
        char other_allele = std::get<2>(afr_site);

        if (find_matching_variant(vcf, rec, afr_pos, skip_reading)) {
            // to include this variant for further analysis, archaic
            //   * has to carry a valid allele (biallelic SNP),
            //   * has to be homozygous at this site,
            //   * can't have a missing genotype,
            //   * must be different than all Africans
            //   * must carry the same allele as ALT seen in 1000 genomes data
	        if (((rec.getNumAlts() == 0) || (rec.getNumAlts() == 1))
                    && has_simple_allele(rec)
                    && (rec.getGT(sample_id, 0) == rec.getGT(sample_id, 1))
                    && (rec.getGT(sample_id, 0) != VcfGenotypeSample::MISSING_GT)
                    && (*rec.getAlleles(rec.getGT(sample_id, 0)) != afr_allele)
                    && (*rec.getAlleles(rec.getGT(sample_id, 0)) == other_allele)) {
                // add this position to the final list of sites
                result.push_back(std::make_tuple(
                    afr_pos,
                    afr_allele,
                    *rec.getAlleles(rec.getGT(sample_id, 0)))
                );
            }
        }
    }

    return result;
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
    std::vector<std::tuple<int, char, char>> fixed_sites = get_fixed_afr_sites(hg1k_vcf_file, "tmp/afr_samples.list");
    std::clog << "[Chromosome " << chr << "] Analysis of the 1000 genomes VCF file DONE (" << fixed_sites.size() << " sites)!\n";

    std::vector<std::tuple<int, char, char>> altai, denisovan;

    std::clog << "[Chromosome " << chr << "] Started scanning the Altai VCF file.\n";
    altai = check_archaic_states(altai_vcf_file, 0, fixed_sites);
    std::clog << "[Chromosome " << chr << "] Analysis of the Altai VCF file DONE (" << altai.size() << " sites)!\n";

    std::clog << "[Chromosome " << chr << "] Started scanning the Denisovan VCF file.\n";
    denisovan = check_archaic_states(denisovan_vcf_file, 0, fixed_sites);
    std::clog << "[Chromosome " << chr << "] Analysis of the Denisovan VCF file DONE (" << denisovan.size() << " sites)!\n";

    // initialize the table of final results
    std::map<int, std::tuple<char, char, bool, bool, bool, bool>> table;
    for (auto & site : altai)     table.emplace(std::get<0>(site), std::make_tuple(std::get<1>(site), std::get<2>(site), false, false, false, false));
    for (auto & site : denisovan) table.emplace(std::get<0>(site), std::make_tuple(std::get<1>(site), std::get<2>(site), false, false, false, false));

    // set boolean-flag at positions where an archaic differs from Africans
    for (auto & site : altai)     std::get<2>(table[std::get<0>(site)]) = true;
    for (auto & site : denisovan) std::get<3>(table[std::get<0>(site)]) = true;

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
