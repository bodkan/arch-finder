#include <iostream>
#include <cstdlib>

#include "VcfFileReader.h"

//
// Check the frequency of the major African allele. If it is higher
// than the specified frequency cutoff, return its value. 
//
std::pair<bool, char>
high_freq_in_africa(VcfRecord& hg1k_rec, float arch_freq_cutoff)
{
    int ref_counter = 0;
    int n_africans = hg1k_rec.getNumSamples();
    
    // since only African samples were loaded from the VCF file previously,
    // scan through all samples and count the number of samples carrying
    // the reference allele
    for (int i = 0; i < n_africans; i++) {
        // get indices (0/1) of alleles of the i-th sample
        int allele_1 = hg1k_rec.getGT(i, 0);
        int allele_2 = hg1k_rec.getGT(i, 1);

        if (allele_1 == 0) ref_counter++;
        if (allele_2 == 0) ref_counter++;
    }

    float ref_freq = (float) ref_counter / (2 * n_africans);

    // if the REF/ALT allele is at high frequency, return true
    // and the value of this allele
    if (ref_freq >= (1 - arch_freq_cutoff)) {
        return std::make_pair(true, *hg1k_rec.getAlleles(0));
    } else if ((1 - ref_freq) >= (1 - arch_freq_cutoff)) {
        return std::make_pair(true, *hg1k_rec.getAlleles(1));
    } else {
        return std::make_pair(false, -1);
    }
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
get_high_freq_afr_sites(const std::string& hg1k_filename,
                    const char* afr_list,
                    float freq_cutoff)
{
    VcfHeader hg1k_header;
    VcfFileReader hg1k_vcf;
    VcfRecord hg1k_rec;

    hg1k_vcf.open(hg1k_filename.c_str(), hg1k_header, afr_list, NULL, NULL);

    std::vector<std::pair<int, char>> fixed_sites;
    while (hg1k_vcf.readRecord(hg1k_rec)) {
        bool freq_above_cutoff = false;
        char major_afr_allele;

        std::tie(freq_above_cutoff, major_afr_allele) = high_freq_in_africa(hg1k_rec, freq_cutoff);

        // consider only alleles which
        //   * are biallelic SNPs
        //   * are above specified frequency in Africans
        if (has_simple_allele(hg1k_rec) && freq_above_cutoff) {
            fixed_sites.push_back(std::make_pair(
                hg1k_rec.get1BasedPosition(),
                major_afr_allele)
            );
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
        int pos = afr_site.first;
        char major_afr_allele = afr_site.second;

        if (find_matching_variant(vcf, rec, pos, skip_reading)) {
            // to include this variant for further analysis, archaic
            //   * has to carry a valid allele (biallelic SNP),
            //   * can't have a missing genotype,
            //   * must carry allele different from majority of Africans
	        if (((rec.getNumAlts() == 0) || (rec.getNumAlts() == 1))
                    && has_simple_allele(rec)
                    && (rec.getGT(sample_id, 0) != VcfGenotypeSample::MISSING_GT)
                    && ((*rec.getAlleles(rec.getGT(sample_id, 0)) != major_afr_allele)
                       || (*rec.getAlleles(rec.getGT(sample_id, 1)) != major_afr_allele))) {
                // add this position to the final list of sites
                char allele_1 = *rec.getAlleles(rec.getGT(sample_id, 0));
                char allele_2 = *rec.getAlleles(rec.getGT(sample_id, 1));
                result.emplace(pos, std::make_pair(
                    major_afr_allele,
                    (allele_1 != major_afr_allele) ? allele_1 : allele_2)
                );
            }
        }
    }

    return result;
}

//
// Get a set of sites which at which Altai and Vindija and Africans
// all differ from each other (i.e. triallelic sites that haven't been
// detected earlier by pairwise comparisons).
//
std::set<int>
get_triallelic_positions(std::map<int, std::pair<char, char>> & altai,
                         std::map<int, std::pair<char, char>> & vindija)
{
    std::set<int> triallelic_positions;

    // iterate through all informative sites from Altai...
    for (auto & altai_site : altai) {
      int pos = altai_site.first;
      char altai_allele = altai_site.second.second;
      // ... and check if Vindija (if also carries an informative allele)
      // has the same allele as Altai
      if ((vindija.count(pos) > 0) && (altai_allele != vindija[pos].second))
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
    if (argc != 6) {
        std::cout << "Usage\n\t./find_informative_sites chr_ID 1000G_VCF Altai_VCF Vindija_VCF arch_freq_cutoff\n";
        return 0;
    }

    std::string chr(argv[1]);

    std::string hg1k_vcf_file(argv[2]);
    std::string altai_vcf_file(argv[3]);
    std::string vindija_vcf_file(argv[4]);

    float arch_freq_cutoff = std::atof(argv[5]);

    std::vector<std::pair<int, char>> fixed_sites = get_high_freq_afr_sites(hg1k_vcf_file, "input/afr_samples.list", arch_freq_cutoff);
    std::clog << "[Chromosome " << chr << "] Analysis of the 1000 genomes VCF file done (" << fixed_sites.size() << " sites)\n";

    std::map<int, std::pair<char, char>> altai, vindija;

    altai = check_archaic_states(altai_vcf_file, 0, fixed_sites);
    std::clog << "[Chromosome " << chr << "] Analysis of the Altai VCF file done (" << altai.size() << " sites)\n";

    vindija = check_archaic_states(vindija_vcf_file, 0, fixed_sites);
    std::clog << "[Chromosome " << chr << "] Analysis of the Vindija VCF file done (" << vindija.size() << " sites)\n";

    // get a set of sites at which Altai, Vindija and Africans
    // all differ from each other
    std::set<int> positions_to_remove = get_triallelic_positions(altai, vindija);

    // remove such triallelic sites
    filter_out_triallelic(altai, positions_to_remove);
    filter_out_triallelic(vindija, positions_to_remove);

    // initialize the table of final results
    std::map<int, std::tuple<char, char, bool, bool, bool, bool>> table;
    for (auto & site : altai) table.emplace(site.first, std::make_tuple(site.second.first, site.second.second, false, false, false, false));
    for (auto & site : vindija) table.emplace(site.first, std::make_tuple(site.second.first, site.second.second, false, false, false, false));

    // set boolean-flag at positions where an archaic differs from Africans
    for (auto & site : altai) std::get<2>(table[site.first]) = true;
    for (auto & site : vindija) std::get<3>(table[site.first]) = true;

    // print out the final table in a BED-like format
    for (auto & site : table) {
        std::cout << chr << "\t" << site.first - 1 << "\t" << site.first << "\t"
            << std::get<0>(site.second) << "\t"
            << std::get<1>(site.second) << "\t"
            << std::get<2>(site.second) << "\t"
            << std::get<3>(site.second) << "\n";
    }
    std::clog << "[Chromosome " << chr << "] Printing out the results done " <<
        "(" << table.size() << " sites in total)\n";

    return 0;
}
