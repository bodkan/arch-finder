# Search for putatively introgressed archaic alleles

The program `find_informative_sites` will scan the 1000 Genomes
Project VCF files and look for sites that are informative of
Neanderthal ancestry. Specifically, it will look for sites at
which:

* Yoruba individuals in the 1000 Genomes Project data are fixed for
  a certain allele
* Altai and Vindija high-coverage Neandertal genomes homozygous for an
  allele that is different from the one seen in Yoruba.

Usage: run `make` and follow the instructions.
