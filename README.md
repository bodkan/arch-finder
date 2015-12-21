# Search for putatively introgressed archaic alleles

The program `find_informative_sites` will scan the 1000 Genomes
Project VCF files and look for sites that are informative of
Neanderthal or Denisovan ancestry. Specifically, it will look for
sites where:

* Africans in 1000 Genomes Project data are fixed for a certain allele
* Altai Neanderthal or Denisovan individuals are homozygous for an
  allele that is different from the one seen in Africans
