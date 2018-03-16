SHELL := /bin/bash

# the following paths can be specified on the command line when invoking make
hg1k_vcf_path := /mnt/sequencedb/1000Genomes/ftp/phase3/20140910
altai_vcf_path := /mnt/454/Vindija/high_cov/genotypes/Altai
vindija_vcf_path := /mnt/454/Vindija/high_cov/genotypes/Vindija33.19

chromosomes := $(shell seq 1 22)
nea_freq := 0.0

src_dir := ./src

input_dir := input
output_bed_dir := output_bed
output_vcf_dir := output_vcf
bin_dir := bin
lib_dir := lib
tmp_dir := tmp
directories := $(bin_dir) $(lib_dir) $(tmp_dir) $(input_dir) $(output_bed_dir) $(output_vcf_dir)

LIBSTATGEN := $(lib_dir)/libStatGen

CXX := g++
CXXFLAGS := -O2 -Wall -Werror -std=c++11

INCLUDES := -I$(LIBSTATGEN)/include
LIBS := -L$(LIBSTATGEN) -lStatGen -lz

hg1k_samples := $(input_dir)/1000genomes_samples.table
all_samples := $(input_dir)/all_samples.list
afr_samples := $(input_dir)/afr_samples.list
non_afr_samples := $(input_dir)/non_afr_samples.list
eur_samples := $(input_dir)/eur_samples.list
eas_samples := $(input_dir)/eas_samples.list
sas_samples := $(input_dir)/sas_samples.list
amr_samples := $(input_dir)/amr_samples.list
all_pops := $(all_samples) $(afr_samples) $(yri_samples) $(non_afr_samples) $(eur_samples) $(eas_samples) $(sas_samples) $(amr_samples)

bin := $(bin_dir)/find_informative_sites

informative_sites_per_chr_bed := $(addprefix $(tmp_dir)/, $(addprefix chr,$(addsuffix _nea_freq_$(nea_freq).bed,$(chromosomes))))

informative_sites_per_chr_vcf := $(addprefix $(output_vcf_dir)/, $(addprefix chr,$(addsuffix _nea_freq_$(nea_freq).vcf.gz,$(chromosomes))))
informative_sites_per_chr_tbi := $(addprefix $(output_vcf_dir)/, $(addprefix chr,$(addsuffix _nea_freq_$(nea_freq).vcf.gz.tbi,$(chromosomes))))

informative_sites_bed := $(output_bed_dir)/informative_sites_nea_freq_$(nea_freq).bed

informative_sites_vcf := $(output_vcf_dir)/informative_sites_nea_freq_$(nea_freq).vcf.gz
informative_sites_tbi := $(output_vcf_dir)/informative_sites_nea_freq_$(nea_freq).vcf.gz.tbi

.PHONY: clean scratch

#.INTERMEDIATE: $(informative_sites_per_chr_bed)

default:
	@echo "Usage:"
	@echo -e "\tmake deps                 -- prepare depencies (binaries, directories etc.)"
	@echo
	@echo -e "\tmake scan [nea_freq=0.0]  -- scan the genome for archaic-like alleles allowing"
	@echo -e "\t                             for a certain frequency of such alleles in Africa"
	@echo -e "\t                             (no archaic-like alleles allowed by default)"
	@echo
	@echo -e "\tNote: adding the -jN argument to make will run N scans in parallel"
	@echo
	@echo -e "\tmake clean_deps           -- clean depencies"
	@echo -e "\tmake clean_results        -- clean results"
	@echo -e "\tmake clean_all            -- clean everything"
	@echo
	@echo -e "\tPaths to directories with 1000 genomes VCFs as well as Altai Neandertal\n \
       and Vindija Neandertal VCFs have to be set using hg1k_vcf_path, altai_vcf_path\n \
       and vinija_vcf_path variables (directly in the Makefile or when invoking\n \
       make on the command-line)"

deps: $(directories) $(bin)

scan: $(output_vcf_dir) $(output_bed_dir) $(informative_sites_bed) $(informative_sites_per_chr_vcf) $(informative_sites_per_chr_tbi) $(informative_sites_vcf) $(informative_sites_tbi)

$(informative_sites_vcf): $(informative_sites_per_chr_vcf)
	bcftools concat $(informative_sites_per_chr_vcf) --output-type z --output $@

$(output_vcf_dir)/%.vcf.gz.tbi: $(output_vcf_dir)/%.vcf.gz
	tabix -f $<

$(output_vcf_dir)/chr%.vcf.gz: $(tmp_dir)/chr%.bed
	chr_id=$(subst chr,,$(subst _nea_freq_$(nea_freq).vcf.gz,,$(notdir $@))); \
	hg1k_vcf_file="$(hg1k_vcf_path)/ALL.chr$${chr_id}.*.vcf.gz"; \
	bcftools norm --multiallelics +snps -R $< $${hg1k_vcf_file} \
	     | bcftools view -M2 -v snps \
	     | bgzip \
	     > $@; \
	# the following hack solves the problem of multiallelic sites
	# being sometimes represented by multiple monoallelic records
	# in the 1000 genomes VCF files -- bcftools norm command above
	# gets rid of these while generating the subsets; however,
	# find_informative_sites.cpp that generates BED files of
	# putatively introgressed sites works on unnormalized VCFs --
	# the following narrows the BED files down to only those sites
	# that are truly biallelic
	mv $< $<_tmp; bedtools intersect -a $<_tmp -b $@ -sorted > $<; rm $<_tmp
	touch $@

$(informative_sites_bed): $(informative_sites_per_chr_bed) $(informative_sites_vcf)
	cat $(informative_sites_per_chr_bed) > $@_tmp; \
	# add a column with a reference allele to the final BED file:
	paste $@_tmp <(bcftools view -H $(informative_sites_vcf) | cut -f4) \
		| awk -v OFS="\t" '{ print $$1, $$2, $$3, $$8, $$4, $$5, $$6, $$7 }' \
		> $@; \
	rm $@_tmp

$(tmp_dir)/chr%.bed: $(bin)
	chr_id=$(subst chr,,$(subst _nea_freq_$(nea_freq),,$(basename $(notdir $@)))); \
	hg1k_vcf_file="$(hg1k_vcf_path)/ALL.chr$${chr_id}.*.vcf.gz"; \
	altai_vcf_file="$(altai_vcf_path)/chr$${chr_id}_mq25_mapab100.vcf.gz"; \
	vindija_vcf_file="$(vindija_vcf_path)/chr$${chr_id}_mq25_mapab100.vcf.gz"; \
	$(bin) $${chr_id} $${hg1k_vcf_file} $${altai_vcf_file} $${vindija_vcf_file} $(nea_freq) > $@

$(bin): src/find_informative_sites.cpp $(all_pops) $(LIBSTATGEN)
	$(CXX) $< -o $@ $(INCLUDES) $(LIBS) $(CXXFLAGS)

scratch:
	$(CXX) src/scratch.cpp -o $(bin_dir)/scratch \
	    $(INCLUDES) $(LIBS) $(CXXFLAGS)

$(LIBSTATGEN):
	git clone --recursive https://github.com/statgen/libStatGen.git $@
	cd $@; make

$(hg1k_samples):
	curl ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel \
	    | tail -n+2 > $@

$(all_samples): $(hg1k_samples)
	cut -f1 $< > $@

$(afr_samples): $(hg1k_samples)
	grep -E "YRI" $< | cut -f1 > $@
#	grep -E "YRI|LWK|GWD|MSL|ESN" $< | cut -f1 > $@

$(non_afr_samples): $(hg1k_samples)
	grep -v "AFR" $< | cut -f1 > $@

$(eur_samples): $(hg1k_samples)
	grep "EUR" $< | cut -f1 > $@

$(eas_samples): $(hg1k_samples)
	grep "EAS" $< | cut -f1 > $@

$(sas_samples): $(hg1k_samples)
	grep "SAS" $< | cut -f1 > $@

$(amr_samples): $(hg1k_samples)
	grep "AMR" $< | cut -f1 > $@

$(directories):
	mkdir -p $@

clean_deps:
	rm -rf $(bin_dir) $(lib_dir) $(tmp_dir) $(input_dir)

clean_results:
	rm -rf $(output_bed_dir) $(output_vcf_dir) $(tmp_dir)

clean_all:
	rm -rf $(directories)
