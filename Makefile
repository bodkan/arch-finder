# the following paths can be specified on the command line when invoking make
hg1k_vcf_path := /mnt/sequencedb/1000Genomes/ftp/phase3/20140910
altai_vcf_path := /mnt/454/HighCovNeandertalGenome/1_Extended_VCF/AltaiNea
denisovan_vcf_path := /mnt/454/HighCovNeandertalGenome/1_Extended_VCF/DenisovaPinky

chromosomes := $(shell seq 1 22)
arch_freq := 0.0

src_dir := ./src

output_bed_dir := output_bed_arch_freq_$(arch_freq)
output_vcf_dir := output_vcf_arch_freq_$(arch_freq)
bin_dir := bin
lib_dir := lib
tmp_dir := tmp
directories := $(output_vcf_dir) $(output_bed_dir) $(bin_dir) $(lib_dir) $(tmp_dir)

LIBSTATGEN := $(lib_dir)/libStatGen

CXX := g++
CXXFLAGS := -O2 -Wall -Werror -std=c++11

INCLUDES := -I$(LIBSTATGEN)/include
LIBS := -L$(LIBSTATGEN) -lStatGen -lz

hg1k_samples := $(tmp_dir)/1000genomes_samples.table
afr_samples := $(tmp_dir)/afr_samples.list
non_afr_samples := $(tmp_dir)/non_afr_samples.list
eur_samples := $(tmp_dir)/eur_samples.list
eas_samples := $(tmp_dir)/eas_samples.list
sas_samples := $(tmp_dir)/sas_samples.list
amr_samples := $(tmp_dir)/amr_samples.list
all_pops := $(afr_samples) $(yri_samples) $(non_afr_samples) $(eur_samples) $(eas_samples) $(sas_samples) $(amr_samples)

bin := $(bin_dir)/find_informative_sites

informative_sites_per_chr_bed := $(addprefix $(output_bed_dir)/, $(addprefix chr,$(addsuffix .bed,$(chromosomes))))

informative_sites_per_chr_vcf := $(addprefix $(output_vcf_dir)/, $(addprefix chr,$(addsuffix .vcf.gz,$(chromosomes))))

arch_informative_sites_bed := $(output_bed_dir)/arch_informative_sites.bed
nea_informative_sites_bed := $(output_bed_dir)/nea_informative_sites.bed
den_informative_sites_bed := $(output_bed_dir)/den_informative_sites.bed

arch_informative_sites_vcf := $(output_vcf_dir)/arch_informative_sites.vcf.gz
nea_informative_sites_vcf := $(output_vcf_dir)/nea_informative_sites.vcf.gz
den_informative_sites_vcf := $(output_vcf_dir)/den_informative_sites.vcf.gz

arch_informative_sites_tbi := $(output_vcf_dir)/arch_informative_sites.vcf.gz.tbi
nea_informative_sites_tbi := $(output_vcf_dir)/nea_informative_sites.vcf.gz.tbi
den_informative_sites_tbi := $(output_vcf_dir)/den_informative_sites.vcf.gz.tbi

.PHONY: clean scratch

.INTERMEDIATE: $(informative_sites_per_chr_vcf)

default:
	@echo "Usage:"
	@echo "\tmake deps [arch_freq=0.0] -- prepare all depencies (binaries, directories for"
	@echo "\t                             output files, etc.)"
	@echo
	@echo "\tmake scan [arch_freq=0.0] -- scan the genome for archaic-like alleles allowing"
	@echo "\t                             for a certain frequency of such alleles in Africa"
	@echo "\t                             (no archaic-like alleles allowed by default)"
	@echo
	@echo "\tNote: adding the -jN argument to make will run N scans in parallel"
	@echo
	@echo "\tmake clean_deps           -- clean depencies"
	@echo "\tmake clean_results        -- clean results"
	@echo "\tmake clean_all            -- clean everything"
	@echo
	@echo "\tPaths to directories with 1000 genomes VCFs as well as Altai Neanderthal\n \
       and Denisovan VCFs have to be set using hg2k_vcf_path, altai_vcf_path\n \
       and denisovan_vcf_path variables (directly in the Makefile or when invoking\n \
       make)"

deps: $(directories) $(bin)

scan: $(arch_informative_sites_bed) $(nea_informative_sites_bed) $(den_informative_sites_bed) $(arch_informative_sites_vcf) $(arch_informative_sites_tbi) $(nea_informative_sites_vcf) $(nea_informative_sites_tbi) $(den_informative_sites_vcf) $(den_informative_sites_tbi)

$(arch_informative_sites_vcf): $(informative_sites_per_chr_vcf)
	bcftools concat $(addprefix $(output_vcf_dir)/, $(addprefix chr,$(addsuffix .vcf.gz,$(chromosomes)))) --output-type z --output $@

$(nea_informative_sites_vcf): $(arch_informative_sites_vcf) $(nea_informative_sites_bed)
	bcftools view $< -R $(nea_informative_sites_bed) | bgzip > $@

$(den_informative_sites_vcf): $(arch_informative_sites_vcf) $(den_informative_sites_bed)
	bcftools view $< -R $(den_informative_sites_bed) | bgzip > $@

$(output_vcf_dir)/%.vcf.gz.tbi: $(output_vcf_dir)/%.vcf.gz
	tabix -f $<

$(output_vcf_dir)/%.vcf.gz: $(output_bed_dir)/%.bed $(non_afr_samples)
	chr_id=$(subst chr,,$(subst .vcf.gz,,$(notdir $@))); \
	hg1k_vcf_file="$(hg1k_vcf_path)/ALL.chr$${chr_id}.*.vcf.gz"; \
	bcftools norm --multiallelics +snps -R $< $${hg1k_vcf_file} \
	     | bcftools view -M2 -v snps \
	     | bgzip \
	     > $@
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

$(arch_informative_sites_bed): $(arch_informative_sites_vcf)
	cat $(informative_sites_per_chr_bed) > $@

$(nea_informative_sites_bed): $(arch_informative_sites_vcf)
	for i in $(chromosomes); do \
	    awk '$$6 == 1' $(output_bed_dir)/chr$${i}.bed >> $@; \
	done

$(den_informative_sites_bed): $(arch_informative_sites_vcf)
	for i in $(chromosomes); do \
	    awk '$$7 ==1' $(output_bed_dir)/chr$${i}.bed >> $@; \
	done

$(output_bed_dir)/%.bed: $(bin)
	chr_id=$(subst chr,,$(basename $(notdir $@))); \
	hg1k_vcf_file="$(hg1k_vcf_path)/ALL.chr$${chr_id}.*.vcf.gz"; \
	altai_vcf_file="$(altai_vcf_path)/AltaiNea.hg19_1000g.$${chr_id}.mod.vcf.gz"; \
	denisovan_vcf_file="$(denisovan_vcf_path)/DenisovaPinky.hg19_1000g.$${chr_id}.mod.vcf.gz"; \
	$(bin) $${chr_id} $${hg1k_vcf_file} $${altai_vcf_file} $${denisovan_vcf_file} $(arch_freq) > $@

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

$(afr_samples): $(hg1k_samples)
	grep -E "YRI|LWK|GWD|MSL|ESN" $< | cut -f1 > $@

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
	rm -rf $(bin_dir) $(lib_dir) $(tmp_dir)

clean_results:
	rm -rf $(output_bed_dir) $(output_vcf_dir)

clean_all:
	rm -rf $(directories)
