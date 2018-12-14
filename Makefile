AFR := YRI

dirs := bin input output tmp

arch_finder := ./bin/arch_finder
sample_info := input/sample_info.csv
info_sites := output/info_sites_$(AFR).bed

default:
	@echo "Usage:"
	@echo "\tmake scan AFR=<1000 Genomes Project African population to use as a reference - YRI by default>"
	@echo
	@echo "Note: adding the -jN will run N scans in parallel, one for each chromosome.\n"

scan: $(info_sites)

$(info_sites): $(addsuffix .bed,$(addprefix tmp/, $(addprefix info_sites_chr,$(shell seq 1 22))))
	cat $^ > $@

tmp/info_sites_chr%.bed: /mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr%.vcf.gz $(dirs) $(sample_info) $(arch_finder)
	$(arch_finder) \
        --vcf $< \
        --bed $@ \
        --outgroups panTro4 panPan1.1 gorGor3 \
        --archaics AltaiNeandertal Vindija33.19 Chagyrskaya-Phalanx \
        --africans $(shell perl -F, -lane 'printf "$$F[0] " if ($$F[1] =~ "$(AFR)")' $(sample_info))

$(arch_finder): $(dirs)
	g++ -std=c++14 -Wall src/arch_finder.cpp -lboost_iostreams -o $@

$(sample_info): $(dirs)
	cp /mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/sample_info.csv $@

clean:
	rm -rf $(dirs)

$(dirs):
	mkdir -p $@
