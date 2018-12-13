dirs := bin input output tmp

arch_finder := ./bin/arch_finder
sample_info := input/sample_info.tsv



output/info_sites_chr%.bed: /mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr%.vcf.gz $(dirs) $(sample_info) $(arch_finder)
	africans=`perl -F, -lane 'printf "$$F[0] " if ($$F[1] == "YRI")' $(sample_info)`; \
	$(arch_finder) \
        --vcf $< \
        --bed $@ \
        --outgroups panTro4 panPan1.1 gorGor3 \
        --archaics AltaiNeandertal Vindija33.19 Chagyrskaya-Phalanx \
        --africans $$africans

$(arch_finder): $(dirs)
	g++ -std=c++14 -Wall src/arch_finder.cpp -lboost_iostreams -o $@

$(sample_info): $(dirs)
	cp /mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/sample_info.csv $@

clean:
	rm -rf $(dirs)

$(dirs):
	mkdir -p $@
