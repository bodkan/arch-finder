AFR := YRI

dirs := input output tmp

sample_info := input/sample_info.csv
info_sites := output/info_sites_$(AFR).tsv

default:
	@echo "Usage:"
	@echo "\tmake scan AFR=<1000 Genomes Project African population to use as a reference - YRI by default>"

scan: $(info_sites)

$(info_sites): $(sample_info)
	./src/arch_finder.R $(AFR) 0.0 $@
	# perl -i -lane 'print $$F[0] . "\t" . ($$F[1] - 1) . "\t" . $$F[1]' $@


$(sample_info): $(dirs)
	cp /mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/sample_info.csv $@

clean:
	rm -rf $(dirs)

$(dirs):
	mkdir -p $@
