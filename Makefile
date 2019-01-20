AFR := YRI

dirs := input output tmp

sample_info := input/sample_info.csv
afr_fixed := output/$(AFR)_fixed.tsv
neand_derived := output/neand_derived.tsv

default:
	@echo "Usage:"
	@echo "\tmake afr_fixed AFR=<1000 Genomes Project African population to use as a reference - YRI by default>"
	@echo "\tmake neand_derived"

afr_fixed: $(afr_fixed)
neand_derived: $(neand_derived)

$(afr_fixed): $(sample_info)
	./src/arch_finder.R $(AFR) 0.0 $@

$(neand_derived): $(sample_info)
	./src/neand_derived.R $@

$(sample_info): $(dirs)
	cp /mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/sample_info.csv $@

clean:
	rm -rf $(dirs)

$(dirs):
	mkdir -p $@
