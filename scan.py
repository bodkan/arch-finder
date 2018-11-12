import os

import pandas as pd

df = pd.read_csv("tmp/sample_info.csv")
africans = " ".join(df.query("superpop == 'EUR'")["sample"])

for chrom in range(1, 23):
    print(f"Processing chromosome {chrom}...", end="")
    os.system(
        "bin/arch_finder "
        f"--vcf /mnt/sequencedb/gendivdata/2_genotypes/giantVcfs/merged_var_nosing_sites_arch_apes_sgdp1_g1000_chr{chrom}.vcf.gz "
        "--archaics AltaiNeandertal Vindija33.19 Chagyrskaya-Phalanx "
        f"--africans {africans} >> yri_sites.bed"
    )
    print("done.")
