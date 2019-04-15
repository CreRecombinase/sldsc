
all_tissues=["PBMC","Adult_Liver","Adipose_Nuclei"]
all_annot=expand("ld/{tiss}_H3K27ac.",tiss=all_tissues)
pref=",".join(all_annot)
rule all:
    input:
        expand("ld/{tiss}_H3K27ac.{chrom}.l2.ldscore.gz",tiss=["Adult_Liver","Adipose_Nuclei","PBMC"],chrom=range(1,23))


rule check_ldsc:
    output:
        temp("installed_ldsc.txt")
    conda:
        "ldsc/environment.yml"
    shell:
        "touch {output}"

rule clean_bed:
    input:
        bedf="bed/{tiss}.H3K27ac.bed",
    output:
        outf="anno_bed/{tiss}.H3K27ac.bed",
    shell:
        "cut -f 1,2,3 {input} > {output}"


rule annot:
    ''' Make ldsc-friendly annotation files'''
    input:
        bedf="anno_bed/{tiss}.H3K27ac.bed",
        bimf="1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}.bim"
    output:
        annof="ld/{tiss}_H3K27ac.{chrom}.annot.gz"
    conda:
        "ldsc/environment.yml"
    shell:
        "python ldsc/make_annot.py --bed-file {input.bedf} --bimfile {input.bimf} --annot-file {output.annof}"

rule ldsc_est:
    input:
        bimf="1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}.bim",
        bedf="1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}.bed",
        famf="1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}.fam",
        annof="ld/{tiss}_H3K27ac.{chrom}.annot.gz",
        hmsnplist="hapmap3_snps/hm.{chrom}.snp"
    params:
        bfile="1000G_EUR_Phase3_plink/1000G.EUR.QC.{chrom}",
        opref="ld/{tiss}_H3K27ac.{chrom}"
    output:
        ldscore="ld/{tiss}_H3K27ac.{chrom}.l2.ldscore.gz",
    conda:
        "ldsc/environment.yml"
    shell:
        "python2 ldsc/ldsc.py --l2 --bfile {params.bfile} --ld-wind-cm 1 --annot {input.annof} --thin-annot --out {params.opref} --print-snps {input.hmsnplist}"

rule munge_sumstats:
    input:
        gwasf="GWAS/{trait}_GWAS.txt",
        hapmap_list="w_hm3.snplist.txt"
    params:
        pref="mg/{trait}"
    output:
        "mg/{trait}.sumstats.gz"
    conda:
        "ldsc/environment.yml"
    shell:
        "python2 ldsc/munge_sumstats.py --sumstats {input.gwasf} --merge-alleles {input.hapmap_list}  --out {params.pref} --a1-inc"

rule part_h2:
    input:
        mungef="mg/{trait}.sumstats.gz",
        weightf=expand("ld/{tiss}_H3K27ac.{chrom}.l2.ldscore.gz",tiss=all_tissues,chrom=range(1,23)),
        basef=expand("weights_hm3_no_hla/weights.{chrom}.l2.ldscore.gz",chrom=range(1,23)),
        frqf=expand("1000G_Phase3_frq/1000G.EUR.QC.{chrom}.frq",chrom=range(1,23)),
        annof=expand("ld/{tissue}_H3K27ac.{chrom}.annot.gz",tissue=all_tissues,chrom=range(1,23))
    output:
        outputf="{trait}.results",
        logf="{trait}.log"
    params:
        pref="{trait}"
    conda:
        "ldsc/environment.yml"
    shell:
        "python ldsc/ldsc.py --h2 {input.mungef} --ref-ld-chr "+pref+" --overlap-annot --frqfile-chr 1000G_Phase3_frq/1000G.EUR.QC. --w-ld-chr weights_hm3_no_hla/weights. --out {params.pref}"
