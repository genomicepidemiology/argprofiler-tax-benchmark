import os
import glob

configfile: "config.yaml"

FILES = glob.glob(os.path.join('data/raw', '*_reads.fq.gz'))
SAMPLES = [f.split(os.sep)[-1].replace('_reads.fq.gz', '') for f in FILES]


databases = ["Silva", "mOTUs", "genomic", "metaphlan", "gtdb"]
bracken_levels = "class order family genus species"

rule all:
     input: 
       expand("data/QC/{sample}_reads_fastqc.html", sample=SAMPLES),
       expand("data/Trimmed/{sample}_reads.trim.fq.gz", sample=SAMPLES),
       expand("data/QC/{sample}_reads.trim_fastqc.html", sample=SAMPLES),
       expand("results/metaphlan/{sample}.bowtie2.bz", sample=SAMPLES),
       expand("results/metaphlan/{sample}_profiled.txt", sample=SAMPLES),
       expand("results/motus/mapped_{sample}.sam", sample=SAMPLES),
       expand("results/motus/mgc_{sample}_table.count", sample=SAMPLES),
       expand("results/motus/taxonomy_{sample}_profile.txt", sample=SAMPLES),
       expand("results/kma/{sample}_Silva.fin",sample=SAMPLES),
       expand("results/kma/{sample}_mOTUs.fin",sample=SAMPLES),
       expand("results/kma/{sample}_genomic.fin",sample=SAMPLES),
       expand("results/kma/{sample}_metaphlan.fin",sample=SAMPLES),
       expand("results/kma/{sample}_gtdb.fin",sample=SAMPLES),
       expand("results/kraken2/{sample}.fromreads.kreport", sample=SAMPLES),
       expand("results/kraken2/{sample}.fromreads.kraken", sample=SAMPLES),
       expand("results/bracken/{sample}_{level}.bracken.report", sample=SAMPLES, level=bracken_levels.split()),
       expand("results/bracken/{sample}_{level}.bracken.tab", sample=SAMPLES, level=bracken_levels.split()),
       expand("results/metaphlan/{sample}_profiled_CAMI.txt", sample=SAMPLES),
       expand("results/motus/taxonomy_{sample}_CAMI.txt", sample=SAMPLES),
       expand("results/opal_cami/{sample}_{database}_KMA.cami", sample=SAMPLES, database=databases),
       #expand("results/kma/{sample}_{database}_KMA.cami", sample=SAMPLES, database=databases)

rule qc_pre:
    """ Quality checking of raw reads """
    input:
        "data/raw/{sample}_reads.fq.gz"
    output:
        "data/QC/{sample}_reads_fastqc.html"
    envmodules:
        "tools",
        "perl/5.30.2",
        "jdk/19",
        "fastqc/0.11.9"
    shell:
        """
            /usr/bin/time -v -o {output}.time fastqc {input} -o data/QC
        """

rule trim:
    """ Trimming of raw reads """
    input:
        "data/raw/{sample}_reads.fq.gz"
    output:
        "data/Trimmed/{sample}_reads.trim.fq.gz"
    params:
        qin="auto",
        k="19",
        ref="adapters",
        mink="11",
        qtrim="r",
        trimq="20",
        minlength="50",
        ziplevel="6",
        overwrite="t",
        statscolumns="5",
        ktrim="r"
    envmodules:
        "tools",
        "jdk/19",
        "bbmap/38.90"
    shell:
        """
            /usr/bin/time -v -o {output}.time bbduk.sh in={input} out={output} qin={params.qin} ref={params.ref} k={params.k} qtrim={params.qtrim} mink={params.mink} trimq={params.trimq} minlength={params.minlength} ziplevel={params.ziplevel} overwrite={params.overwrite} statscolumns={params.statscolumns} ktrim={params.ktrim} tbo
        """
    

rule qc_post:
    """ Quality check of trimmed reads """
    input:
        "data/Trimmed/{sample}_reads.trim.fq.gz"
    output:
        "data/QC/{sample}_reads.trim_fastqc.html"
    envmodules:
        "tools",
        "perl/5.30.2",
        "jdk/19",
        "fastqc/0.11.9"
    shell:
        """
            /usr/bin/time -v -o {output}.time fastqc {input} -o data/QC
        """
    
rule metaphlan:
    """ Running MetaPhlan tool on sample with trimmed reads """
    input: 
        fq="data/Trimmed/{sample}_reads.trim.fq.gz"
    output:
        bw="results/metaphlan/{sample}.bowtie2.bz",
        p="results/metaphlan/{sample}_profiled.txt",
    envmodules:
        "tools",
        "bowtie2/2.5.0",
        "metaphlan/4.0.3"
    shell: 
      """
        /usr/bin/time -v -o {output.p}.time metaphlan --bowtie2out {output.bw} --nproc $PBS_NUM_PPN --input_type fastq -o {output.p} --unclassified_estimation {input.fq} 
      """

rule metaphlan_cami:
    """ Converting metahplan output into following CAMI format"""
    input:
        bw="results/metaphlan/{sample}.bowtie2.bz",
    output:
        p="results/metaphlan/{sample}_profiled_CAMI.txt",
    envmodules:
        "tools",
        "bowtie2/2.5.0",
        "metaphlan/4.0.3"
    shell:
      """
        metaphlan {input.bw} --input_type bowtie2out --CAMI_format_output -o {output.p} --nproc $PBS_NUM_PPN
      """
    

rule motus:
    """ Running mOTUs tool on sample with trimmed reads """
    input:
        fq="data/Trimmed/{sample}_reads.trim.fq.gz"
    output:
        sam="results/motus/mapped_{sample}.sam",
        mgc="results/motus/mgc_{sample}_table.count",
        motu="results/motus/taxonomy_{sample}_profile.txt",
    envmodules:
        "tools",
        "motus/3.0.3"
    shell: 
      """
        /usr/bin/time -v -o {output.sam}.time motus map_tax -s {input.fq} -o {output.sam} -t $PBS_NUM_PPN 
        /usr/bin/time -v -o {output.mgc}.time motus calc_mgc -i {output.sam} -o {output.mgc}
        /usr/bin/time -v -o {output.motu}.time motus calc_motu -i {output.mgc} > {output.motu}
      """

rule motus_cami:
    """ Converting mOTUs output to following CAMI format"""
    input:
        mgc="results/motus/mgc_{sample}_table.count",
    output:
        cmotu="results/motus/taxonomy_{sample}_CAMI.txt"
    envmodules:
        "tools",
        "motus/3.0.3"
    shell: 
     """
        motus calc_motu -i {input.mgc} -C precision > {output.cmotu}
     """

rule kma_silva:
    """ Running KMA tool on sample with trimmed reads using the Silva database"""
    input:
        fq="data/Trimmed/{sample}_reads.trim.fq.gz"
    params:
        db="/home/databases/metagenomics/db/Silva_20200116/Silva_20200116"
    envmodules:
        "tools",
        "kma/1.4.7"
    output:
        kma_finished = "results/kma/{sample}_Silva.fin",
    shell: 
      """
        /usr/bin/time -v -o results/kma/{wildcards.sample}_Silva.time kma -i {input.fq} -nc -na -nf -o results/kma/{wildcards.sample}_Silva -t_db {params.db} -mem_mode -ef -1t1 -apm p -oa

        ./check_status.sh results/kma/{wildcards.sample}_Silva.time {output.kma_finished}
      """

rule kma_mOTUs:
    """ Running KMA tool on sample with trimmed reads using the mOTUs database"""
    input:
        fq="data/Trimmed/{sample}_reads.trim.fq.gz"
    params:
        db="/home/databases/metagenomics/db/mOTUs_20221205/db_mOTU_20221205"
    envmodules:
        "tools",
        "kma/1.4.7"
    output:
        kma_finished = "results/kma/{sample}_mOTUs.fin",

    shell: 
      """
        /usr/bin/time -v -o results/kma/{wildcards.sample}_mOTUs.time kma -i {input.fq} -nc -na -nf -o results/kma/{wildcards.sample}_mOTUs -t_db {params.db} -mem_mode -ef -1t1 -apm p -oa

        ./check_status.sh results/kma/{wildcards.sample}_mOTUs.time {output.kma_finished}
      """

rule kma_genomic:
    """ Running KMA tool on sample with trimmed reads using the genomic database"""
    input:
        fq="data/Trimmed/{sample}_reads.trim.fq.gz"
    params:
        db="/home/databases/metagenomics/kma_db/genomic_20220524/genomic_20220524"
    envmodules:
        "tools",
        "kma/1.4.7"
    output:
        kma_finished = "results/kma/{sample}_genomic.fin",

    shell: 
      """
        /usr/bin/time -v -o results/kma/{wildcards.sample}_genomic.time kma -i {input.fq} -nc -na -nf -o results/kma/{wildcards.sample}_genomic -t_db {params.db} -mem_mode -ef -1t1 -apm f -oa

        ./check_status.sh results/kma/{wildcards.sample}_genomic.time {output.kma_finished}
      """


rule kma_metaphlan:
    """ Running KMA tool on sample with trimmed reads using the MetaPhlan database"""
    input:
        fq="data/Trimmed/{sample}_reads.trim.fq.gz"
    params:
        db="/home/databases/metagenomics/db/metaphlan_20221125/metaphlan_20221125"
    envmodules:
        "tools",
        "kma/1.4.7"
    output:
        kma_finished = "results/kma/{sample}_metaphlan.fin",
    shell: 
      """
        /usr/bin/time -v -o results/kma/{wildcards.sample}_metaphlan.time kma -i {input.fq} -nc -na -nf -o results/kma/{wildcards.sample}_metaphlan -t_db {params.db} -mem_mode -ef -1t1 -apm p -oa

        ./check_status.sh results/kma/{wildcards.sample}_metaphlan.time {output.kma_finished}
      """

rule kma_gtdb:
    """ Running KMA tool on sample with trimmed reads using the GTDB database"""
    input:
        fq="data/Trimmed/{sample}_reads.trim.fq.gz"
    params:
        db="/home/databases/metagenomics/db/gtdb_bac120_20220726/bac120_marker_genes_r207_k16_20220726"
    envmodules:
        "tools",
        "kma/1.4.7"
    output:
        kma_finished = "results/kma/{sample}_gtdb.fin",
    shell: 
      """
        /usr/bin/time -v -o results/kma/{wildcards.sample}_gtdb.time kma -i {input.fq} -nc -na -nf -o results/kma/{wildcards.sample}_gtdb -t_db {params.db} -mem_mode -ef -1t1 -apm p -oa
        
        ./check_status.sh results/kma/{wildcards.sample}_gtdb.time {output.kma_finished} 
      """

#rule kma2cami:
#    input:
#        kma_fin ="results/kma/{sample}_{database}.fin",
#    output:
#        kma_cami = "results/kma/{sample}_{database}_KMA.cami"
#    envmodules:
#        "tools",
#        "anaconda3/2022.10",
#        "mariadb/10.4.12",
#        "mariadb-connector-c/3.1.7"
#    params:
#        kma_mapstat ="results/kma/{sample}_{database}.mapstat",
#    shell:
#     """
#        python to_cami.py -f kma -p {params.kma_mapstat} -o results/kma 
#     """

rule kraken:
    input:
        fq="data/Trimmed/{sample}_reads.trim.fq.gz"
    output:
        kraken_report="results/kraken2/{sample}.fromreads.kreport",
        kraken_out="results/kraken2/{sample}.fromreads.kraken"
    params:
        db="/home/databases/kraken2/"
    envmodules:
        "tools",
        "kraken/2.1.2"
    shell:
        """
            kraken2 --threads $PBS_NUM_PPN --db {params.db} --output {output.kraken_out} --report {output.kraken_report} {input.fq}
        """

rule bracken:
    input:
        kraken_report="results/kraken2/{sample}.fromreads.kreport",
    output:
        br1=expand("results/bracken/{{sample}}_{level}.bracken.report", level=bracken_levels.split()),
        br2=expand("results/bracken/{{sample}}_{level}.bracken.tab", level=bracken_levels.split())
    envmodules:
        "tools",
        "bracken/2.8"
    params:
        db="/home/databases/kraken2/bracken/",
        threshold="0",
        readlength="100",
        levels=bracken_levels
    shell:
        """
        #declare -A levels=( ["species"]="S" ["genus"]="G" ["family"]="F" ["order"]="O" ["class"]="C")
        
        # species S
        bracken -d {params.db} -i {input.kraken_report} -o results/bracken/{wildcards.sample}_species.bracken.report -w results/bracken/{wildcards.sample}_species.bracken.tab -l {params.readlength} -l S -t {params.threshold}

        # genus G
        bracken -d {params.db} -i {input.kraken_report} -o results/bracken/{wildcards.sample}_genus.bracken.report -w results/bracken/{wildcards.sample}_genus.bracken.tab -l {params.readlength} -l G -t {params.threshold}

        # family F
        bracken -d {params.db} -i {input.kraken_report} -o results/bracken/{wildcards.sample}_family.bracken.report -w results/bracken/{wildcards.sample}_family.bracken.tab -l {params.readlength} -l F -t {params.threshold}

        # order O
        bracken -d {params.db} -i {input.kraken_report} -o results/bracken/{wildcards.sample}_order.bracken.report -w results/bracken/{wildcards.sample}_order.bracken.tab -l {params.readlength} -l O -t {params.threshold}

        # class C
        bracken -d {params.db} -i {input.kraken_report} -o results/bracken/{wildcards.sample}_class.bracken.report -w results/bracken/{wildcards.sample}_class.bracken.tab -l {params.readlength} -l C -t {params.threshold}
        """


rule bracken2cami:
    input:
        br2="results/bracken/{sample}_{level}.bracken.tab"
    output:
        cami="results/bracken/{sample}_{level}_bracken.cami"
    envmodules:
        "tools",
        "anaconda3/2022.10",
        "mariadb/10.4.12",
        "mariadb-connector-c/3.1.7"
    shell:
     """
        python to_cami.py -f bracken -p {input.br2} -o results/bracken 
     """

rule kma2cami:
    input:
       mapstat="results/kma/{sample}_{database}.fin",
    output:
       "results/opal_cami/{sample}_{database}_KMA.cami",
    envmodules:
        "tools",
        "anaconda3/2022.10",
        "mariadb/10.4.12",
        "mariadb-connector-c/3.1.7"
    shell:
     """
        python kma2cami.py -o results/opal_cami/ -v -p results/kma/{wildcards.sample}_{wildcards.database}.mapstat 
     """
