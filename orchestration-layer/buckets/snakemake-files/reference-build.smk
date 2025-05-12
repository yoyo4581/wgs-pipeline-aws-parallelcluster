from globals import data_repo

rule all:
    input:
        "%s/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.ann" % data_repo


print(f"{data_repo}")
rule get_genome_ref:
    output:
        "%s/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa" % data_repo
    shell:
        """
        mkdir -p {0}/reference/GRCh38
        wget -P {0}/reference/GRCh38 ftp://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        """.format(data_repo)

rule unzip:
    input:
        "%s/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" % data_repo
    output:
        "%s/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa" % data_repo
    shell:
        """
        gunzip -c {0}/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > {0}/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa"
        """.format(data_repo)

rule index_genome:
    input:
        "%s/reference/GRCh38/{genome}.fa" % data_repo
    output:
        idx=multiext("%s/reference/GRCh38/{genome}.fa" % data_repo, ".ann"),
    log:
        "logs/bwa_index/{genome}.fa.log"
    wrapper:
        "v3.13.2/bio/bwa/index"

