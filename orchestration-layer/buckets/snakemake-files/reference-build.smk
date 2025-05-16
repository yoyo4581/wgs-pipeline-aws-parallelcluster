from globals import data_repo

rule all:
    input:
        "%s/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa" % data_repo,
        "%s/reference/supportingFiles/Homo_sapiens_assembly38.dbsnp138.vcf" % data_repo,
        "%s/reference/supportingFiles/mutect2_files/1000g_pon.hg38.vcf.gz" % data_repo,
        "%s/reference/supportingFiles/mutect2_files/af-only-gnomad.hg38.vcf.gz" % data_repo,
        "%s/reference/supportingFiles/mutect2_files/exome_calling_regions.v1.1.interval_list" % data_repo



print(f"{data_repo}")
rule get_genome_ref:
    output:
        "%s/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" % data_repo
    shell:
        """
        mkdir -p {0}/reference/GRCh38
        wget -P {0}/reference/GRCh38 ftp://ftp.ensembl.org/pub/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
        """.format(data_repo)


rule extract_genome:
    input:
        "%s/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz" % data_repo
    output:
        "%s/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa" % data_repo
    shell:
        """
        gunzip -c {0}/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > {0}/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
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

rule germline_snp_indels:
    output:
        "%s/reference/supportingFiles/Homo_sapiens_assembly38.dbsnp138.vcf" % data_repo
    shell:
        "gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dbsnp138.vcf %s/reference/supportingFiles/" % data_repo

rule FuncotatorSomaticDataSource:
    output:
        directory("%s/reference/supportingFiles/funcotator_dataSources.v1.8.hg38.20230908s" % data_repo)
    shell:
        """
        gatk FuncotatorDataSourceDownloader --somatic --validate-integrity --extract-after-download --hg38
        """

rule germline_allele_freq:
    output:
        "%s/reference/supportingFiles/mutect2_files/af-only-gnomad.hg38.vcf.gz" % data_repo
    shell:
        "gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz %s/reference/supportingFiles/mutect2_files/" % data_repo

rule pon_retrieve:
    output:
        "%s/reference/supportingFiles/mutect2_files/1000g_pon.hg38.vcf.gz" % data_repo
    shell:
        "gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz %s/reference/supportingFiles/mutect2_files/" % data_repo

rule exome_intervals:
    output:
        "%s/reference/supportingFiles/mutect2_files/exome_calling_regions.v1.1.interval_list" % data_repo
    shell:
        "gsutil cp gs://gcp-public-data--broad-references/hg38/v0/exome_calling_regions.v1.1.interval_list %s/reference/supportingFiles/mutect2_files/" % data_repo
        