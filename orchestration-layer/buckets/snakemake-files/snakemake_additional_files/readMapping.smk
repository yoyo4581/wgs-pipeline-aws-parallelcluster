from globals import data_repo, cell_line_sra_pairs

# Rule to generate FastQC reports on each forward/rev of each run
rule preQC_report:
    input:
        forw=f"{data_repo}/temp/{{cell_line}}/{{SRA}}/{{SRA}}_1.fastq.gz",
        reve=f"{data_repo}/temp/{{cell_line}}/{{SRA}}/{{SRA}}_2.fastq.gz"
    output:
        directory(f"{data_repo}/output/{{cell_line}}/{{SRA}}/prefastQC")
    shell:
        """
        mkdir -p {output} && \
        fastqc {input.forw} {input.reve} -o {output}
        """

# Rule to generate MultiQC report of pre-trimmed reads in a results directory organized by cell line
rule preMultiQC:
    input:
        expand(
            f"{data_repo}/output/{{cell_line}}/{{SRA}}/prefastQC/",
            zip,
            cell_line=[pair[0] for pair in cell_line_sra_pairs],
            SRA = [ pair[1] for pair in cell_line_sra_pairs ]
        )
    output:
        f"{data_repo}/output/{{cell_line}}/preMultiqc_report.html"
    shell:
        f"""
        multiqc {{input}} -o {data_repo}/output/{{wildcards.cell_line}}/ --filename preMultiqc_report.html
        """

# Rule to trim paired reads using Trimmomatic
rule trim:
    input:
        forw=f"{data_repo}/temp/{{cell_line}}/{{SRA}}/{{SRA}}_1.fastq.gz",
        reve=f"{data_repo}/temp/{{cell_line}}/{{SRA}}/{{SRA}}_2.fastq.gz"
    output:
        protected(f"{data_repo}/input/{{cell_line}}/{{SRA}}_1_trimmed_paired.fastq.gz"),
        protected(f"{data_repo}/input/{{cell_line}}/{{SRA}}_2_trimmed_paired.fastq.gz"),
        temp(f"{data_repo}/input/{{cell_line}}/{{SRA}}_1_trimmed_unpaired.fastq.gz"),
        temp(f"{data_repo}/input/{{cell_line}}/{{SRA}}_2_trimmed_unpaired.fastq.gz")
    params:
        adapter_file=f"{data_repo}/reference/adapters/TruSeq3-PE.fa"
    shell:
        """
        trimmomatic PE -phred33 \
        {input.forw} {input.reve} \
        {output[0]} {output[2]} \
        {output[1]} {output[3]} \
        ILLUMINACLIP:{params.adapter_file}:2:30:10 \
        SLIDINGWINDOW:4:30 \
        MINLEN:75
        """


# Rule to generate FastQC reports on each forward/rev of each run from trimmed reads
rule postQC_report:
    input:
        forw=f"{data_repo}/input/{{cell_line}}/{{SRA}}_1_trimmed_paired.fastq.gz",
        reve=f"{data_repo}/input/{{cell_line}}/{{SRA}}_2_trimmed_paired.fastq.gz"
    output:
        directory(f"{data_repo}/output/{{cell_line}}/{{SRA}}/postfastQC")
    shell:
        """
        mkdir -p {output} && \
        fastqc {input.forw} {input.reve} -o {output}
        """


# Rule to generate MultiQC report in a results directory organized by cell line
rule postMultiQC:
    input:
        expand(
            f"{data_repo}/output/{{cell_line}}/{{SRA}}/postfastQC/",
            zip,
            cell_line=[pair[0] for pair in cell_line_sra_pairs],
            SRA = [ pair[1] for pair in cell_line_sra_pairs ]
        )
    output:
        f"{data_repo}/output/{{cell_line}}/postMultiqc_report.html"
    shell:
        f"""
        multiqc {{input}} -o {data_repo}/output/{{wildcards.cell_line}}/ --filename postMultiqc_report.html
        """


rule bwa_mem:
    input:
        reads=[f"{data_repo}/input/{{cell_line}}/{{SRA}}_1_trimmed_paired.fastq.gz", 
        f"{data_repo}/input/{{cell_line}}/{{SRA}}_2_trimmed_paired.fastq.gz"],
        idx=f"{data_repo}/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa.ann",
    output:
        protected(f"{data_repo}/input/{{cell_line}}/sorted_{{SRA}}.bam")
    log:
        "logs/bwa_mem/{cell_line}/{SRA}.log",
    params:
        extra=lambda wildcards: f"-R '@RG\\tID:{wildcards.SRA}\\tPL:ILLUMINA\\tLB:{wildcards.cell_line}_lib1\\tSM:{wildcards.cell_line}\\tPI:300'",
        sorting="samtools",
        sort_order="coordinate",
        sort_extra="",
    threads: 8
    wrapper:
        "v5.0.0/bio/bwa/mem"


rule deleteSamples:
    input:
        lambda wildcards: expand(
            "mapped_reads/{cell_line}/sorted_{SRA}.bam",
            cell_line=wildcards.cell_line,
            SRA=config["metadata"][wildcards.cell_line]
        )
    output:
        "data/deletedSamples.txt"
    shell:
        """
        rm -r data/samples
        echo "Samples Deleted after processing" > {output}
        """


rule samtools_index:
    input:
        "mapped_reads/{cell_line}/sorted_{SRA}.bam"
    output:
        "mapped_reads/{cell_line}/sorted_{SRA}.bam.bai"
    shell:
        "samtools index {input}"