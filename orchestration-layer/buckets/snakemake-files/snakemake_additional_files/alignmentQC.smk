#rule to get basic alignment statistics, such as total reads, mapped reads, and duplicates.
rule align_flagstat:
    input:
        bam=f"{data_repo}/input/{{cell_line}}/sorted_{{SRA}}.bam"
    output:
        f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/flagstat_output.txt"
    shell:
        "samtools flagstat {input.bam} > {output}"

#more detailed, read length distribution, insert sizes, and mapping quality scores.
rule align_stats:
    input:
        bam = f"{data_repo}/input/{{cell_line}}/sorted_{{SRA}}.bam"
    output:
        f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/alignment_stats.txt"
    threads: 4
    shell:
        "samtools stats {input.bam} > {output}"

#assess coverage across the genome, especially important in resequencing or variant calling.
rule align_genomeCoverage:
    input:
        bam=f"{data_repo}/input/{{cell_line}}/{{cell_line}}_merged.bam",
    output:
        f"{data_repo}/output/{{cell_line}}/tumor_coverage.bedgraph"
    threads: 4
    log:
        f"{data_repo}/logs/{{cell_line}}/tumor_bedCoverage.log"
    shell:
        """
        bedtools genomecov -ibam {input.bam} -bga > {output}
        """

#comprehensive report including GC distribution, coverage uniformity, insert size distribution
rule align_qualimap:
    input:
        bam=f"{data_repo}/input/{{cell_line}}/sorted_{{SRA}}.bam"
    output:
        directory(f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/qualimap_report/")
    threads: 4
    shell:
        "qualimap bamqc -bam {input.bam} -outdir {output} -c --java-mem-size=24G"

# alignment metrics including mismatch and error rates
rule align_mismatch:
    input:
        fa=f"{data_repo}/reference/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa",
        bam=f"{data_repo}/input/{{cell_line}}/sorted_{{SRA}}.bam",
    output:
        f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/alignment_summary.txt"
    log:
        f"{data_repo}/logs/qualiAlign/{{cell_line}}/{{SRA}}.log"
    threads: 4
    shell:
        "java -Xmx24g -jar picard.jar CollectAlignmentSummaryMetrics R={input.fa} I={input.bam} O={output} > {log}"

# insert size metrics
rule insertMetrics:
    input:
        f"{data_repo}/input/{{cell_line}}/sorted_{{SRA}}.bam"
    output:
        metrics=f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/insert_size_metrics.txt",
        histogram=f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/insert_size_histogram.pdf"
    log:
        f"{data_repo}/logs/insertMetrics/{{cell_line}}/{{SRA}}.log"
    shell:
        "java -Xmx24g -jar picard.jar CollectInsertSizeMetrics I={input} O={output.metrics} H={output.histogram} > {log}"

rule aggregateQC:
    input:
        [
            f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/flagstat_output.txt",
            f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/alignment_stats.txt",
            f"{data_repo}/output/{{cell_line}}/tumor_coverage.bedgraph",
            f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/qualimap_report/",
            f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/alignment_summary.txt",
            f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/insert_size_metrics.txt",
            f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentQC/insert_size_histogram.pdf",
        ]
    output:
        f"{data_repo}/output/{{cell_line}}/{{SRA}}/alignmentMultiqc_report.html"
    threads: 4
    shell:
        f"""
        multiqc {{input}} -o {data_repo}/output/{{wildcards.cell_line}}/{{wildcards.SRA}}/ --filename alignmentMultiqc_report.html
        """

