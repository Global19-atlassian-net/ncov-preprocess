import os


rule upload:
    input:
        "results/masked.fasta.gz",
        "results/aligned.fasta.gz",
        "results/filtered.fasta.gz",
        "results/sequence-diagnostics.tsv",
        "results/flagged-sequences.tsv",
        "results/to-exclude.txt"
    params:
        s3_bucket = config["S3_BUCKET"]
    run:
        for fname in input:
            shell(f"./scripts/upload-to-s3 {fname} {params.s3_bucket}/{os.path.basename(fname)}")




rule download:
    message: "Downloading metadata and fasta files from S3"
    output:
        sequences = config["sequences"],
        metadata = config["metadata"]
    conda: config["conda_environment"]
    shell:
        """
        aws s3 cp s3://nextstrain-ncov-private/metadata.tsv.gz - | gunzip -cq >{output.metadata:q}
        aws s3 cp s3://nextstrain-ncov-private/sequences.fasta.gz - | gunzip -cq > {output.sequences:q}
        """

rule filter:
    message:
        """
        Filtering to
          - excluding strains in {input.exclude}
        """
    input:
        sequences = rules.download.output.sequences,
        metadata = rules.download.output.metadata,
        exclude = config["files"]["exclude"]
    output:
        sequences = "results/filtered.fasta"
    log:
        "logs/filtered.txt"
    params:
        min_length = config["filter"]["min_length"],
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --min-length {params.min_length} \
            --output {output.sequences} 2>&1 | tee {log}
        """

rule align:
    message:
        """
        Aligning sequences to {input.reference}
          - gaps relative to reference are considered real
        """
    input:
        sequences = "results/filtered.fasta",
        reference = config["files"]["alignment_reference"]
    output:
        alignment = "results/aligned.fasta"
    log:
        "logs/align.txt"
    benchmark:
        "benchmarks/align.txt"
    threads: 16
    conda: config["conda_environment"]
    shell:
        """
        mafft \
            --anysymbol \
            --auto \
            --thread {threads} \
            --keeplength \
            --addfragments \
            {input.sequences} \
            {input.reference} > {output} 2> {log}
        """


rule diagnostic:
    message: "Scanning aligned sequences {input.alignment} for problematic sequences"
    input:
        alignment = rules.align.output.alignment,
        metadata = rules.download.output.metadata,
        reference = config["files"]["alignment_reference"]
    output:
        diagnostics = "results/sequence-diagnostics.tsv",
        flagged = "results/flagged-sequences.tsv",
        to_exclude = "results/to-exclude.txt"
    log:
        "logs/diagnostics.txt"
    params:
        mask_from_beginning = config["mask"]["mask_from_beginning"],
        mask_from_end = config["mask"]["mask_from_end"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/diagnostic.py \
            --alignment {input.alignment} \
            --metadata {input.metadata} \
            --reference {input.reference} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --output-flagged {output.flagged} \
            --output-diagnostics {output.diagnostics} \
            --output-exclusion-list {output.to_exclude} 2>&1 | tee {log}
        """

rule refilter:
    message:
        """
        excluding sequences flagged in the diagnostic step in file {input.exclude}
        """
    input:
        sequences = rules.align.output.alignment,
        metadata = rules.download.output.metadata,
        exclude = rules.diagnostic.output.to_exclude
    output:
        sequences = "results/aligned-filtered.fasta"
    log:
        "logs/refiltered.txt"
    conda: config["conda_environment"]
    shell:
        """
        augur filter \
            --sequences {input.sequences} \
            --metadata {input.metadata} \
            --exclude {input.exclude} \
            --output {output.sequences} 2>&1 | tee {log}
        """


rule mask:
    message:
        """
        Mask bases in alignment
          - masking {params.mask_from_beginning} from beginning
          - masking {params.mask_from_end} from end
          - masking other sites: {params.mask_sites}
        """
    input:
        alignment = rules.refilter.output.sequences
    output:
        alignment = "results/masked.fasta"
    log:
        "logs/mask.txt"
    params:
        mask_from_beginning = config["mask"]["mask_from_beginning"],
        mask_from_end = config["mask"]["mask_from_end"],
        mask_sites = config["mask"]["mask_sites"]
    conda: config["conda_environment"]
    shell:
        """
        python3 scripts/mask-alignment.py \
            --alignment {input.alignment} \
            --mask-from-beginning {params.mask_from_beginning} \
            --mask-from-end {params.mask_from_end} \
            --mask-sites {params.mask_sites} \
            --mask-terminal-gaps \
            --output {output.alignment} 2>&1 | tee {log}
        """

rule all:
    input:
        "results/masked.fasta",
        "results/aligned.fasta",
        "results/filtered.fasta"
    output:
        "results/masked.fasta.gz",
        "results/aligned.fasta.gz",
        "results/filtered.fasta.gz"
    shell:
        '''
        gzip {input}
        '''

