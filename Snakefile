#!/usr/bin/env python
import os

configfile: "config.yaml"

ASSAYS = config["assays"]
FIGURES = config["figures"]
# ANNOTATIONS = {k:v["path"] for d in [v["annotations"] for k,v in FIGURES.items()] for k,v in d.items()}

localrules:
    make_stranded_annotations,
    cat_matrices

rule all:
    input:
        expand("datavis/{figure}/{figure}-heatmaps.svg", figure=FIGURES),
        expand("datavis/{figure}/{figure}-metagene-sample-by-anno.svg", figure=FIGURES)

rule make_stranded_annotations:
    input:
        lambda wildcards : FIGURES[wildcards.figure]["annotations"][wildcards.annotation]["path"]
    output:
        "annotations/{figure}/{annotation}-STRANDED.{ext}"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{figure}_{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
        """

rule compute_matrix:
    input:
        annotation = lambda wildcards: "annotations/" + wildcards.figure + "/" + wildcards.annotation + "-STRANDED" + os.path.splitext(FIGURES[wildcards.figure]["annotations"][wildcards.annotation]["path"])[1] if ASSAYS[wildcards.assay]["stranded"] else FIGURES[wildcards.figure]["annotations"][wildcards.annotation]["path"],
        bw = lambda wildcards: ASSAYS[wildcards.assay]["coverage"][wildcards.sample]["path"]
    output:
        dtfile = temp("datavis/{figure}/{annotation}_{assay}_{sample}.mat.gz"),
        matrix = temp("datavis/{figure}/{annotation}_{assay}_{sample}.tsv"),
        melted = temp("datavis/{figure}/{annotation}_{assay}_{sample}-melted.tsv.gz"),
    params:
        group = lambda wildcards: ASSAYS[wildcards.assay]["coverage"][wildcards.sample]["group"],
        upstream = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["upstream"] + ASSAYS[wildcards.assay]["binsize"],
        dnstream = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["downstream"] + ASSAYS[wildcards.assay]["binsize"],
        binsize = lambda wildcards: ASSAYS[wildcards.assay]["binsize"],
        sort = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["sort"],
        binstat = lambda wildcards: ASSAYS[wildcards.assay]["binstat"],
        annolabel = lambda wildcards: FIGURES[wildcards.figure]["annotations"][wildcards.annotation]["label"]
    threads: config["threads"]
    log: "logs/compute_matrix/compute_matrix-{figure}_{annotation}_{sample}_{assay}.log"
    run:
        if FIGURES[wildcards.figure]["parameters"]["type"]=="absolute":
            refpoint = FIGURES[wildcards.figure]["parameters"]["refpoint"]
            if FIGURES[wildcards.figure]["parameters"]["nan_afterend"]:
                shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing region_length --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
            else:
                shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing region_length --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            scaled_length = FIGURES[wildcards.figure]["parameters"]["scaled_length"]
            refpoint = "TSS"
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing region_length --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {refpoint} --group {params.group} -s {wildcards.sample} -a {params.annolabel} -y {wildcards.assay} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        lambda wildcards: expand("datavis/{figure}/{annotation}_{assay}_{sample}-melted.tsv.gz", annotation = [k for k,v in FIGURES[wildcards.figure]["annotations"].items()], sample = [k for k,v in ASSAYS[wildcards.assay]["coverage"].items() if v["group"]==FIGURES[wildcards.figure]["control"] or v["group"]==FIGURES[wildcards.figure]["condition"]], figure=wildcards.figure, assay=wildcards.assay)
    output:
        "datavis/{figure}/{figure}_{assay}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{figure}_{assay}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_heatmaps:
    input:
        lambda wildcards: expand("datavis/{figure}/{figure}_{assay}.tsv.gz", figure=wildcards.figure, assay=FIGURES[wildcards.figure]["include"])
    output:
        "datavis/{figure}/{figure}-heatmaps.svg"
    params:
        cutoffs = lambda wildcards: [ASSAYS[a]["cutoff"] for a in FIGURES[wildcards.figure]["include"]],
        logtxn = lambda wildcards: [ASSAYS[a]["log_transform"] for a in FIGURES[wildcards.figure]["include"]],
        assays = lambda wildcards: [ASSAYS[a]["label"] for a in FIGURES[wildcards.figure]["include"]],
        mtype = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["type"],
        refptlabel = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["refpointlabel"],
        upstream = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["upstream"],
        dnstream = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["downstream"],
        cluster = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["cluster"],
    run:
        if FIGURES[wildcards.figure]["parameters"]["type"]=="scaled":
            scaled_length = FIGURES[wildcards.figure]["parameters"]["scaled_length"]
            endlabel = FIGURES[wildcards.figure]["parameters"]["endlabel"]
        else:
            scaled_length=0
            endlabel = "WINGARDIUM LEVIOSA"
        if FIGURES[wildcards.figure]["parameters"]["cluster"]:
            cluster_assays = [ASSAYS[a]["label"] for a in FIGURES[wildcards.figure]["parameters"]["cluster_using"]]
            k = FIGURES[wildcards.figure]["parameters"]["k"]
        else:
            cluster_assays = "EVERYTHING-seq"
            k = 0
        shell("""Rscript scripts/integrated_heatmaps.R -i {input} -c {params.cutoffs} -l {params.logtxn} -a {params.assays} -t {params.mtype} -r {params.refptlabel} -u {params.upstream} -d {params.dnstream} -s {scaled_length} -e {endlabel} -z {params.cluster} -y {cluster_assays} -k {k} -o {output}""")


rule plot_metagenes:
    input:
        lambda wildcards: expand("datavis/{figure}/{figure}_{assay}.tsv.gz", figure=wildcards.figure, assay=FIGURES[wildcards.figure]["include"])
    output:
        sample_facet_anno = "datavis/{figure}/{figure}-metagene-sample-by-anno.svg",
        sample_facet_group = "datavis/{figure}/{figure}-metagene-sample-by-group.svg",
        group_facet_anno = "datavis/{figure}/{figure}-metagene-group-by-anno.svg",
        group_facet_group = "datavis/{figure}/{figure}-metagene-group-by-group.svg",
    params:
        assays = lambda wildcards: [ASSAYS[a]["label"] for a in FIGURES[wildcards.figure]["include"]],
        trim_pct = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["trim_pct"],
        mtype = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["type"],
        refptlabel = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["refpointlabel"],
        upstream = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["upstream"],
        dnstream = lambda wildcards: FIGURES[wildcards.figure]["parameters"]["downstream"],
    run:
        if FIGURES[wildcards.figure]["parameters"]["type"]=="scaled":
            scaled_length = FIGURES[wildcards.figure]["parameters"]["scaled_length"]
            endlabel = FIGURES[wildcards.figure]["parameters"]["endlabel"]
        else:
            scaled_length=0
            endlabel = "WINGARDIUM LEVIOSA"
        shell("""Rscript scripts/integrated_metagenes.R -i {input} -a {params.assays} -p {params.trim_pct} -t {params.mtype} -r {params.refptlabel} -u {params.upstream} -d {params.dnstream} -s {scaled_length} -e {endlabel} --out1 {output.sample_facet_anno} --out2 {output.sample_facet_group} --out3 {output.group_facet_anno} --out4 {output.group_facet_group}""")







