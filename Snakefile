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
        # expand("datavis/{figure}/{figure}-metagene-sample-by-anno.svg", figure=FIGURES)

rule make_stranded_annotations:
    input:
        lambda wc : FIGURES[wc.figure]["annotations"][wc.annotation]["path"]
    output:
        "annotations/{figure}/{annotation}-STRANDED.{ext}"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{figure}_{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
        """

rule compute_matrix:
    input:
        annotation = lambda wc: "annotations/" + wc.figure + "/" + wc.annotation + "-STRANDED" + os.path.splitext(FIGURES[wc.figure]["annotations"][wc.annotation]["path"])[1] if ASSAYS[wc.assay]["stranded"] else FIGURES[wc.figure]["annotations"][wc.annotation]["path"],
        bw = lambda wc: ASSAYS[wc.assay]["coverage"][wc.sample]["path"]
    output:
        dtfile = temp("datavis/{figure}/{annotation}_{assay}_{sample}.mat.gz"),
        matrix = temp("datavis/{figure}/{annotation}_{assay}_{sample}.tsv"),
        melted = temp("datavis/{figure}/{annotation}_{assay}_{sample}-melted.tsv.gz"),
    params:
        group = lambda wc: ASSAYS[wc.assay]["coverage"][wc.sample]["group"],
        refpoint = lambda wc: "TSS" if FIGURES[wc.figure]["parameters"]["type"]=="scaled" else FIGURES[wc.figure]["parameters"]["refpoint"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"] + ASSAYS[wc.assay]["binsize"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["downstream"] + ASSAYS[wc.assay]["binsize"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        binsize = lambda wc: ASSAYS[wc.assay]["binsize"],
        binstat = lambda wc: ASSAYS[wc.assay]["binstat"],
        nan_afterend = lambda wc: [] if FIGURES[wc.figure]["parameters"]["type"]=="scaled" or not FIGURES[wc.figure]["parameters"]["nan_afterend"] else "--nanAfterEnd",
        anno_label = lambda wc: FIGURES[wc.figure]["annotations"][wc.annotation]["label"]
    threads: config["threads"]
    log: "logs/compute_matrix/compute_matrix-{figure}_{annotation}_{sample}_{assay}.log"
    run:
        if FIGURES[wildcards.figure]["parameters"]["type"]=="absolute":
            shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} {params.nan_afterend} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {params.scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {params.refpoint} -g {params.group} -s {wildcards.sample} -a {params.anno_label} -y {wildcards.assay} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        lambda wc: expand("datavis/{figure}/{annotation}_{assay}_{sample}-melted.tsv.gz", annotation = [k for k,v in FIGURES[wc.figure]["annotations"].items()], sample = [k for k,v in ASSAYS[wc.assay]["coverage"].items() if v["group"]==FIGURES[wc.figure]["control"] or v["group"]==FIGURES[wc.figure]["condition"]], figure=wc.figure, assay=wc.assay)
    output:
        "datavis/{figure}/{figure}_{assay}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{figure}_{assay}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_figures:
    input:
        matrices = lambda wc: expand("datavis/{figure}/{figure}_{assay}.tsv.gz", figure=wc.figure, assay=FIGURES[wc.figure]["include"]),
        annotations = lambda wc: [v["path"] for k,v in FIGURES[wc.figure]["annotations"].items()]
    output:
        heatmap = "datavis/{figure}/{figure}-heatmaps.svg",
        sample_facet_anno = "datavis/{figure}/{figure}-metagene-sample-facet-anno.svg",
        group_facet_anno = "datavis/{figure}/{figure}-metagene-group-facet-anno.svg",
        group_facet_group = "datavis/{figure}/{figure}-metagene-group-facet-group.svg",
    params:
        annotations_out = lambda wc: ["datavis/" + wc.figure + "/" + annotation + "_cluster-" + str(cluster) + ".bed" for annotation in FIGURES[wc.figure]["annotations"] for cluster in range(1, FIGURES[wc.figure]["annotations"][annotation]["n_clusters"]+1)],
        clusters_out = lambda wc: ["datavis/" + wc.figure + "/" + annotation + ".pdf" for annotation in FIGURES[wc.figure]["annotations"]],
        conditions = lambda wc: [FIGURES[wc.figure]["control"], FIGURES[wc.figure]["condition"]],
        cutoffs = lambda wc: [ASSAYS[a]["cutoff"] for a in FIGURES[wc.figure]["include"]],
        trim_pct = lambda wc: [ASSAYS[a]["trim_pct"] for a in FIGURES[wc.figure]["include"]],
        # logtxn = lambda wc: [str(ASSAYS[a]["log_transform"]).upper() for a in FIGURES[wc.figure]["include"]],
        logtxn = lambda wc: [ASSAYS[a]["log_transform"] for a in FIGURES[wc.figure]["include"]],
        pcount = lambda wc: [0 if not ASSAYS[a]["log_transform"] else ASSAYS[a]["pseudocount"] for a in FIGURES[wc.figure]["include"]],
        assays = lambda wc: [ASSAYS[a]["label"] for a in FIGURES[wc.figure]["include"]],
        mtype = lambda wc: FIGURES[wc.figure]["parameters"]["type"],
        refptlabel = lambda wc: FIGURES[wc.figure]["parameters"]["refpointlabel"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["downstream"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        endlabel = lambda wc:  "HAIL SATAN" if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["endlabel"],
        cmap = lambda wc: FIGURES[wc.figure]["parameters"]["heatmap_colormap"],
        sortmethod = lambda wc: FIGURES[wc.figure]["parameters"]["arrange"],
        cluster_using = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else [k + "-"  + assay for k,v in FIGURES[wc.figure]["parameters"]["cluster_using"].items() for assay in v],
        cluster_five = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_five"],
        cluster_three = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_three"],
        k = lambda wc: [v["n_clusters"] for k,v in FIGURES[wc.figure]["annotations"].items()]
    # run:
    #     print(input)
    #     print(params)
    script:
        "scripts/integrated_figures.R"

# rule plot_heatmaps:
#     input:
#         lambda wc: expand("datavis/{figure}/{figure}_{assay}.tsv.gz", figure=wc.figure, assay=FIGURES[wc.figure]["include"])
#     output:
#         "datavis/{figure}/{figure}-heatmaps.svg"
#     params:
#         cutoffs = lambda wc: [ASSAYS[a]["cutoff"] for a in FIGURES[wc.figure]["include"]],
#         logtxn = lambda wc: [ASSAYS[a]["log_transform"] for a in FIGURES[wc.figure]["include"]],
#         assays = lambda wc: [ASSAYS[a]["label"] for a in FIGURES[wc.figure]["include"]],
#         mtype = lambda wc: FIGURES[wc.figure]["parameters"]["type"],
#         refptlabel = lambda wc: FIGURES[wc.figure]["parameters"]["refpointlabel"],
#         upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"],
#         dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["downstream"],
#         cluster = lambda wc: FIGURES[wc.figure]["parameters"]["cluster"],
#     run:
#         if FIGURES[wc.figure]["parameters"]["type"]=="scaled":
#             scaled_length = FIGURES[wildcards.figure]["parameters"]["scaled_length"]
#             endlabel = FIGURES[wildcards.figure]["parameters"]["endlabel"]
#         else:
#             scaled_length=0
#             endlabel = "WINGARDIUM LEVIOSA"
#         if FIGURES[wildcards.figure]["parameters"]["cluster"]:
#             cluster_assays = [ASSAYS[a]["label"] for a in FIGURES[wildcards.figure]["parameters"]["cluster_using"]]
#             k = FIGURES[wildcards.figure]["parameters"]["k"]
#         else:
#             cluster_assays = "EVERYTHING-seq"
#             k = 0
#         shell("""Rscript scripts/integrated_heatmaps.R -i {input} -c {params.cutoffs} -l {params.logtxn} -a {params.assays} -t {params.mtype} -r {params.refptlabel} -u {params.upstream} -d {params.dnstream} -s {scaled_length} -e {endlabel} -z {params.cluster} -y {cluster_assays} -k {k} -o {output}""")


# rule plot_metagenes:
#     input:
#         lambda wc: expand("datavis/{figure}/{figure}_{assay}.tsv.gz", figure=wc.figure, assay=FIGURES[wc.figure]["include"])
#     output:
#         sample_facet_anno = "datavis/{figure}/{figure}-metagene-sample-by-anno.svg",
#         sample_facet_group = "datavis/{figure}/{figure}-metagene-sample-by-group.svg",
#         group_facet_anno = "datavis/{figure}/{figure}-metagene-group-by-anno.svg",
#         group_facet_group = "datavis/{figure}/{figure}-metagene-group-by-group.svg",
#     params:
#         assays = lambda wc: [ASSAYS[a]["label"] for a in FIGURES[wc.figure]["include"]],
#         trim_pct = lambda wc: FIGURES[wc.figure]["parameters"]["trim_pct"],
#         mtype = lambda wc: FIGURES[wc.figure]["parameters"]["type"],
#         refptlabel = lambda wc: FIGURES[wc.figure]["parameters"]["refpointlabel"],
#         upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"],
#         dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["downstream"],
#     run:
#         if FIGURES[wildcards.figure]["parameters"]["type"]=="scaled":
#             scaled_length = FIGURES[wildcards.figure]["parameters"]["scaled_length"]
#             endlabel = FIGURES[wildcards.figure]["parameters"]["endlabel"]
#         else:
#             scaled_length=0
#             endlabel = "WINGARDIUM LEVIOSA"
#         shell("""Rscript scripts/integrated_metagenes.R -i {input} -a {params.assays} -p {params.trim_pct} -t {params.mtype} -r {params.refptlabel} -u {params.upstream} -d {params.dnstream} -s {scaled_length} -e {endlabel} --out1 {output.sample_facet_anno} --out2 {output.sample_facet_group} --out3 {output.group_facet_anno} --out4 {output.group_facet_group}""")

