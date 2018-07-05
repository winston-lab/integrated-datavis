library(cowplot)
library(tidyverse)
library(forcats)
library(viridis)
library(seriation)
library(psych)
library(ggthemes)
library(gtable)

format_xaxis = function(refptlabel, upstream, dnstream){
    function(x){
        if (first(upstream)>500 | first(dnstream)>500){
            return(if_else(x==0, refptlabel, as.character(x)))
        } else {
            return(if_else(x==0, refptlabel, as.character(x*1000)))
        }
    }
}

theme_heatmap = theme_minimal() +
    theme(text = element_text(size=12, color="black", face="bold"),
          axis.text.y = element_blank(),
          axis.text.x = element_text(size=12, color="black", face="bold", margin = unit(c(0,0,0,0), "cm")),
          axis.title.y = element_blank(),
          strip.text = element_text(size=12, color="black", face="bold"),
          strip.text.y = element_blank(),
          # strip.background = element_rect(fill="white", size=0),
          strip.background = element_blank(),
          legend.position="top",
          legend.text = element_text(size=10, face="plain"),
          legend.margin = margin(0,0,0,0),
          legend.box.margin = margin(0,0,0,0),
          panel.grid.major.x = element_line(color="black", size=1.5),
          panel.grid.minor.x = element_line(color="black"),
          panel.grid.major.y = element_blank(),
          panel.grid.minor.y = element_blank(),
          panel.spacing.x = unit(.25, "cm"))

theme_metagene = theme_light() +
    theme(text = element_text(size=12, color="black", face="bold"),
          axis.text = element_text(size=12, color="black", face="bold"),
          axis.text.y = element_text(size=10, face="plain"),
          axis.title = element_text(size=10, face="plain"),
          strip.text = element_text(size=12, color="black", face="bold"),
          strip.background = element_rect(fill="white", size=0),
          strip.placement="outside",
          legend.position = "top",
          legend.title = element_blank(),
          legend.text = element_text(size=12),
          panel.spacing.x = unit(0.5, "cm"))

main = function(inputs, anno_paths, conditions, cutoff_pcts, trim_pcts, logtxn, pcount, assays, ptype, refptlabel,
                upstream, dnstream, scaled_length, sortmethod, cluster_assays, cluster_five, cluster_three, k,
                cmap, endlabel, anno_out, cluster_out, heatmap_out, meta_sample_byannotation_out,
                meta_group_byannotation_out, meta_group_bycondition_out) {
    n_assays = length(assays)

    dflist = list()
    cluster_df = tibble()

    #import data
    for (i in 1:n_assays){
        dflist[[assays[[i]]]] =
            read_tsv(inputs[[i]], col_names=c('group', 'sample', 'annotation',
                                              'assay', 'index', 'position', 'signal')) %>%
            group_by(annotation) %>%
            mutate(annotation_labeled=paste(n_distinct(index), annotation)) %>%
            ungroup() %>%
            mutate(annotation=annotation_labeled) %>%
            select(-annotation_labeled) %>%
            mutate_at(vars('group', 'sample', 'annotation'), funs(fct_inorder(., ordered=TRUE)))

        #cluster on data averaged over group only
        #scale data 0 to 1 for each annotation to cluster more on position and less on signal levels,
        if (sortmethod=="cluster"){
            cluster_df = dflist[[i]] %>%
                filter(paste0(group, "-", assay) %in% cluster_assays &
                           between(position, cluster_five/1000, cluster_three/1000)) %>%
                group_by(group, annotation, assay, index, position) %>%
                summarise(mean=mean(signal)) %>%
                group_by(group, annotation, assay, index) %>%
                mutate(mean=scales::rescale(mean)) %>%
                ungroup() %>%
                bind_rows(cluster_df, .)
        }
    }

    #import annotation information
    n_anno = n_distinct(dflist[[1]][["annotation"]], na.rm=TRUE)
    annotations = dflist[[1]] %>% distinct(annotation) %>% filter(!is.na(annotation)) %>%  pull(annotation)
    bed = tibble()
    for (i in 1:n_anno){
        bed = read_tsv(anno_paths[i], col_names=c('chrom', 'start', 'end', 'name', 'score', 'strand')) %>%
            mutate(annotation=annotations[i],
                   score = as.character(score)) %>%
            rowid_to_column(var="index") %>%
            bind_rows(bed, .)
    }

    #sort (by clustering, annotation length, or keep original order)
    if (sortmethod=="cluster"){
        reorder = tibble()

        for (i in 1:n_anno){
            rr = cluster_df %>% filter(annotation==annotations[i]) %>%
                select(-annotation) %>%
                unite(cid, c("group", "assay", "position"), sep="~") %>%
                spread(cid, mean, fill=0) %>%
                remove_rownames() %>%
                column_to_rownames(var="index")

            d = dist(rr, method="euclidean")
            l = kmeans(d, k[i])[["cluster"]]

            pdf(file=cluster_out[i], width=6, height=6)
            unsorted = dissplot(d, method=NA, newpage=TRUE,
                                main=paste0(annotations[i], " unsorted"),
                                options=list(zlim=c(0, quantile(d, probs=0.97)),
                                             col=viridis(100, direction=-1)))

            if (k[i] > 1){
                seriated = dissplot(d, labels=l, method="OLO",
                                    options=list(silhouettes=TRUE,
                                                 zlim=c(0, quantile(d, probs=0.97)),
                                                 col=viridis(100, direction=-1)))
                dev.off()
                sub_reorder = tibble(annotation=annotations[i],
                                     cluster=seriated[["labels"]],
                                     og_index=seriated[["order"]]) %>%
                    mutate(new_index = row_number())
            } else if (k[i]==1){
                seriated = seriate(d, method="OLO")
                sub_reorder = tibble(annotation=annotations[i],
                                     cluster=as.integer(1),
                                     og_index=get_order(seriated)) %>%
                    mutate(new_index = row_number())
                dev.off()
            }

            reorder = reorder %>% bind_rows(sub_reorder)

            sorted = sub_reorder %>% left_join(bed, by=c("annotation", "og_index"="index")) %>%
                select(-c(annotation, og_index, new_index))
            for (j in 1:k[i]){
                sorted %>% filter(cluster==j) %>%
                    select(-cluster) %>%
                    write_tsv(anno_out[sum(k[0:(i-1)])+j], col_names=FALSE)
            }
        }
        for (m in 1:n_assays){
            dflist[[m]] = reorder %>%
                group_by(annotation, cluster) %>%
                mutate(cluster_labeled = paste0("cluster ", cluster, " (n=", n_distinct(new_index), ")")) %>%
                ungroup() %>% mutate(cluster = cluster_labeled) %>% select(-cluster_labeled) %>%
                left_join(dflist[[m]], ., by=c("annotation", "index"="og_index")) %>%
                group_by(annotation, cluster) %>%
                mutate(new_index = as.integer(new_index+1-min(new_index))) %>%
                ungroup() %>%
                arrange(annotation, cluster, new_index)
        }
    } else if (sortmethod=="length"){
        sorted = bed %>% group_by(annotation) %>%
            arrange(end-start, .by_group=TRUE) %>%
            rowid_to_column(var="new_index") %>%
            mutate(new_index = as.integer(new_index+1-min(new_index))) %>%
            ungroup()

        for (i in 1:n_anno){
            sorted %>% filter(annotation==annotations[i]) %>%
                select(-c(new_index, index, annotation)) %>%
                write_tsv(path=anno_out[i], col_names=FALSE)
        }

        for (m in 1:n_assays){
            dflist[[m]] = sorted %>%
                select(index, new_index, annotation) %>%
                right_join(dflist[[m]], by=c("annotation", "index")) %>%
                mutate(cluster = as.integer(1))
        }
    } else {
        for (i in 1:n_anno){
            bed %>% filter(annotation==annotations[i]) %>%
                select(-c(index, annotation)) %>%
                write_tsv(path=anno_out[i], col_names=FALSE)
        }

        for (m in 1:n_assays){
            dflist[[m]] = dflist[[m]] %>%
                mutate(new_index = index, cluster=as.integer(1))
        }
    }

    heatmaps = list()

    for (i in 1:n_assays){
        hmap_df = dflist[[i]]

        #account for datasets missing conditions
        if (! all(conditions %in% hmap_df[["group"]])) {
            hmap_df = hmap_df %>% complete(group=conditions,
                                           fill=list(annotation=hmap_df %>% slice(1) %>% pull(annotation),
                                                     cluster=hmap_df %>% slice(1) %>% pull(cluster))) %>%
                mutate(group = ordered(group, levels=conditions))
        }

        hmap_df = hmap_df %>%
            group_by(group, annotation, assay, position, cluster, new_index) %>%
            summarise(mean = mean(signal))

        #if sortmethod isn't length, fill in missing data with minimum signal
        if (sortmethod != "length"){
            hmap_df = hmap_df %>% group_by(group, annotation, assay, cluster) %>%
                complete(new_index, position, fill=list(mean=min(hmap_df[["mean"]], na.rm=TRUE)))
        }

        climit = hmap_df %>% filter(mean > 0) %>% pull(mean) %>% quantile(probs=cutoff_pcts[i], na.rm=TRUE)

        if (logtxn[i]) {
            heatmaps[[assays[i]]] =
                ggplot(data = hmap_df, aes(x=position, y=new_index, fill=log2(mean+pcount[i])))
        } else {
            heatmaps[[assays[i]]] =
                ggplot(data = hmap_df, aes(x=position, y=new_index, fill=mean))
        }
        heatmaps[[i]] = heatmaps[[i]] +
            geom_raster() +
            scale_fill_viridis(option=cmap, na.value = "#FFFFFF00",
                               limits = c(NA, climit), oob=scales::squish,
                               guide=guide_colorbar(title.position="top", barwidth=12,
                                                    barheight=1, title.hjust=0.5),
                               name=if (logtxn[[i]]){ bquote(bold(log[2] ~ .(assays[[i]]) ~ "signal"))}
                                    else {bquote(bold(.(assays[[i]]) ~ "signal"))}) +
            scale_y_reverse(expand=c(0.02, 0)) +
            theme_heatmap
        if(max(k)>1){
            heatmaps[[i]] = heatmaps[[i]] +
                facet_grid(annotation+cluster~group, scale="free_y", space="free_y", switch="y", drop=FALSE)
        } else {
            heatmaps[[i]] = heatmaps[[i]] +
                facet_grid(annotation~group, scale="free_y", space="free_y", switch="y", drop=FALSE)
        }

        if (ptype=="absolute"){
            heatmaps[[i]] = heatmaps[[i]] +
                scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                                   labels = format_xaxis(refptlabel = refptlabel,
                                                         upstream = upstream,
                                                         dnstream = dnstream),
                                   name= paste("distance from", refptlabel,
                                               if_else(upstream>500 | dnstream>500, "(kb)", "(nt)")),
                                   expand=c(0.05, 0))
        } else {
            heatmaps[[i]] = heatmaps[[i]] +
                scale_x_continuous(breaks = c(0, (scaled_length/2)/1000, scaled_length/1000),
                                   labels = c(refptlabel, "", endlabel),
                                   name= "scaled distance",
                                   expand=c(0.05, 0))
        }
    }

    #EXTREMELY JANKY WAY TO GET FACET LABELS...build an invisible plot...
    #TODO: facet by clusters also and nicer labeling for clusters
    facet_label = ggplot(data = hmap_df, aes(x=0, y=new_index, fill=mean)) +
        geom_raster() +
        scale_fill_gradient(low="#FFFFFF00", high="#FFFFFF00",
                            guide=guide_colorbar(title.position="top", barwidth=0.1,
                                                barheight=1, title.hjust=0.5),
                           name=bquote(bold("."))) +
        scale_y_reverse(expand=c(0.02, 0)) +
        facet_grid(annotation~group, scale="free_y", space="free_y") +
        theme_heatmap +
        theme(text = element_text(color="#FFFFFF00"),
              strip.text.x = element_text(color="#FFFFFF00"),
              axis.text.x = element_text(color="#FFFFFF00"),
              panel.grid = element_blank(),
              strip.text.y = element_text(size=12, color="black", face="bold", angle=0, hjust=1))

    for (i in 0:((n_assays-1) %/% 4)){
        heatmaps = append(heatmaps, list(facet_label), after=4*i+i)
    }

    all_heatmaps = plot_grid(plotlist = heatmaps, align="h", ncol=min(n_assays+1, 5), axis="trbl")

    ggplot2::ggsave(heatmap_out, plot = all_heatmaps,
           width=2+2/15*max(nchar(as.character(annotations)))+min(n_assays,4)*16,
           height=ceiling(n_assays/4)*25, units="cm", limitsize=FALSE)

    metadf_sample = tibble()
    metadf_group = tibble()

    for (i in 1:n_assays){
        metadf_sample = dflist[[i]] %>%
            group_by(group, sample, annotation, assay, position, cluster) %>%
            summarise(mid = median(signal, na.rm=TRUE),
                      low = quantile(signal, trim_pcts[i], na.rm=TRUE),
                      high = quantile(signal, 1-trim_pcts[i], na.rm=TRUE)) %>%
            ungroup() %>%
            mutate(y_min = min(low, na.rm = TRUE),
                   y_range = max(high, na.rm=TRUE)-y_min) %>%
            mutate_at(vars(low, mid, high), funs((.-y_min)/y_range)) %>%
            select(-c(y_min, y_range)) %>%
            bind_rows(metadf_sample, .)

        metadf_group = dflist[[i]] %>%
            mutate(signal = scales::rescale(signal)) %>%
            group_by(group, annotation, assay, position, cluster) %>%
            summarise(mid = median(signal),
                      low = quantile(signal, trim_pcts[i]),
                      high = quantile(signal, 1-trim_pcts[i])) %>%
            ungroup() %>%
            mutate(y_min = min(low, na.rm = TRUE),
                   y_range = max(high, na.rm=TRUE)-y_min) %>%
            mutate_at(vars(low, mid, high), funs((.-y_min)/y_range)) %>%
            select(-c(y_min, y_range)) %>%
            bind_rows(metadf_group, .)
    }

    if (sortmethod != "cluster"){
        metadf_sample = metadf_sample %>% mutate(cluster = paste("cluster", cluster))
        metadf_group = metadf_group %>% mutate(cluster = paste("cluster", cluster))
    }
    metadf_sample = metadf_sample %>% mutate_at(vars(group, sample, assay, cluster), funs(fct_inorder(., ordered=TRUE)))
    metadf_group = metadf_group %>% mutate_at(vars(group, assay, cluster), funs(fct_inorder(., ordered=TRUE)))

    nest_right_facets = function(ggp, level=2, outer="replicate", inner="annotation"){
        og_grob = ggplotGrob(ggp)
        strip_loc = grep("strip-r", og_grob[["layout"]][["name"]])
        strip = gtable_filter(og_grob, "strip-r", trim=FALSE)
        strip_heights = gtable_filter(og_grob, "strip-r")[["heights"]]

        strip_top = min(strip[["layout"]][["t"]])
        strip_bot = max(strip[["layout"]][["b"]])
        strip_x = strip[["layout"]][["r"]][1]

        mat = matrix(vector("list", length=(length(strip)*2-1)*level), ncol=level)
        mat[] = list(zeroGrob())

        facet_grob = gtable_matrix("rightcol", grobs=mat,
                                   widths=unit(rep(1,level), "null"),
                                   heights=strip_heights)

        if(level==3){
            rep_grob_indices = seq(1, length(strip_loc), sum(k))
            for (rep_idx in 1:max_reps){
                #add replicate facet label
                facet_grob = facet_grob %>%
                        gtable_add_grob(grobs = og_grob$grobs[[strip_loc[rep_grob_indices[rep_idx]]]]$grobs[[level]],
                                        t = ((sum(k)*2))*(rep_idx-1)+1,
                                        b = ((sum(k)*2))*(rep_idx)-1,
                                        l = level, r = level)
                #for each annotation within each replicate
                for (anno_idx in 1:n_anno){
                    t = ((sum(k)*2))*(rep_idx-1)+1+sum(k[1:anno_idx])-k[1]+2*(anno_idx-1)
                    b = t + k[anno_idx]
                    facet_grob = facet_grob %>%
                            gtable_add_grob(grobs = og_grob$grobs[[strip_loc[rep_grob_indices[rep_idx]]+
                                                                       sum(k[1:anno_idx])-k[1]]]$grobs[[2]],
                                            t = t, b = b, l = 2, r = 2)
                }
            }
        } else if(level==2){
            if (outer=="annotation"){
                outer_grob_indices = 1+lag(k, default=0)
                n_outer = n_anno
            } else if (outer=="replicate"){
                outer_grob_indices = seq(1, length(strip_loc), sum(k))
                n_outer = max_reps
            }
            for (idx in 1:n_outer){
                if (outer=="annotation"){
                    t=((k[idx]*2))*(idx-1)+1
                    b=((k[idx]*2))*(idx)-1
                } else if (outer=="replicate"){
                    if (inner=="cluster"){
                        t=((k*2))*(idx-1)+1
                        b=((k*2))*(idx)-1
                    } else{
                        t = (n_anno*2)*(idx-1)+1
                        b = (n_anno*2)*(idx)-1
                    }
                }
                facet_grob = facet_grob %>%
                    gtable_add_grob(grobs = og_grob$grobs[[strip_loc[outer_grob_indices[idx]]]]$grobs[[2]],
                                    t=t, b=b, l=2, r=2)
            }
        }
        new_grob = gtable_add_grob(og_grob, facet_grob, t=strip_top, r=strip_x, l=strip_x, b=strip_bot, name='rstrip')
        return(new_grob)
    }

    meta = function(ggp, nest_right=TRUE){
        ggp = ggp +
            geom_vline(xintercept=0, size=1, color="grey65")
        if (ptype=="scaled"){
            ggp = ggp +
                geom_vline(xintercept=scaled_length, size=1, color="grey65")
        }
        ggp = ggp +
            geom_ribbon(alpha=0.4, size=0) +
            geom_line() +
            scale_y_continuous(name="relative signal", breaks = scales::pretty_breaks(n=3)) +
            scale_color_ptol(guide=guide_legend(keywidth = unit(2, "cm"), label.position="top", label.hjust=0.5)) +
            scale_fill_ptol() +
            theme_metagene
        if (ptype=="absolute"){
            ggp = ggp +
                scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                                   labels= function(x){if_else(x==0, refptlabel,
                                                               if(upstream>500 | dnstream>500){as.character(x)} else {as.character(x*1000)})},
                                   name=paste("distance from", refptlabel, if(upstream>500 | dnstream>500){"(kb)"} else {"(nt)"}),
                                   limits = c(-upstream/1000, dnstream/1000),
                                   expand=c(0,0))
        } else {
            ggp = ggp +
                scale_x_continuous(breaks=c(0, (scaled_length/2)/1000, scaled_length/1000),
                                   labels=c(refptlabel, "", endlabel),
                                   name="scaled distance",
                                   limits = c(-upstream/1000, (scaled_length+dnstream)/1000),
                                   expand=c(0,0))
        }
        if (nest_right){
            ggp = ggp %>% nest_right_facets(level=2, outer="annotation", inner="cluster")
        }
        return(ggp)
    }

    meta_sample_byannotation = ggplot(data = metadf_sample,
                                      aes(x=position, y=mid, ymin=low, ymax=high,
                                          group=sample, color=group, fill=group))

    meta_group_byannotation = ggplot(data = metadf_group,
                                     aes(x=position, y=mid, ymin=low, ymax=high,
                                         group=group, color=group, fill=group))
    if (n_anno==1 && max(k)==1){
        meta_sample_byannotation = meta(meta_sample_byannotation + facet_wrap(~assay, ncol=4),
                                        nest_right=FALSE)
        meta_group_byannotation = meta(meta_group_byannotation + facet_wrap(~assay, ncol=4),
                                       nest_right=FALSE)
    } else {
        meta_sample_byannotation = meta(meta_sample_byannotation + facet_grid(annotation+cluster~assay))
        meta_group_byannotation = meta(meta_group_byannotation + facet_grid(annotation+cluster~assay))
    }

    ggplot2::ggsave(meta_sample_byannotation_out,
                    plot=meta_sample_byannotation,
                    width= if(n_anno==1 && max(k)==1){42} else {2+n_assays*10},
                    height= if(n_anno==1 && max(k)==1){2+8.5*ceiling(n_assays/4)} else {3+sum(k)*4.5},
                    units="cm", limitsize=FALSE)
    ggplot2::ggsave(meta_group_byannotation_out,
                    plot=meta_group_byannotation,
                    width= if(n_anno==1 && max(k)==1){42} else {2+n_assays*10},
                    height= if(n_anno==1 && max(k)==1){2+8.5*ceiling(n_assays/4)} else {3+sum(k)*4.5},
                    units="cm", limitsize=FALSE)

    meta_group_bycondition = (ggplot(data = metadf_group %>% mutate(clabel = fct_inorder(paste0(annotation, ", ", cluster), ordered=TRUE)),
                aes(x=position, y=mid, ymin=low, ymax=high,
                    group=clabel, fill=clabel, color=clabel)) +
        facet_grid(group~assay, switch="y")) %>% meta(nest_right=FALSE) +
        scale_color_ptol(guide=guide_legend(ncol=1)) +
        scale_fill_ptol() +
        theme(legend.position="bottom",
              strip.text.y = element_text(angle=180, hjust=1),
              plot.margin = margin(0.5, 1, 0.5, 0.5, "cm"))

    ggplot2::ggsave(meta_group_bycondition_out, plot=meta_group_bycondition, width=2+n_assays*10, height=1.5*sum(k)+4+4.5*n_distinct(dflist[[1]][["group"]]), units="cm", limitsize=FALSE)
}

main(inputs = snakemake@input[["matrices"]],
     anno_paths = snakemake@input[["annotations"]],
     conditions = snakemake@params[["conditions"]],
     cutoff_pcts = snakemake@params[["cutoffs"]],
     trim_pcts = snakemake@params[["trim_pct"]],
     logtxn = snakemake@params[["logtxn"]],
     pcount = snakemake@params[["pcount"]],
     assays = snakemake@params[["assays"]],
     ptype = snakemake@params[["mtype"]],
     refptlabel = snakemake@params[["refptlabel"]],
     upstream = snakemake@params[["upstream"]],
     dnstream= snakemake@params[["dnstream"]],
     scaled_length= snakemake@params[["scaled_length"]],
     endlabel = snakemake@params[["endlabel"]],
     sortmethod = snakemake@params[["sortmethod"]],
     cluster_assays = snakemake@params[["cluster_using"]],
     cluster_five = snakemake@params[["cluster_five"]],
     cluster_three = snakemake@params[["cluster_three"]],
     k = snakemake@params[["k"]],
     cmap= snakemake@params[["cmap"]],
     anno_out = snakemake@params[["annotations_out"]],
     cluster_out = snakemake@params[["cluster_out"]],
     heatmap_out = snakemake@output[["heatmap"]],
     meta_sample_byannotation_out = snakemake@output[["sample_facet_anno"]],
     meta_group_byannotation_out = snakemake@output[["group_facet_anno"]],
     meta_group_bycondition_out = snakemake@output[["group_facet_group"]])
