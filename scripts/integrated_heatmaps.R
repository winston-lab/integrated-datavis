library(argparse)
library(cowplot)
library(tidyverse)
library(forcats)
library(viridis)
library(dendsort)

parser = ArgumentParser()
parser$add_argument('-i', '--inputs', type='character', nargs='+')
parser$add_argument('-c', '--cutoffs', type='double', nargs='+')
parser$add_argument('-l', '--logtxn', type='character', nargs='+')
parser$add_argument('-a', '--assays', type='character', nargs='+')
parser$add_argument('-t', '--type', type='character')
parser$add_argument('-r', '--refptlabel', type='character', nargs='+')
parser$add_argument('-u', '--upstream', type='integer')
parser$add_argument('-d', '--dnstream', type='integer')
parser$add_argument('-s', '--scaled_length', type='integer')
parser$add_argument('-e', '--endlabel', type='character', nargs='+')
parser$add_argument('-z', '--cluster', type='character')
parser$add_argument('-y', '--cluster_assays', type='character', nargs='+')
parser$add_argument('-k', dest='k', type='integer', nargs='+')
parser$add_argument('-o', '--output', type='character')

args = parser$parse_args()

format_xaxis = function(refptlabel, upstream, dnstream){
    function(x){
        if (first(upstream)>500 | first(dnstream)>500){
            return(if_else(x==0, refptlabel, as.character(x)))                                                                                                                           
        }
        else {
            return(if_else(x==0, refptlabel, as.character(x*1000)))
        } 
    }                                                                                                                                                                                    
}

main = function(inputs, cutoffs, logtxn, assays, type, refptlabel,
                upstream, dnstream, scaled_length, cluster, cluster_assays, k, endlabel, outpath) {
    nassays = length(assays) 
    dflist = list()
    
    dfcluster = tibble()
    
    for (i in 1:nassays){
        cutoff = cutoffs[[i]]
        dflist[[assays[[i]]]] =
            read_tsv(inputs[[i]], col_names=c('group', 'sample', 'annotation',
                                              'assay', 'index', 'position', 'signal')) %>% 
            mutate(index = paste(annotation, index)) %>% 
            mutate_at(vars('group', 'sample', 'annotation', 'index'), funs(fct_inorder(., ordered=TRUE))) %>% 
            group_by(group, annotation, assay, index, position) %>% 
            summarise(mean = mean(signal)) %>% 
            ungroup()
        
        #for putting number of indices in facet labels
        dd = dflist[[i]] %>% group_by(annotation) %>%
            summarise(n = n_distinct(index)) %>% 
            transmute(annotation=annotation,
                      label = paste(n, annotation))
        
        dflist[[i]] = dflist[[i]] %>% left_join(dd, by="annotation") %>% 
            mutate(annotation = fct_inorder(label, ordered=TRUE)) %>% 
            select(-label)
        
        if (cluster=="True" && assays[[i]] %in% cluster_assays){
            dfcluster = dfcluster %>%
                bind_rows(dflist[[i]] %>%
                              mutate(mean = (mean-min(mean, na.rm=TRUE))/
                                         (max(mean, na.rm=TRUE)-min(mean, na.rm=TRUE))))
        }

        cutoff_val = quantile(dflist[[i]]$mean, probs=cutoff, na.rm=TRUE)
        
        dflist[[i]] = dflist[[i]] %>% mutate(mean = pmin(mean, cutoff_val))
        
        if (logtxn[[i]]=="True"){
            pcount = 0.1
            dflist[[i]] = dflist[[i]] %>% mutate(mean = log2(mean+pcount))    
        }
    }
    
    if (cluster=="True") {
        annotations = unique(dflist[[1]]$annotation)
        roworder = tibble()
        
        for (j in 1:length(annotations)){
            i = annotations[[j]]
            #k-means clustering on non-logtransformed data scaled 0 to 1
            rr = dfcluster %>% filter(annotation==i) %>% 
                select(-annotation) %>% 
                unite(cid, c("group", "assay", "position"), sep="~") %>% 
                spread(cid, mean, fill=0) %>% 
                remove_rownames() %>%
                column_to_rownames(var="index")
            rclust = kmeans(rr, centers=k[[j]])
            
            #hierarchical clustering of k-means centers
            centerclust = rclust$centers %>% dist() %>% hclust() %>% dendsort(isReverse=TRUE)
            
            reorder = rclust$cluster %>% as_tibble() %>%
                rownames_to_column(var="index") %>% 
                mutate(index = ordered(index, levels=levels(dflist[[1]]$index)),
                       value = ordered(value, levels=centerclust$order)) %>% 
                arrange(value, index)
            roworder = bind_rows(reorder, roworder)
        }
        for (i in 1:nassays){
            dflist[[i]]$index = ordered(dflist[[i]]$index, levels=roworder$index)
        }
    }
    
    plotlist = list()
    
    theme_default = theme_minimal() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=12, color="black", face="bold", margin = unit(c(0,0,0,0), "cm")),
              axis.title.y = element_blank(),
              strip.text = element_text(size=12, color="black", face="bold"),
              strip.text.y = element_blank(),
              strip.background = element_blank(),
              legend.position="top",
              legend.text = element_text(size=10, face="plain"),
              legend.margin = margin(0,0,0,0),
              legend.box.margin = margin(0,0,0,0),
              panel.grid.major.x = element_line(color="black"),
              panel.grid.minor.x = element_line(color="black"),
              panel.grid.major.y = element_blank(),
              panel.grid.minor.y = element_blank(),
              panel.spacing.x = unit(.25, "cm"))
    
    for (i in 1:nassays){
        plotlist[[assays[[i]]]] =
            ggplot(data = dflist[[i]], aes(x=position, y=fct_rev(index), fill=mean)) +
            geom_raster() +
            scale_fill_viridis(option='inferno',
                               na.value = "#FFFFFF00",
                               guide=guide_colorbar(title.position="top", barwidth=12,
                                                    barheight=1, title.hjust=0.5),
                               name=if (logtxn[[i]]){ bquote(bold(log[2] ~ .(assays[[i]]) ~ "signal"))}
                                    else {bquote(bold(.(assays[[i]]) ~ "signal"))}) +
            scale_y_discrete(expand=c(0.02, 0)) +
            facet_grid(annotation~group, scale="free_y", space="free_y", switch="y") +
            theme_default
        if (type=="absolute"){
            plotlist[[i]] = plotlist[[i]] +
                scale_x_continuous(breaks = scales::pretty_breaks(n=3),
                                   labels = format_xaxis(refptlabel = refptlabel,
                                                         upstream = upstream,
                                                         dnstream = dnstream),
                                   name= paste("distance from", refptlabel, 
                                               if_else(upstream>500 | dnstream>500, "(kb)", "(nt)")),
                                   expand=c(0.05, 0))
        }
        else {
            plotlist[[i]] = plotlist[[i]] +
                scale_x_continuous(breaks = c(0, (scaled_length/2)/1000, scaled_length/1000),
                                   labels = c(refptlabel, "", endlabel),
                                   name= "scaled distance",
                                   expand=c(0.05, 0))
        }
    }
    
    #EXTREMELY JANKY WAY TO GET FACET LABELS...build an invisible plot lmao...
    facet_label = ggplot(data = dflist[[1]], aes(x=0, y=fct_rev(index), fill=mean)) +
        geom_raster() +
        scale_fill_gradient(low="#FFFFFF00", high="#FFFFFF00",
                            guide=guide_colorbar(title.position="top", barwidth=0.1,
                                                barheight=1, title.hjust=0.5),
                           name=bquote(bold("."))) +
        scale_y_discrete(expand=c(0.02, 0)) +
        facet_grid(annotation~group, scale="free_y", space="free_y") +
        theme_default +
        theme(text = element_text(color="#FFFFFF00"),
              strip.text.x = element_text(color="#FFFFFF00"),
              axis.text.x = element_text(color="#FFFFFF00"),
              panel.grid = element_blank(),
              strip.text.y = element_text(size=12, color="black", face="bold", angle=0, hjust=1))
    
    for (i in 0:((nassays-1) %/% 4)){
        plotlist = append(plotlist, list(facet_label), after=4*i+i)
    }
   
    allplots = plot_grid(plotlist = plotlist, align="h", ncol=min(nassays+1, 5), axis="trbl")
    
    ggplot2::ggsave(outpath, plot = allplots,
           width=2+2/15*max(nchar(as.character(unique(dflist[[1]]$annotation))))+min(nassays,4)*12,
           height=ceiling(nassays/4)*20, units="cm", limitsize=FALSE)
}

main(inputs = args$inputs,
     cutoffs = args$cutoffs,
     logtxn = args$logtxn,
     assays = args$assays,
     type = args$type,
     refptlabel = paste(args$refptlabel, collapse=" "),
     upstream = args$upstream,
     dnstream= args$dnstream,
     scaled_length= args$scaled_length,
     endlabel = paste(args$endlabel, collapse=" "),
     cluster = args$cluster,
     cluster_assays = args$cluster_assays,
     k = args$k,
     outpath = args$output)
