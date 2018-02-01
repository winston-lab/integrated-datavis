library(argparse)
library(psych)
library(tidyverse)
library(forcats)
library(ggthemes)

parser = ArgumentParser()
parser$add_argument('-i', '--inputs', type='character', nargs='+')
parser$add_argument('-a', '--assays', type='character', nargs='+')
parser$add_argument('-p', '--trim_pct', type='double')
parser$add_argument('-t', '--type', type='character')
parser$add_argument('-r', '--refptlabel', type='character', nargs='+')
parser$add_argument('-u', '--upstream', type='integer')
parser$add_argument('-d', '--dnstream', type='integer')
parser$add_argument('-s', '--scaled_length', type='integer')
parser$add_argument('-e', '--endlabel', type='character', nargs='+')
parser$add_argument('--out1', type='character')
parser$add_argument('--out2', type='character')
parser$add_argument('--out3', type='character')
parser$add_argument('--out4', type='character')

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

main = function(inputs, assays, trim_pct, type, refptlabel,
                upstream, dnstream, scaled_length, endlabel,
                out1, out2, out3, out4) {

    theme_default = theme_light() +
        theme(text = element_text(size=12, color="black", face="bold"),
              axis.text = element_text(size=12, color="black"),
              axis.text.y = element_text(size=10, face="plain"),
              axis.title = element_text(size=10, face="plain"),
              strip.text = element_text(size=12, color="black", face="bold"),
              strip.text.y = element_text(angle=-180, hjust=1),
              strip.background = element_blank(),
              strip.placement="outside",
              legend.position="top",
              legend.title = element_blank(),
              panel.spacing.x = unit(0.5, "cm"),
              plot.margin = margin(.25, .5, .25, .25, "cm"))
    
    x_label = function(ggp){
        if (type=="absolute"){
            ggp = ggp +
                scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                                   label = format_xaxis(refptlabel = refptlabel,
                                                        upstream=upstream,
                                                        dnstream=dnstream),
                                   name = paste("distance from", refptlabel,
                                                if_else(upstream>500 | dnstream>500, "(kb)", "(nt)")),
                                   expand=c(0,0))
        }
        else {
            ggp = ggp +
                scale_x_continuous(breaks = c(0, (scaled_length/2)/1000, scaled_length/1000),
                                   labels = c(refptlabel, "", endlabel),
                                   name= "scaled distance",
                                   expand=c(0,0))
        }
        return(ggp)
    }
    
    format_meta = function(ggp){
        ggp = ggp + geom_vline(xintercept=0, size=1, color="grey65")
        if (type=="scaled"){
            ggp = ggp +
                geom_vline(xintercept=scaled_length/1000, size=1, color="grey65")
        }
        ggp = ggp +
            geom_ribbon(alpha=0.4, size=0) +
            geom_line() +
            scale_fill_ptol(guide=guide_legend(label.position="top", label.hjust=0.5)) +
            scale_color_ptol() +
            scale_y_continuous(name = "normalized signal", breaks=c(0,0.5,1)) +
            theme_default
        ggp = x_label(ggp)
        return(ggp)
    }
    
    sampledf = tibble()
    groupdf = tibble()

    nassays = length(assays)
    
    for (i in 1:nassays){
        tempdf = read_tsv(inputs[[i]],
                          col_names=c('group', 'sample', 'annotation', 'assay',
                                      'index', 'position', 'signal')) %>% 
            mutate_at(vars(group, sample, annotation), funs(fct_inorder(., ordered=TRUE)))

        #put the number of each annotation into the annotation name
        ff = tempdf %>% group_by(annotation) %>% summarise(nanno=n_distinct(index))
        tempdf = tempdf %>% left_join(ff, by="annotation") %>%
            mutate(annotation = fct_inorder(paste(nanno, annotation), ordered=TRUE)) %>% 
            select(-nanno)
        
        tempsample = tempdf %>% group_by(group, sample, annotation, assay, position) %>% 
            summarise(mean = winsor.mean(signal, trim=trim_pct),
                      sem = winsor.sd(signal, trim=trim_pct)/sqrt(n())) %>% 
            ungroup() %>%
            mutate(sem = sem/(max(mean, na.rm=TRUE)-min(mean, na.rm=TRUE)),
                   mean = (mean-min(mean, na.rm=TRUE))/(max(mean, na.rm=TRUE)-min(mean, na.rm=TRUE)))
        sampledf = sampledf %>% bind_rows(tempsample)
        rm(tempsample)
        
        tempgroup = tempdf %>% group_by(group, annotation, assay, position) %>%
            summarise(mean = winsor.mean(signal, trim=trim_pct),
                      sem = winsor.sd(signal, trim=trim_pct)/sqrt(n())) %>% 
            ungroup() %>%
            mutate(sem = sem/(max(mean, na.rm=TRUE)-min(mean, na.rm=TRUE)),
                   mean = (mean-min(mean, na.rm=TRUE))/(max(mean, na.rm=TRUE)-min(mean, na.rm=TRUE)))
        groupdf = groupdf %>% bind_rows(tempgroup)
        rm(tempdf, tempgroup)
    }
    sampledf = sampledf %>% mutate(assay = fct_inorder(assay, ordered=TRUE))
    groupdf = groupdf %>% mutate(assay = fct_inorder(assay, ordered=TRUE))
        
    nannotations = n_distinct(sampledf$annotation)
    ngroups = n_distinct(sampledf$group)
    facet_anno_width = 2+2*max(nchar(unique(as.character(sampledf$annotation))))/15+6*nassays
    facet_group_width = 2+2*max(nchar(unique(as.character(sampledf$group))))/15+6*nassays
    
    sample_facet_anno = ggplot(data = sampledf, aes(x=position, y=mean, ymax=mean+1.96*sem, ymin=mean-1.96*sem,
                                                    group=sample, color=group, fill=group)) +
        facet_grid(annotation~assay, switch="y")
    sample_facet_anno = format_meta(sample_facet_anno)
    ggsave(out1, plot=sample_facet_anno, width=facet_anno_width, height=4+4*nannotations, units="cm", limitsize=FALSE)
    
    sample_facet_group = ggplot(data = sampledf %>% mutate(zz=fct_inorder(paste(sample, annotation), ordered=TRUE)), aes(x=position, y=mean, ymax=mean+1.96*sem, ymin=mean-1.96*sem,
                                group=zz, color=annotation, fill=annotation)) +
        facet_grid(group~assay, switch="y")
    sample_facet_group = format_meta(sample_facet_group)
    ggsave(out2, plot=sample_facet_group, width=facet_group_width, height=4+4*ngroups, units="cm", limitsize=FALSE)
    
    group_facet_anno = ggplot(data = groupdf, aes(x=position, y=mean, ymax=mean+1.96*sem, ymin=mean-1.96*sem,
                                                  color=group, fill=group)) +
        facet_grid(annotation~assay, switch="y")
    group_facet_anno = format_meta(group_facet_anno)
    ggsave(out3, plot=group_facet_anno, width=facet_anno_width, height=4+4*nannotations, units="cm", limitsize=FALSE)
    
    group_facet_group = ggplot(data = groupdf, aes(x=position, y=mean, ymax=mean+1.96*sem, ymin=mean-1.96*sem,
                                group=annotation, color=annotation, fill=annotation)) +
        facet_grid(group~assay, switch="y")
    group_facet_group = format_meta(group_facet_group)
    ggsave(out4, plot=group_facet_group, width=facet_group_width, height=4+4*ngroups, units="cm", limitsize=FALSE)
}

main(inputs = args$inputs,
     assays = args$assays,
     trim_pct = args$trim_pct,
     type = args$type,
     refptlabel = paste(args$refptlabel, collapse=" "),
     upstream = args$upstream,
     dnstream= args$dnstream,
     scaled_length= args$scaled_length,
     endlabel = paste(args$endlabel, collapse=" "),
     out1 = args$out1,
     out2 = args$out2,
     out3 = args$out3,
     out4 = args$out4)
