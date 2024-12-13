# BAMM analysis
# https://www.biorxiv.org/content/10.1101/2023.09.26.559479v1.full
# https://github.com/padpadpadpad/myxo_diversification/tree/main

# Analysis of the BAMM (Bayesian Analysis of Macroevolutionary Mixtures) run looking for variation in speciation, extinction, and diversification rates on the phylogenetic tree.
# follows lots of the work here: http://bamm-project.org/postprocess.html#bammtools


#--------------------------#
# what this script does ####
#--------------------------#

# plots the tree with most likely rate shift configuration
# runs a phylogenetic regression on tip diversification rate across floral architecture

#------------------------------#
# load in packages and data ####
#------------------------------#

# load packages
# librarian::shelf(caper, ggtree, ggnewscale, RColorBrewer, patchwork, ape, 
#                  phytools, BAMMtools, coda, MetBrewer, nlme, emmeans, tidyverse, 
#                  magick)

library(BAMMtools)
library(diversitree)
library(tidyverse)
library(phytools)
library(treeio) 
library(ggtree)
library(dplyr)
library(tidyr)
library(ape)
library(phangorn)
library(RColorBrewer)
library(MetBrewer)
library(ggnewscale)
library(castor)
library(ggtreeExtra)
library(magick)
library(plotrix)
library(patchwork) # https://patchwork.data-imaginist.com/reference/plot_layout.html
library(openxlsx)

library(hisse)
library(gridExtra)
library(viridis)
library(ggplot2)

#----------------------------------#
# Create folder to save results ####
#----------------------------------#

dir <- "results_BAMM"
todaydate <- format(Sys.time(), "%d%b%Y")
folder_name <- paste0(dir, "/", todaydate)
if (!dir.exists(dir)) {
  dir.create(dir)
}
if (!dir.exists(folder_name)) {
  dir.create(folder_name)
}


#--------------#
# Load data ####
#--------------#

# Note that I have corrected the subfamily classification for "Micklethwaitia"
# in all input files
list.files("data")
tree <- read.tree("data/legume_matK.BEAST_treePL.tre")
keel.misse.data = read.csv("data/keel_misse_4plot.csv", row.names = 1)

# load in bamm event data
edata <- getEventData(tree, eventdata = "data/event_data.txt", burnin = 0.4)


#----------------------------------#
# Look at number of rate shifts ####
#----------------------------------#

# get the same thing from the event data
shift_probs <- summary(edata)

# plot the probabilities
plot_rateshifts <- ggplot(shift_probs, aes(shifts, prob)) +
  geom_col(col = 'black', fill = 'light grey') +
  theme_bw(base_size = 10) +
  labs(x = 'Number of shifts',
       y = 'Probability')


#---------------------------------#
# Plot results from bamm model ####
#---------------------------------#

# calculate credible shift set - the distinct set of shift configurations that account for 95% of the probability of the data
d_css <- credibleShiftSet(edata, expectedNumberOfShifts = 500, threshold = 5, set.limit = 0.95)

# number of distinct configurations in the data
d_css$number.distinct

# view more information about the credible set
summary(d_css)
# even the single best shift configuration has a very very low posterior probablity

# extract the shift configuration that has the maximum marginal probability and plot
# calculate max shift credibility
msc_set <- maximumShiftCredibility(edata, maximize='product')

# grab the best configuration and plot it
msc_config <- subsetEventData(edata, index = msc_set$sampleindex)
plot.bammdata(msc_config, lwd=2)
addBAMMshifts(msc_config, cex = 2)


# phylorate plot
pdf(paste0(folder_name, '/phylorate.pdf'), width = 10, height = 10)
plot.bammdata(edata, lwd = 2,  tau = 0.001, breaksmethod = 'jenks', method = 'polar')
addBAMMshifts(edata, par.reset = FALSE, cex = 1)
#addBAMMlegend(q, location=c(0, 1, 140, 220))
dev.off()



# plot this configuration using ggtree

# get mean phylorates that underly the colorised plot produced by plot.bammdata
# from here: https://groups.google.com/g/bamm-project/c/W6s38xzm6OU/m/LALF47xVS54J
#mbt <- getMeanBranchLengthTree(edata, rate = "speciation")

# get the mean branch lengths from the best tree configuration as identified from maximumShiftCredibility 
mbt <- getMeanBranchLengthTree(msc_config, rate = "ndr")

# get shift nodes from "best model"
shiftnodes <- getShiftNodesFromIndex(edata, index = msc_set$sampleindex)

# get tree
tree_bamm <- mbt$phy

tree_bamm$tip.label[grepl("CAESAL", tree_bamm$tip.label)]

# get the edge lengths in a dataframe
d_tree_bamm <- data.frame(tree_bamm$edge, edge_num=1:nrow(tree_bamm$edge), edge_length = tree_bamm$edge.length)
colnames(d_tree_bamm)=c("parent", "node", "edge_num", 'edge_length')

# transform these to log for the colour scale
d_tree_bamm <- mutate(d_tree_bamm, log_edge_length = log(edge_length))

# constrained families
constrained_groups <- c('POLYGALACEAE', 'PAPILIONOIDEAE', 'CERCIDOIDEAE', 
                        'DETARIOIDEAE', 'CAESALPINIOIDEAE')

# reorder d_meta so the tip labels link to the order of the tips in the tree
d_meta <- tibble(tip_label = tree_bamm$tip.label,
                 keel = keel.misse.data[tree$tip.label, 1],
                 misse = keel.misse.data[tree$tip.label, 2],
                 group = NA)
d_meta$group[grepl("POLYGALACEAE|Polygalaceae", d_meta$tip_label)] <- "POLYGALACEAE"
d_meta$group[grepl("PAPILIONOIDEAE", d_meta$tip_label)] <- "PAPILIONOIDEAE"
d_meta$group[grepl("CERCIDOIDEAE", d_meta$tip_label)] <- "CERCIDOIDEAE"
d_meta$group[grepl("DETARIOIDEAE", d_meta$tip_label)] <- "DETARIOIDEAE"
d_meta$group[grepl("CAESALPINIOIDEAE", d_meta$tip_label)] <- "CAESALPINIOIDEAE"

# find the mrca of each of the constrained families
d_meta2 <- filter(d_meta, group %in% constrained_groups) %>%
  dplyr::select(group, tip_label) %>%
  group_by(group) %>%
  nest() %>%
  mutate(mrca = NA)

for(i in 1:nrow(d_meta2)){
  d_meta2$mrca[i] <- findMRCA(tree_bamm, tips = d_meta2$data[[i]]$tip_label)
}

d_meta2 <- dplyr::select(d_meta2, group2 = group, mrca) %>%
  mutate(blank_label = '')

# add colour for the different groups
cols <- c(colorRampPalette(brewer.pal(11, "Spectral"))(nrow(d_meta2)))
names(cols) <- sort(d_meta2$group2)

# remove any tip labels not in our tree from d_meta
d_meta <- filter(d_meta, tip_label %in% tree_bamm$tip.label)

d_meta$keel <- as.character(d_meta$keel)

th = max(branching.times(tree))
rootAge = max(node.depth.edgelength(tree))

palaeocene <- th - c(66, 56)
oligocene <- th - c(33.9, 23)
pliocene <- th - c(5.333, 2.58)
segments <- data.frame(
  xmin = c(palaeocene[2], oligocene[2], pliocene[2]),
  xmax = c(palaeocene[1], oligocene[1], pliocene[1]),
  ymin = -Inf,
  ymax = Inf
)

# plot tree using ggtree
# first colour branches and add rate shifts

p1 <- ggtree(tree, 
             layout = 'fan', 
             open.angle = 10,
             aes(col = log_edge_length)) %<+% 
  d_tree_bamm +
  annotate("rect", 
           xmin = segments$xmin, xmax = segments$xmax, 
           ymin = segments$ymin, ymax = segments$ymax, 
           alpha = 0.2, fill = "#A0E0F9") +
  scale_color_gradientn('Net diversification\n(branch colours)', 
                        colors = met.brewer(name='Hiroshige', direction=-1, override.order = F), 
                        breaks=c(min(d_tree_bamm$log_edge_length, na.rm = TRUE) + 
                                   abs(min(d_tree_bamm$edge_length, na.rm = TRUE))*0.2, 
                                 max(d_tree_bamm$log_edge_length, na.rm = TRUE) * 0.95), 
                        labels=c("Slow","Fast")) +
  geom_point2(aes(subset=(node %in% shiftnodes)), color="black", size=3, alpha = 0.8)+
  NULL

# next add tip points

p2 <- p1 %<+% d_meta +
  new_scale_color() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=keel),
    width=2,
    offset=0.02
  ) +
  #geom_tippoint(aes(x=x+x*0.02, col = keel), size = 0.1) +
  scale_fill_manual('Floral architecture\n(tip points)', 
                    values = c("#E6E7E8", "#7F62CE"), 
                    labels = c('non-keel', 'keel')) +

  new_scale_fill() +
  geom_fruit(
    geom=geom_tile,
    mapping=aes(fill=misse),
    width=2,
    offset=0.022
  ) +
  scale_fill_viridis_c(option="inferno", name="MiSSE")
  
  # guides(fill = guide_legend(override.aes = list(size = 5)),
  #        size = 'none')

tree_plot <- p2 +
  new_scale_color() +
  geom_cladelab(data = d_meta2,
                mapping = aes(node = mrca,
                              color = group2,
                              label = blank_label),
                offset = castor::get_all_distances_to_root(tree) %>% max() * 0.06,
                barsize = 4) +
  scale_color_manual('Group (outer bar)', values = cols) +
  guides(color = guide_legend(override.aes = list(size = 0.1, shape = 1)))


tree_plot

# save plot out
ggsave(paste0(folder_name, '/bamm_tree_with_clade_colors.pdf'), tree_plot, height = 9, width = 12)


tree_plot <- p2 +
  new_scale_color() +
  geom_cladelab(data = d_meta2,
                mapping = aes(node = mrca,
                              label = blank_label),
                barcolour = "black",
                textcolour = "white",
                offset = castor::get_all_distances_to_root(tree) %>% max() * 0.06,
                barsize = 4)
tree_plot

# save plot out
ggsave(paste0(folder_name, '/bamm_tree.pdf'), tree_plot, height = 9, width = 12)



#----------------------------------------#
# look at rate variation through time ####
#----------------------------------------#

# write function to get rate through time into the correct format
get_rate_through_time_df <- function(ephy, ...){
  rtt <- getRateThroughTimeMatrix(ephy, ...)
  
  # get dataframe of each part
  # first speciation
  rtt_sp <- rtt$lambda %>%
    data.frame() %>%
    mutate(sample = 1:n()) %>%
    pivot_longer(starts_with('X'), names_to = 'time_point', values_to = 'speciation') %>%
    mutate(time_point = parse_number(time_point)) %>%
    group_by(sample) %>%
    mutate(time = unname(rtt$times)) %>%
    ungroup()
  
  # second extinction
  rtt_ex <- rtt$mu %>%
    data.frame() %>%
    mutate(sample = 1:n()) %>%
    pivot_longer(starts_with('X'), names_to = 'time_point', values_to = 'extinction') %>%
    mutate(time_point = parse_number(time_point)) %>%
    group_by(sample) %>%
    mutate(time = unname(rtt$times)) %>%
    ungroup()
  
  rtt_comb <- left_join(rtt_sp, rtt_ex)
  
  rtt_comb <- mutate(rtt_comb, net_div = speciation - extinction) %>%
    pivot_longer(cols = c(speciation, extinction, net_div), names_to = 'process', values_to = 'rate')
  
  return(rtt_comb)
}

# get rate through time estimates for the whole tree
rtt_all <- get_rate_through_time_df(ephy = edata)

# create means
rtt_combine_means <- group_by(rtt_all, time_point, process, time) %>%
  summarise(ave_rate = mean(rate), .groups = 'drop',
            lower_ci = quantile(rate, 0.025),
            upper_ci = quantile(rate, 0.975))

# get rate through time estimates for each shift node
rtt_shift <- tibble(shift_node = shiftnodes,
                    n = 1:length(shiftnodes)) %>%
  nest(data = shift_node) %>%
  mutate(temp = purrr::map(data, possibly(~get_rate_through_time_df(ephy = edata, node = .x$shift_node, nodetype = 'include'), otherwise = NA_real_)),
         is_tib = purrr::map_dbl(temp, is_tibble))

rtt_shift2 <- filter(rtt_shift, is_tib == 1) %>%
  dplyr::select(data, temp) %>%
  unnest(data) %>%
  unnest(temp)

# create means
rtt_shift_means <- group_by(rtt_shift2, time_point, process, time, shift_node) %>%
  summarise(ave_rate = mean(rate), .groups = 'drop',
            lower_ci = quantile(rate, 0.025),
            upper_ci = quantile(rate, 0.975))

# create a plot
ggplot() +
  geom_line(aes(time, ave_rate, group = shift_node), rtt_shift_means, col = 'dark grey') +
  geom_ribbon(aes(time, ymin = lower_ci, ymax = upper_ci, group = shift_node), alpha = 0.01, rtt_shift_means, fill = 'dark grey') +
  geom_line(aes(time, ave_rate), rtt_combine_means) +
  geom_ribbon(aes(time, ymin = lower_ci, ymax = upper_ci), alpha = 0.1, rtt_combine_means) +
  facet_wrap(~process, scales = 'free') +
  theme_bw(base_size = 14) +
  labs(x = 'Relative time',
       y = 'rate')

# save plot out
ggsave(paste0(folder_name, '/bamm_rate_through_time.png'), last_plot(), height = 5, width = 12)

# create a plot for just net diversification
p_rtt <- ggplot() +
  geom_line(aes(time, ave_rate, group = shift_node), 
            filter(rtt_shift_means, process == 'net_div'), col = 'dark grey') +
  geom_ribbon(aes(time, ymin = lower_ci, ymax = upper_ci, group = shift_node), 
              alpha = 0.05, filter(rtt_shift_means, process == 'net_div'), 
              fill = 'dark grey') +
  geom_ribbon(aes(time, ymin = lower_ci, ymax = upper_ci), alpha = 0.5, 
              filter(rtt_combine_means, process == 'net_div')) +
  geom_line(aes(time, ave_rate), filter(rtt_combine_means, process == 'net_div')) +
  theme_bw(base_size = 10) +
  labs(x = 'Relative time',
       y = 'Net diversification rate')

#--------------------------------------------#
# look at tip-specific evolutionary rates ####
#--------------------------------------------#
# grab out tip rates 
tip_rates <- data.frame(tip_label2 = edata$tip.label,
                        speciation = edata$meanTipLambda,
                        extinction = edata$meanTipMu) %>%
  add_column(keel = d_meta$keel) %>%
  group_by(keel) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(floral_axis = gsub(':', '/ ', keel),
         floral_axis = gsub('_', ' ', keel),
         net_diversification = speciation - extinction)

labels <- gsub("-", "\\\u00ad", c('non-keel', 'keel'))

p1 <- ggplot(tip_rates, aes(forcats::fct_reorder(floral_axis, n), speciation)) +
  geom_boxplot(aes(fill = keel), show.legend = FALSE) +
  geom_point(shape = 21, aes(fill = keel), position = position_jitter(width = 0.2), 
             size = 0.8, stroke = 0.1, alpha = 0.8, show.legend = FALSE) +
  theme_bw(base_size = 12) +
  scale_x_discrete(labels = labels) +
  labs(x = 'Floral architecture',
       y = 'Speciation rate') +
  scale_fill_manual('Floral architecture', values = c("#E6E7E8", "#7F62CE"), labels = labels) 


p2 <- ggplot(tip_rates, aes(forcats::fct_reorder(floral_axis, n), extinction)) +
  geom_boxplot(aes(fill = keel), show.legend = FALSE) +
  geom_point(shape = 21, aes(fill = keel), position = position_jitter(width = 0.2), 
             size = 0.8, stroke = 0.1, alpha = 0.8, show.legend = FALSE) +
  theme_bw(base_size = 12) +
  scale_x_discrete(labels = labels) +
  labs(x = 'Floral architecture',
       y = 'Extinction rate') +
  scale_fill_manual('Floral architecture', values = c("#E6E7E8", "#7F62CE"), labels = labels) 


p3 <- ggplot(tip_rates, aes(forcats::fct_reorder(floral_axis, n), net_diversification)) +
  geom_boxplot(aes(fill = keel), show.legend = FALSE) +
  geom_point(shape = 21, aes(fill = keel), position = position_jitter(width = 0.2), 
             size = 0.8, stroke = 0.1, alpha = 0.8, show.legend = FALSE) +
  theme_bw(base_size = 12) +
  scale_x_discrete(labels = labels) +
  labs(x = 'Floral architecture',
       y = 'Tip specific\nnet diversification rate') +
  scale_fill_manual('Floral architecture', values = c("#E6E7E8", "#7F62CE"), labels = labels) 


p1 + p2 + p3

# save plot out
ggsave(paste0(folder_name, '/bamm_tip_rates.png'), last_plot(), height = 5, width = 17)

# make plot of just net diversification rate
p_tiprates <- tip_rates %>%
  ggplot(., aes(forcats::fct_reorder(floral_axis, n), net_diversification)) +
  geom_boxplot(aes(fill = keel), outliers = FALSE, show.legend = FALSE) +
  geom_point(shape = 21, aes(fill = keel), position = position_jitter(width = 0.2), 
             size = 0.9, stroke = 0.1, alpha = 0.8, show.legend = FALSE) +
  theme_bw(base_size = 10) +
  scale_x_discrete(labels = labels) +
  labs(x = 'Floral architecture',
       y = 'BAMM net diversification rate') +
  scale_fill_manual('Floral architecture', values = c("#E6E7E8", "#7F62CE"), labels = labels) 
  #ylim(c(-2, 30))


p_misse_tiprates <- ggplot(d_meta) +
  geom_boxplot(aes(x=keel, y=misse, fill = keel), outliers = FALSE, show.legend = FALSE) +
  scale_fill_manual('Floral architecture', values = c("#E6E7E8", "#7F62CE"), labels = labels) +

  new_scale_fill() +
  geom_point(shape = 21, aes(x=keel, y=misse, fill = misse), position = position_jitter(width = 0.2), 
             size = 0.9, stroke = 0.1, alpha = 0.8, show.legend = FALSE) +
  scale_fill_viridis_c(option="inferno") +
  
  theme_bw(base_size = 10) +
  scale_x_discrete(labels = labels) +
  labs(x = 'Floral architecture',
       y = 'MiSSE tip rate')



#----------------------#
# Assemble Figure 2 ####
#----------------------#

layout <- c(
  'AAAAAAAA
   CCDDEEFF'
)

tree_plot2 <- image_read(paste0(folder_name, '/bamm_tree_ai.png'), density = 300)
tree_plot2 <- image_trim(tree_plot2)

# make plot
image_ggplot(tree_plot2) + 
  p_rtt + 
  plot_rateshifts + 
  p_tiprates + 
  p_misse_tiprates +
  plot_layout(design = layout, heights = c(0.8, 0.2)) + 
  plot_annotation(tag_levels = list(c('(a)', '(b)', '(c)', '(d)', '(e)'))) &
  theme(plot.tag = element_text(size = 14))

ggsave(paste0(folder_name, '/BAMM_paps_to_add_image.pdf'), last_plot(), width = 12, height = 12)




