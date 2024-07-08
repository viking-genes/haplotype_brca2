#### calculating haplotype length between all carriers
rm(list = ls(all = TRUE))
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

args = commandArgs(TRUE) 
target_hapl = args[1]

source("scripts/haplotype_analysis_functions.r")

haps = fread(paste0(phased_gt_data_path, "cohort.gsa.chip.b37_chr13_HRC.r1-1_phased.haps"))
samples = fread(paste0(phased_gt_data_path, "cohort.gsa.chip.b37_chr13_HRC.r1-1_phased.sample"))
carriers = fread("results/cohort_gsa_brca2_carrier_haplotype_ids.txt")

##### defining haplotype more precisely. 
## reformatting hap data a bit
names(haps)[1:5] = c("chr", "variant", "pos", "A1", "A2")
### include sample IDs in column names
# keep the order of IDs from the original file
haps_ids = data.frame(ID = samples$ID_2[2:nrow(samples)])
haps_ids$ord = 1:nrow(haps_ids)
# since for every person i have 2 haplotypes in .haps file I need to duplicate ids
haps_ids = rbind(haps_ids, haps_ids) %>% arrange(ord)
haps_ids$ID = paste(haps_ids$ID, 1:2, sep = "_")
## rename columns in the data
names(haps)[6:ncol(haps)] = haps_ids$ID
haps$ord = 1:nrow(haps)



# i want to run everyone against everyone, so I need exact haplotype id + just a list of all other ids, and that for everyone
carriers = carriers %>% select(hapl, id)
carr_exp = ddply(carriers, .(hapl), function(d) {
	d = data.frame(target = d$hapl, candidate = carriers$id)
	# making sure to not analyse target vs itself
	target_id = gsub("(.*)_.*", "\\1", unique(d$target))
	d = d %>% filter(candidate != target_id)
})
# keeping only the current target
carr_exp = carr_exp %>% filter(target == target_hapl)

shared_haps = ddply(carr_exp, .(target, candidate), function(d) {
	haplotype_definition_fine(32900634, data.frame(haps), d$target, d$candidate)
})

fwrite(shared_haps, 
	sprintf("results/shared_haplotypes/%s_brca2_gsa_shared_haplotypes.txt", target_hapl),
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)


