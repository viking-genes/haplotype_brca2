#### the comparison of whether the carriers from UKB share the same haplotype as each other and as our cohort
rm(list = ls(all = TRUE))
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

source("scripts/haplotype_analysis_functions.r")

# The carriers file has two columns, separated by space. Each row contains ID of the carrier, eg 123456 123456 (the first and the second column are the same). They were found in the exome sequencing data
carriers = fread("data/brca2.variant_13.32900634.A.G_ukbb_200k_carriers_forQCTOOLS.txt")
haps = fread("data/ukbb_brca2.variant_13.32900634.A.G_carrier.hap.haps", na.strings = "-")
samples = fread("data/ukbb_brca2.variant_13.32900634.A.G_carrier.hap.sample")

#### creating the list of SNPs in the file
#### You can use whichever data/tool to create this file
#### the output file needs to be tab separated, without a header and contain columns chr rsid position reference_allele alternate_allele
## done in bash - needs to be run only once ##########################
#cd /data

#cut -d ' ' -f 1-5 ukbb_brca2.variant_13.32900634.A.G_carrier.hap.haps > ukbb_chr13_phased_SNPs.txt
## \bash #################################

d = fread("data/ukbb_chr13_phased_SNPs.txt")

### include sample IDs in column names
# keep the order of IDs from the original file
haps_ids = data.frame(ID = samples$ID_2[2:nrow(samples)])
haps_ids$ord = 1:nrow(haps_ids)
# since for every person i have 2 haplotypes in .haps file I need to duplicate ids
haps_ids = rbind(haps_ids, haps_ids) %>% arrange(ord)
haps_ids$ID = paste(haps_ids$ID, 1:2, sep = "_")

## rename columns in the data
names(haps) = c("chr", "variant", "pos", "A1", "A2", haps_ids$ID)
### Creating a list of SNPs that are in the region surrounding the variant
# my SNP is at 32900634 (b37; 13:32326497 on b38); looking around that position (250kb in this case, but take as big window is needed to get ~100 SNPs).
# actually, in this case I need fewer SNPs as the sharing might be quite short; looking at getting ~50 SNPs. I also know from the initial viking analysis that
# the haplotype is not symmetrical around the rare variant, so taking more SNPs upstream. If this region is too big, the haplotype_search function won't work 
# as the region will already include breaking of the haplotype
# i get 53 SNPs in this way
d_snp = d %>% filter((V3 > (32900634-5000)) & (V3 < (32900634+200000)))
names(d_snp) = c("chr", "variant", "pos", "A1", "A2")

haps_brca2 = inner_join(haps, d_snp, by = c("chr", "variant", "pos", "A1", "A2"))
# remove missing values
haps_brca2 = haps_brca2[complete.cases(haps_brca2),]

write.table(haps_brca2, 
  "data/ukbb_brca2_haps.haps", 
  dec = ".", sep = " ", quote = FALSE, row.names = FALSE)

### Visualising
# adding our variant in
variant = data.frame(chr = 13, variant = "rs81002858", pos = 32900634, A1 = "A", A2 = "G")
haps_brca2 = full_join(haps_brca2, variant, by = c("chr", "variant", "pos", "A1", "A2"))

haps_brca2_l = haps_brca2 %>%
  arrange(pos) %>%
  mutate(snpOrd = 1:nrow(haps_brca2)) %>%
  gather(hapl, value, contains("_"))

pdf("results/ukbb_brca2_hapl_visual.pdf", width = 10)
  ggplot(haps_brca2_l, aes(x = hapl, y = snpOrd)) + 
    geom_tile(aes(fill = value)) + 
    theme(legend.position = "top", legend.direction = "horizontal")
dev.off()

png("results/ukbb_brca2_hapl_visual.png", width = 1000)
  ggplot(haps_brca2_l, aes(x = hapl, y = snpOrd)) + 
    geom_tile(aes(fill = value)) + 
    theme(legend.position = "top", legend.direction = "horizontal")
dev.off()

### looking into the plot above I can tell that these haplotypes are likely shared: 1153382_2, 1958217_2, 3839424_1, 4395333_1, 4917538_2
ukb_aff_seq = define_affected_seq(haps_brca2_l, "1153382_2")
# pasting it all together in one string
ukb_aff_seq_seq = paste(ukb_aff_seq$aff_seq, collapse = "")

### now searching for this exact same sequence in the 9 carriers, to find which exact haplotypes carry the variant
ukb_hapl_carriers = haplotype_search(haps_brca2_l, ukb_aff_seq_seq)

write.table(ukb_hapl_carriers, 
            "results/ukb_brca2_carrier_haplotype_ids.txt", 
            sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)

### visualising only affected individuals 
ukb_hapl_carriers = ukb_hapl_carriers %>% mutate(id = as.character(id))

haps_aff_l = haps_brca2 %>%
  arrange(pos) %>%
  mutate(snpOrd = 1:nrow(haps_brca2)) %>%
  gather(hapl, value, contains("_"))
haps_aff_l = haps_aff_l %>% mutate(id = gsub("(.*)_.*", "\\1", hapl))

haps_aff_l = right_join(haps_aff_l, select(ukb_hapl_carriers, id), by = "id")
haps_aff_l = haps_aff_l %>% select(-id)

### visualise them
# make sure that snps are in the right order
haps_aff_l = haps_aff_l %>% 
  spread(hapl, value)
haps_aff_l = haps_aff_l %>% 
  arrange(pos) %>% 
  mutate(snpOrd = 1:nrow(haps_aff_l))
haps_aff_l = haps_aff_l %>% 
  gather(hapl, value, contains("_")) %>%
  mutate(id = gsub("(\\d+)_\\d", "\\1", hapl)) %>%
  mutate(seq = ifelse(value == 0, A1, A2))

pdf("results/ukbb_200k_brca2_hapl_affected_visual.pdf")
  ggplot(haps_aff_l, aes(x = hapl, y = snpOrd)) + 
    geom_tile(aes(fill = value)) + 
    geom_text(aes(label = seq), colour = "white", size = 2) +
    #facet_wrap(~id, ncol = 3, scales = "free") + 
    theme(legend.position = "none", axis.text.x = element_text(angle = -90))
dev.off()



### haplotype length

### i want to run everyone against everyone, so I need exact haplotype id + just a list of all other ids, and that for everyone
### i need a data.frame that will contain three columns: hapl, target, candidate. For each target haplotype I need all other haplotypes and their IDs

carriers = ukb_hapl_carriers %>% select(hapl, id)
carr_exp = ddply(carriers, .(hapl), function(d) {
	d = data.frame(target = d$hapl, candidate = carriers$id)
	# making sure to not analyse target vs itself
	target_id = gsub("(.*)_.*", "\\1", unique(d$target))
	d = d %>% filter(candidate != target_id)
})
# because column names can't start with numbers R has added X infront of every ID, which then brakes the hapl_def_fine function
# as there's a gsub there that doesn't replace the whole column name and then cnames = Xcandidate, which doesn't work
# to avoid this I'm adding X in front of every ID
carr_exp = carr_exp %>% mutate(target = paste0("X", target), candidate = paste0("X", candidate))

haps = haps %>% arrange(pos)
haps = haps[complete.cases(haps), ]
haps = haps %>% mutate(ord = 1:nrow(haps))

shared_haps = ddply(carr_exp, .(target, candidate), function(d) {
	haplotype_definition_fine(32900634, data.frame(haps), d$target, d$candidate)
})

fwrite(shared_haps, 
	"results/shared_haplotypes/ukbb_brca2_shared_haplotypes.txt",
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)

# calculate lengths
ukb_length = ddply(shared_haps, .(target, candidate), function(a) {
	data.frame(length_mb = (max(a$pos, na.rm = T) - min(a$pos, na.rm = T))/1000000)
})

fwrite(ukb_length, 
	"results/ukb_brca2_shared_haplotypes_lengths.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)


###### visualise the whole region that's share
d_snp = d %>% filter((V3 >  min(shared_haps$pos) & (V3 < max(shared_haps$pos))))
names(d_snp) = c("chr", "variant", "pos", "A1", "A2")

haps_brca2_region = inner_join(haps, d_snp, by = c("chr", "variant", "pos", "A1", "A2"))
# selecting only carrier haplotypes
haps_brca2_region = haps_brca2_region %>% select(chr, variant, pos, A1, A2, any_of(ukb_hapl_carriers$hapl))

# remove missing values
haps_brca2_region = haps_brca2_region[complete.cases(haps_brca2_region),]

## Visualising
# adding our variant in
variant = data.frame(chr = 13, variant = "rs81002858", pos = 32900634, A1 = "A", A2 = "G")
haps_brca2_region = full_join(haps_brca2_region, variant, by = c("chr", "variant", "pos", "A1", "A2"))

haps_brca2_region_l = haps_brca2_region %>%
  arrange(pos) %>%
  mutate(snpOrd = 1:nrow(haps_brca2_region)) %>%
  gather(hapl, value, contains("_"))

pdf("results/ukbb_brca2_hapl_visual_whole_region.pdf", width = 10)
  ggplot(haps_brca2_region_l, aes(x = hapl, y = snpOrd)) + 
    geom_tile(aes(fill = value)) + 
    theme(legend.position = "top", legend.direction = "horizontal")
dev.off()

###### Do UKB people share haplotype with COHORT people?
### take the smallest shared haplotype from UKB, take just this region in the COHORT, merge the two and have a look
ukb_shortest_pair = ukb_length %>% slice(which.min(length_mb))
# take their shared haplotype and min and max pos
ukb_shortest = shared_haps %>% filter(target == ukb_shortest_pair$target, candidate == ukb_shortest_pair$candidate)

## extract this region from the ukb haplotypes data
ukb_snps = d %>% filter((V3 >= min(ukb_shortest$pos)) & (V3 <= max(ukb_shortest$pos)))
names(ukb_snps) = c("chr", "variant", "pos", "A1", "A2")

ukb_haps_region = inner_join(haps, ukb_snps, by = c("chr", "variant", "pos", "A1", "A2"))
# extract only 1 haplotype, from the shortest haplotype pair
ukb_haps_region = ukb_haps_region %>% select(chr, variant, pos, A1, A2, any_of(gsub("X", "", ukb_shortest_pair$target)))
# convert 0|1 into sequence
# it turns out that in the UKB (or COHORT) A1 and A2 are the other way around compared to COHORT so need to swap them around
ukb_haps_region = ukb_haps_region %>% mutate(new_A1 = A2, new_A2 = A1) %>% select(-A1, -A2)
ukb_haps_region = ukb_haps_region %>% rename(A1 = new_A1, A2 = new_A2)
ukb_haps_region = ukb_haps_region %>% mutate(ukb_seq = ifelse(`1153382_2` == 0, A1, A2))

## now take this region out of the cohort data
cohort_snps = fread("data/cohort_gsa_chr13_phased_SNPs.txt")
cohort_region = cohort_snps %>% filter((V3 >= min(ukb_shortest$pos)) & (V3 <= max(ukb_shortest$pos)))

## prepare cohort haplotypes
cohort_haps = fread(paste0(hased_gt_data_path, "cohort.gsa.chip.b37.regeneron_chr13_HRC.r1-1_phased.haps"))
# list of samples, in the same order as in the haps dataframe
cohort_samples = fread(paste0(hased_gt_data_path, "cohort.gsa.chip.b37.regeneron_chr13_HRC.r1-1_phased.sample"))

cohort_haps_region = inner_join(cohort_haps, cohort_region, by = c("V1", "V2", "V3", "V4", "V5"))
## add column names
# first 5 columns are easy 
names(cohort_haps_region)[1:5] = c("chr", "variant", "pos", "A1", "A2")

## the remaining columns are sample IDs
# keep the order of IDs from the original file
cohort_haps_ids = data.frame(ID = cohort_samples$ID_2[2:nrow(cohort_samples)])
cohort_haps_ids$ord = 1:nrow(cohort_haps_ids)
# since for every person i have 2 haplotypes in .haps file I need to duplicate ids
cohort_haps_ids = rbind(cohort_haps_ids, cohort_haps_ids) %>% arrange(ord)
cohort_haps_ids$ID = paste(cohort_haps_ids$ID, 1:2, sep = "_")

## rename columns in the haps data to contain sample IDs
names(cohort_haps_region)[6:ncol(cohort_haps_region)] = cohort_haps_ids$ID

## now extract selected COHORT id for comparison
# aff_haplotype_id is ID of the affected haplotype - in the form of ID1234_2 = so instead of aff_haplotype_id write ID1234 (no quotes needed)
cohort_haps_region = cohort_haps_region %>% select(chr, variant, pos, A1, A2, aff_haplotype_id)
# convert 0|1 into sequence
cohort_haps_region = cohort_haps_region %>% mutate(cohort_seq = ifelse(aff_haplotype_id == 0, A1, A2))


### merge the two datasets together
d_merged = full_join(cohort_haps_region, ukb_haps_region, by = c("chr", "pos", "A1", "A2"), suffix = c(".cohort", ".ukb"))
d_merged = d_merged %>% arrange(pos)
d_merged$snpOrd = 1:nrow(d_merged)
d_merged_l = d_merged %>% gather(id, seq, contains("_seq"))

pdf("results/ukbb_cohort_hapl_comparison_visual.pdf")
  ggplot(d_merged_l, aes(x = id, y = snpOrd)) + 
    geom_tile(aes(fill = seq)) + 
    theme(legend.position = "top", legend.direction = "horizontal") + 
    theme_bw()
dev.off()

pdf("results/ukbb_cohort_hapl_comparison_visual_complete_cases.pdf")
  ggplot(d_merged_l[which(complete.cases(d_merged_l)),], aes(x = id, y = snpOrd)) + 
    geom_tile(aes(fill = seq)) + 
    theme(legend.position = "top", legend.direction = "horizontal") + 
    theme_bw()
dev.off()
	
### check if the other ID, the one that didn't match within the UKB matches cohort
ukb_haps_region_odd = inner_join(haps, ukb_snps, by = c("chr", "variant", "pos", "A1", "A2"))
# extract only 1 haplotype, from the shortest haplotype pair
ukb_haps_region_odd = ukb_haps_region_odd %>% select(chr, variant, pos, A1, A2, `4438134_1`, `4438134_2`, `1153382_2`)
# convert 0|1 into sequence
# it turns out that in the UKB (or cohort) A1 and A2 are the other way around compared to cohort so need to swap them around
ukb_haps_region_odd = ukb_haps_region_odd %>% mutate(new_A1 = A2, new_A2 = A1) %>% select(-A1, -A2)
ukb_haps_region_odd = ukb_haps_region_odd %>% rename(A1 = new_A1, A2 = new_A2)
ukb_haps_region_odd = ukb_haps_region_odd %>% mutate(`4438134_1_seq` = ifelse(`4438134_1` == 0, A1, A2))
ukb_haps_region_odd = ukb_haps_region_odd %>% mutate(`4438134_2_seq` = ifelse(`4438134_2` == 0, A1, A2))
ukb_haps_region_odd = ukb_haps_region_odd %>% mutate(`1153382_2_seq` = ifelse(`1153382_2` == 0, A1, A2))

d_merged_odd = full_join(cohort_haps_region, ukb_haps_region_odd, by = c("chr", "pos", "A1", "A2"), suffix = c(".cohort", ".ukb"))
d_merged_odd = d_merged_odd %>% arrange(pos)
d_merged_odd$snpOrd = 1:nrow(d_merged_odd)
d_merged_odd_l = d_merged_odd %>% gather(id, seq, contains("_seq"))

pdf("results/ukbb_odd_cohort_hapl_comparison_visual.pdf")
  ggplot(d_merged_odd_l, aes(x = id, y = snpOrd)) + 
    geom_tile(aes(fill = seq)) + 
    theme(legend.position = "top", legend.direction = "horizontal") + 
    theme_bw()
dev.off()

pdf("results/ukbb_odd_cohort_hapl_comparison_visual_complete_cases.pdf")
  ggplot(d_merged_odd_l[complete.cases(d_merged_odd),], aes(x = id, y = snpOrd)) + 
    geom_tile(aes(fill = seq)) + 
    theme(legend.position = "top", legend.direction = "horizontal") + 
    theme_bw()
dev.off()

# aff_haplotype_id is ID of the affected haplotype - in the form of ID1234_2 = so instead of aff_haplotype_id write ID1234 (no quotes needed)
# ukb_haplotype_id is the same for the UKB
d_merged_odd_binary_l = d_merged_odd %>% gather(id, value, ukb_haplotype_id, aff_haplotype_id)

pdf("results/ukbb_odd_cohort_hapl_comparison_visual_binary.pdf")
  ggplot(d_merged_odd_binary_l, aes(x = id, y = snpOrd)) + 
    geom_tile(aes(fill = value)) + 
    theme(legend.position = "top", legend.direction = "horizontal") + 
    theme_bw()
dev.off()

pdf("results/ukbb_odd_cohort_hapl_comparison_visual_binary_complete_cases.pdf")
  ggplot(d_merged_odd_binary_l[complete.cases(d_merged_odd_binary_l),], aes(x = id, y = snpOrd)) + 
    geom_tile(aes(fill = value)) + 
    theme(legend.position = "top", legend.direction = "horizontal") + 
    theme_bw()
dev.off()
