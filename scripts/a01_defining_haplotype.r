#### defining BRCA2 rs81002858 haplotype
#### variant (build 38): rs81002858    13:32326497:A:G    c.517-2A>G
#### variant (build 37): 13:32900634
library(ggplot2)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)

source("scripts/haplotype_analysis_functions.r")

#### run in bash - needs to be run only once
# cd ${phased_gt_data_path}
# cut -d ' ' -f 1-5 cohort.gsa.chip.b37_chr13_HRC.r1-1_phased.haps > ${brca2_haplotype_project_path}/data/cohort_gsa_chr13_phased_SNPs.txt
#### \run in bash
# haps data frame contains actual genotypes; for each individual there are two haplotypes
haps = fread(paste0(phased_gt_data_path, "cohort.gsa.chip.b37_chr13_HRC.r1-1_phased.haps"))
# list of samples, in the same order as in the haps dataframe
samples = fread(paste0(phased_gt_data_path, "cohort.gsa.chip.b37_chr13_HRC.r1-1_phased.sample")

d = fread("data/cohort_gsa_chr13_phased_SNPs.txt")

### Creating a list of SNPs that are in the region surrounding the variant
# my SNP is at 32900634 (b37; 13:32326497 on b38); looking around that position (250kb in this case, but take as big window is needed to get ~100 SNPs)
# i get 106 SNPs in this way
d_snp = d %>% filter((V3 > (32900634-250000)) & (V3 < (32900634+250000)))

### extracting these SNPs from the haplotype file; this df contains 106 rows (SNPs) and 4193 columns (two DNA strands per individual)
haps_brca2 = inner_join(haps, d_snp, by = c("V1", "V2", "V3", "V4", "V5"))
### add column names
## first 5 columns are easy 
names(haps_brca2)[1:5] = c("chr", "variant", "pos", "A1", "A2")

## the remaining columns are sample IDs
# keep the order of IDs from the original file
haps_ids = data.frame(ID = samples$ID_2[2:nrow(samples)])
haps_ids$ord = 1:nrow(haps_ids)
# since for every person i have 2 haplotypes in .haps file I need to duplicate ids
haps_ids = rbind(haps_ids, haps_ids) %>% arrange(ord)
haps_ids$ID = paste(haps_ids$ID, 1:2, sep = "_")

## rename columns in the haps data to contain sample IDs
names(haps_brca2)[6:ncol(haps_brca2)] = haps_ids$ID

write.table(haps_brca2, 
  "data/cohort_gsa_brca2_haps.haps", 
  dec = ".", sep = " ", quote = FALSE, row.names = FALSE)


### Visualising
haps_brca2_l = haps_brca2 %>%
  mutate(snpOrd = 1:nrow(haps_brca2)) %>%
  gather(hapl, value, contains("_"))

pdf("results/cohort_gsa_brca2_hapl_visual.pdf")
  ggplot(haps_brca2_l, aes(x = hapl, y = variant)) + 
    geom_tile(aes(fill = value)) + 
    theme(legend.position = "top", legend.direction = "horizontal")
dev.off()


##### visualising only affected individuals. These were found in the ACMG analysis. The file consist only of carrier IDs, one in each row. Column name is id
carriers = fread("data/brca2.variant_13.32326497.A.G.txt")
# it does not matter that this .fam file is from b38; I only need sample IDs, which are the same on both builds
gsa_ids = fread(paste0(gt_data_path, "/cohort.gsa.chip.b38.fam"))
gsa_ids = gsa_ids %>% rename(id = V2)

gsa_carr = inner_join(carriers, gsa_ids, by = "id") # 9 samples

### extracting these people from the haplotype file
# creating the id column that is the same as the hapl column, but with _1 or _2 removed
# COHORT needs to be replaced with whatever the root of the sample IDs is in your cohort. 
# eg if samples in the cohort are called ABCD0001, ABCD0002, et, COHORT would be replaced with ABCD
haps_aff_l = haps_brca2_l %>%
		mutate(id = gsub("(COHORT\\d+)_\\d", "\\1", hapl))
# mergeing with the carrier list - this way we'll only have the carriers in haps_aff_l
haps_aff_l = right_join(haps_aff_l, select(gsa_carr, id), by = "id")
haps_aff_l = haps_aff_l %>% select(-id)

## visualise them
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

# From these plots one can detect which haplotypes are matching
pdf("results/cohort_gsa_brca2_hapl_affected_initial_visual.pdf")
  ggplot(haps_aff_l, aes(x = hapl, y = snpOrd)) + 
    geom_tile(aes(fill = value)) + 
    #geom_text(aes(label = seq), colour = "white", size = 2) +
    #facet_wrap(~id, ncol = 3, scales = "free") + 
    theme(legend.position = "none", axis.text.x = element_text(angle = -90))
dev.off()

png("results/cohort_gsa_brca2_hapl_affected_initial_visual.png", res = 150, width = 1000, height = 1000)
  ggplot(haps_aff_l, aes(x = hapl, y = snpOrd)) + 
    geom_tile(aes(fill = value)) + 
    #geom_text(aes(label = seq), colour = "white", size = 2) +
     #facet_wrap(~id, ncol = 3, scales = "free") + 
    theme(legend.position = "none", axis.text.x = element_text(angle = -90))
dev.off()

##### Searching for the affected haplotype in carriers
haps_aff_l = haps_aff_l %>% select(-id)
haps_aff_l = haps_aff_l %>% select(-seq)

### this haplotype was found by visual inspection of the above plots. 
# this data frame has the following columns: variant, pos, aff_seq. aff_seq is the affected haplotype sequence, in allelic notation (eg A A A G C etc)
# aff_haplotype_id is ID of the affected haplotype - in the form of ID1234_2
gsa_aff_seq = define_affected_seq(haps_aff_l, aff_haplotype_id)
# pasting it all together in one string
gsa_aff_seq_seq = paste(gsa_aff_seq$aff_seq, collapse = "")

### now searching for this exact same sequence in the carriers, to find which exact haplotypes carry the variant
gsa_hapl_carriers = haplotype_search(haps_aff_l, gsa_aff_seq_seq)

write.table(gsa_hapl_carriers, 
            "results/cohort_gsa_brca2_carrier_haplotype_ids.txt", 
            sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)

# saving also in the format needed for the next script (just haplotype IDs)  - it's used in the a03_haplotype_length_summary script
gsa_hapl_carriers_out = gsa_hapl_carriers %>% select(hapl)
write.table(gsa_hapl_carriers_out, 
  "brca2_gsa_carrier_haplotypes.txt", 
  sep = "\t", dec = ".", quote = FALSE, row.names = FALSE, col.names = FALSE)

# visualise these
haps_aff = haps_aff_l %>% select(-snpOrd) %>% spread(hapl, value)
hapl_carriers_vis = prepare_haplotype_visualisation(haps_aff, gsa_hapl_carriers)
# leaving only the carrier haplotype in
hapl_carriers_vis = left_join(select(gsa_hapl_carriers, hapl), hapl_carriers_vis, by = "hapl")
pdf("results/cohort_gsa_brca2_hapl_carrier_haplotype.pdf")
  haplotype_visualisation(hapl_carriers_vis)
dev.off()

png("results/cohort_gsa_brca2_hapl_carrier_haplotype.png", res = 150, width = 1000, height = 1000)
  haplotype_visualisation(hapl_carriers_vis)
dev.off()


