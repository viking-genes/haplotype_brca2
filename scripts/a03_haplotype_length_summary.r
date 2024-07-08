### summarising the haplotype length data
rm(list = ls(all = TRUE))
library(scales)
library(data.table)
library(plyr)
library(dplyr)
library(tidyr)
library(ggplot2)

source("scripts/haplotype_analysis_functions.r")

haps_path = "results/shared_haplotypes/"
gsa_haps_files = list.files(haps_path, pattern = "gsa_shared_haplotypes")

gsa_haps = lapply(gsa_haps_files, function(f) {
	d = fread(paste0(haps_path, f))
})
gsa_haps = do.call(rbind, gsa_haps)

gsa_length = ddply(gsa_haps, .(target, candidate), function(a) {
	data.frame(length_mb = (max(a$pos) - min(a$pos))/1000000)
})

fwrite(gsa_length, 
	"results/cohort_brca2_gsa_shared_haplotypes_lengths.txt", 
	sep = "\t", dec = ".", quote = FALSE, row.names = FALSE)

##### create a figure of carriers and select max positions.
## first extracting the carriers from the haplotype file of the whole cohort
haps = fread(paste0(phased_gt_data_path, "cohort.gsa.chip.b37_chr13_HRC.r1-1_phased.haps"))
# list of samples, in the same order as in the haps dataframe
samples = fread(paste0(phased_gt_data_path, "cohort.gsa.chip.b37.regeneron_chr13_HRC.r1-1_phased.sample"))

# renaming the columns in haps files
names(haps)[1:5] = c("chr", "variant", "pos", "A1", "A2")
# adding sample names as column names
# keep the order of IDs from the original file
haps_ids = data.frame(ID = samples$ID_2[2:nrow(samples)])
haps_ids$ord = 1:nrow(haps_ids)
# since for every person i have 2 haplotypes in .haps file I need to duplicate ids
haps_ids = rbind(haps_ids, haps_ids) %>% arrange(ord)
haps_ids$ID = paste(haps_ids$ID, 1:2, sep = "_")

## rename columns in the haps data to contain sample IDs
names(haps)[6:ncol(haps)] = haps_ids$ID


### selecting just the SNPs that were shared between any of the haplotype carriers
d = fread("data/cohorts_gsa_chr13_phased_SNPs.txt")
# here I want min and max positions from the gsa_haps
d_snp = d %>% filter((V3 > min(gsa_haps$pos)) & (V3 < max(gsa_haps$pos)))
names(d_snp)[1:5] = c("chr", "variant", "pos", "A1", "A2")

haps_brca2 = left_join(d_snp, haps, by = c("chr", "variant", "pos", "A1", "A2"))

### now selecting the haplotypes
carriers = fread("brca2_gsa_carrier_haplotypes.txt", header = FALSE)
names(carriers) = "hapl"

brca2_carriers = haps_brca2 %>% select(chr, variant, pos, A1, A2, any_of(carriers$hapl))
# adding in our SNP of interest that will just be a grey line
variant = data.frame(chr = 13, variant = "rs81002858", pos = 32900634, A1 = "A", A2 = "G")
brca2_carriers = full_join(brca2_carriers, variant, by = c("chr", "variant", "pos", "A1", "A2"))
brca2_carriers = brca2_carriers %>% arrange(pos)
brca2_carriers$snpOrd = 1:nrow(brca2_carriers)


### reformatting into a long format for visualisation
# cohort_id is root of the ID of samples in the cohort. For example, if samples in the cohort are called COHORT0001, COHORT0002, et, cohort_id would be "COHORT"
brca2_carriers_l = brca2_carriers %>% gather(hapl, value, contains(cohort_id))
# create the sequence
brca2_carriers_l = brca2_carriers_l %>% mutate(seq = ifelse(value == 0, A1, A2))

# creating y-axis breaks for plotting, using snpOrd to plot, but pos to label the ticks
snpOrd_deciles <- round(quantile(brca2_carriers$snpOrd, probs = seq(0, 1, by = 0.1)), digits = 0)
# Extract rows corresponding to deciles
y_labels <- brca2_carriers[brca2_carriers$snpOrd %in% snpOrd_deciles, ]

pdf("results/cohort_gsa_brca2_hapl_carrier_whole_haplotype.pdf")
	ggplot(brca2_carriers_l, aes(x = hapl, y = snpOrd))  + 
        geom_tile(aes(fill = value)) + 
        theme(legend.position = "none", axis.text.x = element_text(angle = -90)) +
        scale_y_continuous(breaks = y_labels$snpOrd, labels = y_labels$pos) +
        labs(y = "chr 13 position")
dev.off()

#### for the figure I want to select a smaller region. I'll take the sharing between target_halp_id and candidate_hapl_id and then few SNPs upstream to show the breaking
# target_halp_id and candidate_hapl_id need to be replace with the exact haplotypes, eg COHORT0001_1 and COHORT0012_2
selected_pos = gsa_haps %>% filter(target == target_hapl_id & candidate == candidate_hapl_id) 
selected_pos = selected_pos %>% slice(c(which.min(pos), which.max(pos)))  

# filtering to keep only positions within the selected range (+8kb upstream to make it to 8Mb region)
brca2_carriers_sel_l = brca2_carriers_l %>% filter(pos >= (min(selected_pos$pos)+1500000) & pos <= (max(selected_pos$pos+2300000 )))
# remaking the labels for y axis
y_brakes_sel = brca2_carriers_sel_l %>% select(snpOrd, pos) %>% distinct()
snpOrd_quintiles <- round(quantile(y_brakes_sel$snpOrd, probs = seq(0, 1, by = 0.25)), digits = 0)
# Extract rows corresponding to deciles
y_labels <- y_brakes_sel[y_brakes_sel$snpOrd %in% snpOrd_quintiles, ]
# divide pos by 1Mb so it's less busy for visualisation
y_labels = y_labels %>% mutate(pos = round(pos/1000000, digits = 1))

pdf("results/cohort_gsa_brca2_hapl_carrier_whole_haplotype_ordered_selected_region.pdf")
    ggplot(brca2_carriers_sel_l, aes(x = hapl, y = snpOrd))  + 
        geom_tile(aes(fill = value)) + 
        scale_y_continuous(breaks = y_labels$snpOrd, labels = y_labels$pos) +
        labs(y = "chr 13 position (Mb)") + 
        theme(legend.position = "none", axis.text.x = element_text(angle = -90, vjust = 0.5), axis.title.x = element_blank()) 
dev.off()

# renaming COHORTID with ID1:9 to retain anonymity
new_ids = data.frame(id = unique(brca2_carriers_sel_l$id), new_id = paste0("ID", 1:9))
brca2_carriers_sel_l = left_join(brca2_carriers_sel_l, new_ids, by = "id")
# keeping the haplotype number
brca2_carriers_sel_l = brca2_carriers_sel_l %>% mutate(hapl_n = gsub(".*_(.*)", "\\1", hapl))
brca2_carriers_sel_l = brca2_carriers_sel_l %>% mutate(new_hapl = paste(new_id, hapl_n, sep = "_"))

### this region contains ~8Mb
pdf("results/cohort_gsa_brca2_hapl_carrier_whole_haplotype_ordered_selected_region_ids_renamed.pdf")
    ggplot(brca2_carriers_sel_l, aes(x = new_hapl, y = snpOrd))  + 
        geom_tile(aes(fill = value)) + 
        scale_y_continuous(breaks = y_labels$snpOrd, labels = y_labels$pos) +
        labs(y = "chr 13 position (Mb)") + 
        scale_fill_gradient(low = "black", high = "white") +  # Set colors for gradient
  		theme_bw() +
        theme(legend.position = "none", axis.text.x = element_text(angle = -90, vjust = 0.5), axis.title.x = element_blank()) 
dev.off()
