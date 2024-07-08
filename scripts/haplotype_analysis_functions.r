

define_affected_seq = function(haps_aff_l, carr_hapl_id) {
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

  haps_aff_seq = haps_aff_l %>%
    select(-id, -value) %>%
    spread(hapl, seq)

  ### defining the sequence of interest
  aff_seq_df = haps_aff_seq %>%
    arrange(snpOrd) %>%
    select(variant, pos, carr_hapl_id) %>%
    rename(aff_seq = carr_hapl_id)
 
  return(aff_seq_df)
}

haplotype_search = function(haps_var_l, aff_seq) {
    ### reformat all haplotypes to sequences
    haps_var_l = haps_var_l %>%
        mutate(seq = ifelse(value == 0, A1, A2))

    ### make every person's haplotype one long string
    haps_var_seq = ddply(haps_var_l, .(hapl), function(d) {
        d = d %>% 
            arrange(snpOrd)
            out = data.frame(hapl_seq = paste(d$seq, collapse = ""))
    })

    ### filter out those that have the same affected haplotype
    haps_aff_all = haps_var_seq %>% 
        filter(hapl_seq == aff_seq) %>% 
        mutate(id = gsub("(\\d+)_\\d", "\\1", hapl))
        
    return(haps_aff_all)

}
### visualising haplotypes
prepare_haplotype_visualisation = function(haps_selected_all, haps_candidates) {
    # adding in the snp of interest
    snp = data.frame(chr = 17, variant = "rs45553935", pos = 41209139)
    haps_selected_all = merge(haps_selected_all, snp, by = c("variant", "chr", "pos"), all = T)
    haps_var_l = haps_selected_all %>%
        mutate(snpOrd = 1:nrow(haps_selected_all)) %>%
        gather(hapl, value, contains("_")) %>% 
        mutate(id = gsub("(\\w+)_\\d", "\\1", hapl))
    
    haps_aff_var = merge(haps_var_l, select(haps_candidates, id), by = "id", all.y = TRUE)

    haps_aff_var = haps_aff_var %>% 
        mutate(seq = ifelse(value == 0, A1, A2)) %>%
        arrange(snpOrd)
    
    return(haps_aff_var)
}

haplotype_visualisation = function(haps_aff_var) {
  p = ggplot(haps_aff_var, aes(x = hapl, y = snpOrd)) + 
        geom_tile(aes(fill = value)) + 
        #geom_text(aes(label = seq), colour = "white", size = 2) +
        #facet_wrap(~id, ncol = 19, scales = "free_x") + 
        theme(legend.position = "none", axis.text.x = element_text(angle = -90))
    
    return(p)
  
}

haplotype_definition_fine = function(var_pos, haps_qt, target_hapl, candidate) {    
    #starting point - just making sure I have at least some variants to start with
    hap_lower  = haps_qt %>% filter((pos > (var_pos-100000) & (pos <= var_pos)))
    lower_start = which(haps_qt$variant == hap_lower$variant[which.max(hap_lower$pos)])
    hap_lower = hap_lower[, c("variant", "pos", "ord", "A1", "A2", grep(target_hapl, names(haps_qt), value = T), grep(candidate, names(haps_qt), value = T))]
    names(hap_lower) = gsub(candidate, "candidate", names(hap_lower))
    names(hap_lower) = gsub(target_hapl, "carrier_hapl", names(hap_lower))
    # keeping only 1 variant for the start
    hap_lower = hap_lower %>% filter(ord == lower_start)
    hap_lower_l = hap_lower %>% gather(id, hapl, contains("_"))
    hap_lower_l = hap_lower_l %>% mutate(seq = ifelse(hapl == 0, A1, A2))
    hap_lower = hap_lower_l %>% select(-hapl) %>% spread(id, seq)   
    hap_lower = hap_lower %>% arrange(ord)
    # calculate the number of SNPs on the chromosome that are below the start position
    lower_chr_tot_n = haps_qt %>% filter(pos <= hap_lower$pos) %>% nrow() 

    count =1
        # run this while either of the candidate haplotypes match with the target haplotype or when it hits the end of a chromosome
        while ((count <= lower_chr_tot_n & (sum(hap_lower$carrier_hapl != hap_lower$candidate_1) == 0) & hap_lower$ord[count] > 1) | 
                ((count <= lower_chr_tot_n & sum(hap_lower$carrier_hapl != hap_lower$candidate_2) == 0) & hap_lower$ord[count] > 1) ) {
            # selecting everyting from the start position to the count position
            hap_lower  = data.frame(haps_qt[(lower_start - count):lower_start, ])
            # selecting the target haplotype and both haplotypes of the candindate
            hap_lower = hap_lower[, c("variant", "pos", "ord", "A1", "A2", grep(target_hapl, names(haps_qt), value = T), grep(candidate, names(haps_qt), value = T))]
            names(hap_lower) = gsub(candidate, "candidate", names(hap_lower))
            names(hap_lower) = gsub(target_hapl, "carrier_hapl", names(hap_lower))
            # transforming into long format so that we can rename the alleles from 0 and 1 to actual alleles
            hap_lower_l = hap_lower %>% gather(id, hapl, contains("_"))
            hap_lower_l = hap_lower_l %>% mutate(seq = ifelse(hapl == 0, A1, A2))
            hap_lower = hap_lower_l %>% select(-hapl) %>% spread(id, seq)   
            hap_lower = hap_lower %>% arrange(ord)
            writeLines(sprintf("%s %s downstream iteration no %s", target_hapl, candidate, count))
            count=count+1
        }
    # the last iteration there was a mismatch and that's why while loop finished, but we still have that last variant that's actually a mismatch 
    # in the table. in this case the last variant will be the most downstream variant, so in the end i want the table minus the first row
    hap_lower = hap_lower %>% arrange(ord)
    hap_lower = hap_lower[2:nrow(hap_lower) ,]
    
    hap_upper  = haps_qt %>% filter((pos < (var_pos+100000) & (pos > var_pos)))
    upper_start = which(haps_qt$variant == hap_upper$variant[which.min(hap_upper$pos)])
    hap_upper = hap_upper[, c("variant", "pos", "ord", "A1", "A2", grep(target_hapl, names(haps_qt), value = T), grep(candidate, names(haps_qt), value = T))]
    names(hap_upper) = gsub(candidate, "candidate", names(hap_upper))
    names(hap_upper) = gsub(target_hapl, "carrier_hapl", names(hap_upper))
    # keeping only 1 variant for the start
    hap_upper = hap_upper %>% filter(ord == upper_start)
    hap_upper_l = hap_upper %>% gather(id, hapl, contains("_"))
    hap_upper_l = hap_upper_l %>% mutate(seq = ifelse(hapl == 0, A1, A2))
    hap_upper = hap_upper_l %>% select(-hapl) %>% spread(id, seq)   
  
    # run this while either of the candidate haplotypes match with the target haplotype or when it hits the end of a chromosome
    count =1
    while (((sum(hap_upper$carrier_hapl != hap_upper$candidate_1) == 0) & hap_upper$ord[count] < nrow(haps_qt)) | 
              ((sum(hap_upper$carrier_hapl != hap_upper$candidate_2) == 0) & hap_upper$ord[count] < nrow(haps_qt))) {
        hap_upper  = data.frame(haps_qt[upper_start:(upper_start + count), ])
        hap_upper = hap_upper[, c("variant", "pos", "ord", "A1", "A2", grep(target_hapl, names(haps_qt), value = T), grep(candidate, names(haps_qt), value = T))]
        names(hap_upper) = gsub(candidate, "candidate", names(hap_upper))
        names(hap_upper) = gsub(target_hapl, "carrier_hapl", names(hap_upper))
        hap_upper_l = hap_upper %>% gather(id, hapl, contains("_"))
        hap_upper_l = hap_upper_l %>% mutate(seq = ifelse(hapl == 0, A1, A2))
        hap_upper = hap_upper_l %>% select(-hapl) %>% spread(id, seq)   
        writeLines(sprintf("%s %s upstream iteration no %s", target_hapl, candidate, count))
        count=count+1
    }
    # the last iteration there was a mismatch and that's why while loop finished, but we still have that last variant that's actually a mismatch 
    # in the table. so in the end i want the table minus the last row
    hap_upper = hap_upper %>% arrange(pos)
    hap_upper = hap_upper[1:(nrow(hap_upper)-1) ,]
    haplotype = rbind(hap_lower, hap_upper)
    return(haplotype)
}   

king_pairs_segments_plotting = function(all_segments, pair_data, title) {
    p = ggplot() +
        geom_rect(data = all_segments, aes(xmin = start, xmax = end, ymin = 0, max = 0.9), fill = 'white', color = "black", size = 0.45) + 
        geom_rect(data = pair_data, aes(xmin = start, xmax = end, ymin = 0, ymax = 0.9, fill = IBDType)) + 
        #geom_rect(data = pair_data, aes(xmin = start, xmax = end, ymin = 0, max = 0.9), color = "black", alpha = 0, size = 0.45) + 
        scale_fill_manual(values = c("IBD0" = "white", "IBD1" = "dodgerblue2", "IBD2" = "firebrick2"), drop = FALSE) + 
        facet_grid(chr ~ .) +
        labs(x = "Position (Mb)", y = "", title = title) +
        theme(
            legend.position = "bottom", legend.key = element_rect(color = "black"),
        panel.background = element_rect(fill = 'grey80', color = 'grey80'), panel.border = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        axis.text.y = element_blank(), axis.ticks.y = element_blank()
    )
    return(p)
}

