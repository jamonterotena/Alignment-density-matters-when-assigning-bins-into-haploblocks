#############################################################################
##### THE ALIGNMENT NUMBER MATTERS WHEN ASSIGNING BINS INTO HAPLOBLOCKS #####
#############################################################################
#The following lines will set up the requirements for running the script. The focus of this script is in the sections 'Feedback' and 'Suggestions'. The functions were developed from your publication to serve further particular purpose, which I will explain in the section that I mentioned.
##### Required packages ##### 
for (pkg in c("plyr", "dplyr", "magrittr", "GenomicRanges", "ggplot2", "tidyr", "viridis", "stringr", "ggdendro", "reshape2", "grid", "imager")) {
  if (!eval(bquote(requireNamespace(.(pkg), quietly = TRUE)))) {
    eval(bquote(install.packages(.(pkg))))
  }
  eval(bquote(library(.(pkg))))
  rm(pkg)
}
##### Required functions ##### 
plot_line_bin_median <- function(data, bin_start = 0, bin_end, bin_size = 1e+07, cut_off = 99.99, ymin = 97, ymax = 100, reference_name = data$rid, query_name = data$qid, x_label_gap){
  comparison_filt_bin <- edited_bin_data(data = data, bin_size = bin_size, bin_start = bin_start, bin_end = bin_end)
  comparison_medians <- aggregate(perc_id ~ bin, data = comparison_filt_bin, FUN=median)
  colnames(comparison_medians)[2] <- "perc_id_median"
  comparison_medians$cut_off <- NA
  comparison_medians$cut_off <- ifelse(comparison_medians$perc_id >= cut_off, paste0(">=",cut_off), paste0("<",cut_off))
  comparison_medians$cut_off <- as.factor(comparison_medians$cut_off)
  
  ggplot(comparison_medians, aes(x=bin, y = perc_id_median, colour = cut_off)) +
    geom_line(colour = "grey") + 
    geom_point(size = 1)  + 
    ylim(ymin,ymax) +
    scale_colour_manual(values = c("73D055FF", "440154FF")) +
    labs(colour = "% id cutoff") +
    # scale_colour_viridis() +
    xlab(paste0 ("alignment position in ", reference_name, "'s chr", unique(data$chrom), " (bp)")) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ylab("Median bin percentage of identity") +
    scale_x_continuous(breaks = round( seq (min(bin_start), max(bin_end), by = x_label_gap), 1), position = "top") +
    
    ggtitle(paste0("Median line: ", reference_name, " vs ", query_name, ", zoom region: ", bin_start, "-", bin_end, ", bin size: ", bin_size/1000000, "-Mbp"))
}
block_summary <- function(median_cutoffs_copy, bin_size, reference_name = "NA", query_name ="NA", show_only_coords = FALSE){
  blocks <- unique(median_cutoffs_copy$block_no[!grepl("NO_BLOCK", median_cutoffs_copy$block_no)])
  blocks <- blocks[complete.cases(blocks)]
  
  block_positions <- data.frame(bin_size = character(), comparison = character(), block_no = numeric(), block_start = numeric(), block_end = numeric())
  
  for(block in blocks){
    bin_size_rep <- paste0(bin_size/1000000, "-Mbp", sep = "")
    block_data <- subset(median_cutoffs_copy, block_no == block)
    block_start <- (min(block_data$bin) - bin_size)
    block_end <- (max(block_data$bin))
    comparison <- paste0(reference_name, "->", query_name)
    to_add <- data.frame(bin_size = bin_size_rep,comparison = comparison, block_no = block, block_start = block_start, block_end = block_end)
    block_positions <- rbind(block_positions, to_add)
  }
  coords <- data.frame("start" = block_positions$block_start, "end" = block_positions$block_end)
  #print(paste0("BLOCK SUMMARY AT ", bin_size, "-MBP BIN SIZE"))
  ifelse(show_only_coords==TRUE, return(coords), return(block_positions))
}
plot_aln_pid_and_length <- function(data, xmin = 0, xmax = max(data$re), ymin = 97, ymax = 100, reference_name = data$rid, query_name = data$qid, x_label_gap = 5000000, dot_size = 2){
  ggplot(data[data$r_mid > xmin & data$r_mid < xmax,], aes(x=r_mid, y=perc_id, colour=r_length)) +
    theme_bw() + 
    xlab(paste0 ("alignment position in ", reference_name, "'s chr", unique(data$chrom), " (bp)")) + 
    ylab(paste0('alignment percentage of identity vs ', query_name)) +
    geom_point(alpha = .5, size = dot_size) +
    ylim(ymin, ymax+0.05) +
    scale_colour_viridis() +
    geom_hline(yintercept = 99.99, linetype = "dashed", color = "red") +
    annotate ("text", colour = "red", size = 4.5, x = xmin, y = 100.025, label = "cutoff 99.99%" ) +
    scale_x_continuous(breaks = seq( xmin, xmax, by = x_label_gap), position = "top") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    ggtitle(paste0("Percentage of identity of individual alignments: ", " ", reference_name, " vs ", query_name, ", zoom region: ", xmin/1e06, "-", xmax/1e06, " Mbp"))
}
edited_bin_data <- function(data, bin_size, bin_start = 0, bin_end){
  bins <- seq(bin_start, bin_end, by = bin_size)
  data$bin <- NA
  for (i in bins){
    data$bin <- ifelse(((data$r_mid > (i-bin_size)) & (data$r_mid < (i-1))), i, data$bin)
  }
  return(data)
}
##### New functions ##### 
plot_aln_and_bins <- function(aln_subset = data.frame(), bin_size = bin_size, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "target region", fill_target = "orange", color_target_text = "black", fill_predictions = "green", color_prediction_text = "black",  ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000, print_tables = TRUE, prediction_text = TRUE, aln_text = TRUE){
  medians_zoom <- plot_line_bin_median(data = aln_subset, bin_size = bin_size, bin_start = zoom_start, bin_end = zoom_end, ymin = ymin, cut_off = cut_off, reference_name = reference_assembly, query_name = query_assembly, x_label_gap = x_label_gap)
  medians_zoom_bin_info <- assign_blocks_mummer(medians_zoom[["data"]], original_file = aln_subset)
  aln_number_per_bin <- medians_zoom_bin_info$aln_number[2:length(medians_zoom_bin_info$aln_number)]
  mid_point <- (medians_zoom_bin_info$bin_start[2:length(medians_zoom_bin_info$aln_number)]+medians_zoom_bin_info$bin_end[2:length(medians_zoom_bin_info$aln_number)])/2
  medians_zoom_block_sum <- block_summary(medians_zoom_bin_info, bin_size = bin_size, show_only_coords = FALSE, reference_name = reference_assembly, query_name = query_assembly)
  if (print_tables) { 
    print(paste0("BINS AT ", bin_size/1000000, "-MBP BIN SIZE"))
    print(medians_zoom_bin_info)
    print(paste0("BLOCK SUMMARY AT ", bin_size/1000000, "-MBP BIN SIZE"))
    print(medians_zoom_block_sum)
  }
  medians_zoom_block_coords <- block_summary(medians_zoom_bin_info, bin_size = bin_size, show_only_coords = TRUE, reference_name = reference_assembly, query_name = query_assembly)
  ifelse(prediction_text, haploblock_prediction_text <- paste0("HAPLOBLOCK PREDICTIONS AT ", bin_size/1e06, "-MBP BIN SIZE"), haploblock_prediction_text <- "")
  ifelse(aln_text, aln_text <- paste0(aln_number_per_bin, " aln"), aln_text <- "")
  coords = data.frame()
  if(nrow(medians_zoom_block_coords) == 0){
    coords <- data.frame(start = zoom_start , end = zoom_end)
    fill_predictions <- "white"
    haploblock_prediction_text <- ""
  } else {
    coords <- block_summary(medians_zoom_bin_info, bin_size = bin_size, show_only_coords = TRUE, reference_name = reference_assembly, query_name = query_assembly)
  }
  graph <- plot_aln_pid_and_length(aln_subset, xmin = zoom_start, ymin = ymin, xmax = zoom_end, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = x_label_gap, dot_size = dot_size)
  if (is.data.frame(highlighted_target) == FALSE){
    graph_common <- graph + geom_vline(xintercept = seq(zoom_start, zoom_end, by = bin_size)) + geom_rect(data = coords, inherit.aes = FALSE, aes(xmin = coords[[1]], xmax = coords[[2]], ymin = ymin+1/4*(100-ymin), ymax = 100), color = "transparent", fill = fill_predictions, alpha = 0.3) +  annotate ("text", x = coords[[1]]+(coords[[2]]-coords[[1]])/2, y = (ymin+1/4*(100-ymin)+2/12*(100-ymin)), label = haploblock_prediction_text, size = 8, fontface = "bold", colour = color_prediction_text) +  annotate("text", x = mid_point, y = (ymin+1/4*(100-ymin)+1/12*(100-ymin)), size = 5, label = aln_text)
    print(graph_common)
  } else {
    graph_onlyiftarget <- graph + geom_vline(xintercept = seq(zoom_start, zoom_end, by = bin_size)) + geom_rect(data = coords, inherit.aes = FALSE, aes(xmin = coords[[1]], xmax = coords[[2]], ymin = ymin+1/4*(100-ymin), ymax = 100), color = "transparent", fill = fill_predictions, alpha = 0.3) +  annotate ("text", x = coords[[1]]+(coords[[2]]-coords[[1]])/2, y = (ymin+1/4*(100-ymin)+2/12*(100-ymin)), label = haploblock_prediction_text, size = 8, fontface = "bold", colour = color_prediction_text) +  annotate("text", x = mid_point, y = (ymin+1/4*(100-ymin)+1/12*(100-ymin)), size = 5, label = aln_text) + geom_rect(data = highlighted_target, inherit.aes = FALSE, aes(xmin = highlighted_target[[1]], xmax = highlighted_target[[2]], ymin = ymin, ymax = (ymin+1/4*(100-ymin))), color = "transparent", fill = fill_target, alpha = 0.3) + annotate ("text", x = highlighted_target[[1]]+(highlighted_target[[2]] - highlighted_target[[1]])/2, y = (ymin+1/6*(100-ymin)), label = target_text, size = 6, fontface = "bold", colour = color_target_text) + annotate ("text", x = highlighted_target[[1]]+(highlighted_target[[2]] - highlighted_target[[1]])/2, y = (ymin+1/12*(100-ymin)), label = paste0(nrow(aln_subset[(aln_subset$r_mid>=highlighted_target[[1]])&(aln_subset$r_mid<=highlighted_target[[2]]),]), " aln"), size = 5) 
    print(graph_onlyiftarget)
  }
}
assign_blocks_mummer <-function(median_cutoffs, original_file){
  median_cutoffs_copy <- median_cutoffs
  median_cutoffs_copy$block_no <- NA
  block_no = 1
  
  for (i in seq(1, nrow(median_cutoffs_copy))){
    if(median_cutoffs_copy[i, "perc_id_median"] < 99.99){
      median_cutoffs_copy[i, "block_no"] <- "NO_BLOCK"
    } else if (median_cutoffs_copy[i, "perc_id_median"] >= 99.99){
      median_cutoffs_copy[i, "block_no"] <- block_no
      if (i > (nrow(median_cutoffs_copy)-3)){
      } else if ((median_cutoffs_copy[i+1, "perc_id_median"] < 99.99) & (median_cutoffs_copy[i+2, "perc_id_median"] < 99.99) & (median_cutoffs_copy[i+3, "perc_id_median"] < 99.99)){
        block_no <- block_no + 1
      }
    }
  }
  binSize <- (median_cutoffs_copy$bin[2] - median_cutoffs_copy$bin[1])
  median_cutoffs_copy$bin_start <- (median_cutoffs_copy$bin - binSize)
  median_cutoffs_copy$bin_end <- (median_cutoffs_copy$bin)
  number_block <- numeric()
  for (i in 1:nrow(median_cutoffs_copy)){
    tempo <- original_file[original_file$r_mid>=as.numeric(median_cutoffs_copy$bin_start[i]) & original_file$r_mid<=as.numeric(median_cutoffs_copy$bin_end[i]), ]
    ntempo <- nrow(tempo)
    number_block <- append(number_block, ntempo)
  }
  median_cutoffs_copy$aln_number <- number_block
  # print(paste0("BINS AT ", bin_size, "-MBP BIN SIZE"))
  # median_cutoffs$bin_size <- rep(binSize, nrow(median_cutoffs))
  return(median_cutoffs_copy)
}
##### New functions considering alignment number #####
assign_blocks_mummer_edited2 <- function (aln_subset = aln_subset, bin_size = bin_size, aln_threshold = aln_threshold, bin_start = zoom_start, bin_end = zoom_end){
  ##### CREATING THE NEW DATA FRATE THAT ASSIGNS MUMMER BLOCKS #####
  median_cutoffs <- plot_line_bin_median(data = aln_subset, bin_start = bin_start, bin_end = bin_end, bin_size = bin_size, 
                                         cut_off = 99.99, ymin = 97, ymax = 100, x_label_gap = 50000000,
                                         reference_name = reference_assembly, query_name = query_assembly)
  median_cutoffs <- median_cutoffs$data
  median_cutoffs$block_no <- NA
  
  ##### ADDING THE NUMBER OF ALIGNMENTS #####
  median_cutoffs$bin_start <- (median_cutoffs$bin - bin_size)
  median_cutoffs$bin_end <- (median_cutoffs$bin)
  aln_number_block <- numeric()
  for (i in seq(1, nrow(median_cutoffs))){
    aln <- aln_subset[aln_subset$r_mid>=as.numeric(median_cutoffs$bin_start[i]) 
                      & aln_subset$r_mid<=as.numeric(median_cutoffs$bin_end[i]), ]
    aln_number <- nrow(aln)
    aln_number_block <- append(aln_number_block, aln_number)
  }
  median_cutoffs$aln_number <- aln_number_block
  
  ##### ADDING THE MUMMER BLOCK NUMBER #####
  block_no = 1
  for (i in seq(1, nrow(median_cutoffs))){
    if ((median_cutoffs[i, "perc_id_median"] < 99.99)){
      median_cutoffs[i, "block_no"] <- "NO_BLOCK"
    } else {
      if (median_cutoffs[i, "aln_number"] >= aln_threshold) {
        median_cutoffs[i, "block_no"] <- block_no
      } else {
        if   (   (i>2 && ((median_cutoffs[i-2, "perc_id_median"] >= 99.99) && (median_cutoffs[i-2, "aln_number"] >= aln_threshold))) 
                 || (i>1 && ((median_cutoffs[i-1, "perc_id_median"] >= 99.99) && (median_cutoffs[i-1, "aln_number"] >= aln_threshold))) 
                 || (i<nrow(median_cutoffs) && ((median_cutoffs[i+1, "perc_id_median"] >= 99.99) && (median_cutoffs[i+1, "aln_number"] >= aln_threshold))) 
                 || (i<(nrow(median_cutoffs)-1) && ((median_cutoffs[i+2, "perc_id_median"] >= 99.99) && (median_cutoffs[i+2, "aln_number"] >= aln_threshold))) ) {
          median_cutoffs[i, "block_no"] <- block_no
        } else {
          median_cutoffs[i, "block_no"] <- "NO_BLOCK"
        }
      }
    } 
    if (i > (nrow(median_cutoffs)-3)){
    } else if  ( ((median_cutoffs[i+1, "perc_id_median"] < 99.99) | (median_cutoffs[i+1, "aln_number"] < aln_threshold)) 
                 & ((median_cutoffs[i+2, "perc_id_median"] < 99.99) | (median_cutoffs[i+2, "aln_number"] < aln_threshold))){
      block_no <- block_no + 1
    }
  }
  return(median_cutoffs)
}
plot_aln_and_bins_edited2 <- function(aln_subset = aln_subset, aln_threshold = aln_threshold, bin_size = bin_size, zoom_start = zoom_start, 
                                      zoom_end = zoom_end, highlighted_target = target, target_text = "target region", fill_target = "orange", 
                                      color_target_text = "black", fill_predictions = "green", color_prediction_text = "black",  ymin = 98.5, 
                                      cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, 
                                      vline = TRUE, x_label_gap = 1000000, print_tables = TRUE, prediction_text = TRUE, aln_text = TRUE){
  #APPLY FUNCTION TO ASSIGN BLOCK
  blocks_dataframe <- assign_blocks_mummer_edited2(aln_subset = aln_subset, bin_size = bin_size, aln_threshold = aln_threshold, bin_start = zoom_start, bin_end = zoom_end)
  #EXTRACT ALIGNMENT NUMBER PER BIN
  aln_number_per_bin <- blocks_dataframe$aln_number[2:length(blocks_dataframe$aln_number)]
  #CALCULATE BIN MIDPOINT
  mid_point <- (blocks_dataframe$bin_start[2:length(blocks_dataframe$aln_number)]+blocks_dataframe$bin_end[2:length(blocks_dataframe$aln_number)])/2
  #CREATE A SUMMARY OF THE BLOCKS
  medians_zoom_block_sum <- block_summary(blocks_dataframe, bin_size = bin_size, show_only_coords = FALSE, reference_name = reference_name, query_name = query_name)
  #RETURN THE SUMMARY OF THE BLOCKS (ONLY IF OPTION print_tables IS TRUE)
  if (print_tables) { 
    print(paste0("BINS AT ", bin_size/1000000, "-MBP BIN SIZE"))
    print(blocks_dataframe)
  }
  medians_zoom_block_coords <- block_summary(blocks_dataframe, bin_size = bin_size, show_only_coords = TRUE, reference_name = reference_name, query_name = query_name)
  ifelse(prediction_text, haploblock_prediction_text <- paste0("HAPLOBLOCK PREDICTIONS AT ", bin_size/1e06, "-MBP BIN SIZE"), haploblock_prediction_text <- "")
  ifelse(aln_text, aln_text <- paste0(aln_number_per_bin, " aln"), aln_text <- "")
  coords = data.frame()
  if(nrow(medians_zoom_block_coords) == 0){
    coords <- data.frame(start = zoom_start , end = zoom_end)
    fill_predictions <- "white"
    haploblock_prediction_text <- ""
  } else {
    coords <- block_summary(blocks_dataframe, bin_size = bin_size, show_only_coords = TRUE, reference_name = reference_name, query_name = query_name)
  }
  graph <- plot_aln_pid_and_length(aln_subset, xmin = zoom_start, ymin = ymin, xmax = zoom_end, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = x_label_gap, dot_size = dot_size)
  if (is.data.frame(highlighted_target) == FALSE){
    graph <- graph + geom_rect(data = coords, inherit.aes = FALSE, aes(xmin = coords[[1]], xmax = coords[[2]], ymin = ymin+1/4*(100-ymin), ymax = 100), color = "transparent", fill = fill_predictions, alpha = 0.3) +  annotate ("text", x = coords[[1]]+(coords[[2]]-coords[[1]])/2, y = (ymin+1/4*(100-ymin)+2/12*(100-ymin)), label = haploblock_prediction_text, size = 8, fontface = "bold", colour = color_prediction_text) +  annotate("text", x = mid_point, y = (ymin+1/4*(100-ymin)+1/12*(100-ymin)), size = 5, label = aln_text)
  } else {
    graph <- graph + geom_rect(data = coords, inherit.aes = FALSE, aes(xmin = coords[[1]], xmax = coords[[2]], ymin = ymin+1/4*(100-ymin), ymax = 100), color = "transparent", fill = fill_predictions, alpha = 0.3) +  annotate ("text", x = coords[[1]]+(coords[[2]]-coords[[1]])/2, y = (ymin+1/4*(100-ymin)+2/12*(100-ymin)), label = haploblock_prediction_text, size = 8, fontface = "bold", colour = color_prediction_text) +  annotate("text", x = mid_point, y = (ymin+1/4*(100-ymin)+1/12*(100-ymin)), size = 5, label = aln_text) + geom_rect(data = highlighted_target, inherit.aes = FALSE, aes(xmin = highlighted_target[[1]], xmax = highlighted_target[[2]], ymin = ymin, ymax = (ymin+1/4*(100-ymin))), color = "transparent", fill = fill_target, alpha = 0.3) + annotate ("text", x = highlighted_target[[1]]+(highlighted_target[[2]] - highlighted_target[[1]])/2, y = (ymin+1/6*(100-ymin)), label = target_text, size = 6, fontface = "bold", colour = color_target_text) + annotate ("text", x = highlighted_target[[1]]+(highlighted_target[[2]] - highlighted_target[[1]])/2, y = (ymin+1/12*(100-ymin)), label = paste0(nrow(aln_subset[(aln_subset$r_mid>=highlighted_target[[1]])&(aln_subset$r_mid<=highlighted_target[[2]]),]), " aln"), size = 5) 
  }
  ifelse(vline, print(graph + geom_vline(xintercept = seq(zoom_start, zoom_end, by = bin_size))), print(graph))
}
##### Parameters ##### 
chromosome <- "5B" 
reference_assembly <- "Lancer" 
query_assembly <- "Paragon" 
zoom_start <- 655000000 
zoom_end <- 665000000 
target_start <- 655700000 
target_end <- 656600000 
target <- data.frame(target_start, target_end)
##### Database ##### 
aln_library <- readRDS(file = paste0("all_20_kb_filtered_delta_", chromosome,"_tables.rds"))
aln_subset <- aln_library[grepl(paste0("^", tolower(reference_assembly), sep = ""), aln_library$comparison) & grepl(tolower(query_assembly), aln_library$comparison),]

##### Feedback ##### 
#Hi! My name is Jose Antonio Montero-Tena and I am a PhD student in bioinformatics in the Department of Plant Breeding of JLU Giessen, Germany. In October 2021, I presented my master thesis, in which I analysed a SNP haplotype involved in higher root dry mass in wheat by using the bioinformatic resources described in your article 'A haplotype-led approach to increase the precision of wheat breeding', focusing mainly on the website crop-haplotypes.com and the R script to assign blocks from mummer alignments. I would like to give you now a brief feedback of the conclusions of my thesis and some suggestions that I thought could help overcome pitfalls originated from the alignment databases.  
#My object of research was the haplotype RDMa-h2. This genome region had been discovered in wheat chr5B and located by BLAST between 655700000 and 656600000. crop-haplotypes.com suggested that the haplotype was shared between Lancer and Paragon. I used the R script to map the region within the chromosome and try to modify the method so that the region could be narrowed down in order to link the associated traits with causal genes or discover recombination breakpoints. However, I observed that the alignment coverage in this database was not as high as desirable when zooming into a 10-Mbp spanning my target haploblock:
#Since my worked focused particularly on a 1-Mbp segment in Chromosome 5B, found to be shared by Lancer and Paragon, I developed some new functions. The particularity of these functions is that they 
plot_aln_pid_and_length(data = aln_subset, xmin = zoom_start, xmax = zoom_end, ymin = 98.5, reference_name = reference_assembly, query_name = query_assembly , x_label_gap = 1000000, dot_size = 4)
#This region contains few alignments and some sections specially limited numbers (~663-664Mbp) that were shown to hamper haploblock assignment:
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 5.0e06, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
#The haploblock is here assigned regardless of this 'gap'. Although this wrong assignment could be overcome by simply reducing the bin size, at 1-Mbp we still face similar issues:
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 1.0e06, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
#These untrustworthy haploblocks are a product of the fact that the method does not take the number of alignments per bin into consideration. Therefore, if some alignments within a bin do not manage to drop the pident low enough below the cutoff, the bin will despite be assigned to a haploblock. This could be the case when some sections of low alignment density are assigned. Also, bins could be assigned even on the base of one alignment.
#In my thesis, I asked myself about the reason for this low alignment coverage and proved that this has to do with the comparison type and the chromosome. Comparisons between chromosome-level assemblies and scaffold-level assemblies have significantly lower alignment coverage, alignment length and number of alignments per bin. 
image1 <- load.image("chr 5B vs ALL.png")
plot(image1, axes=FALSE)
#Also, alignment coverage is not homogeneously distributed between chromosomes, with the chromosomes of my study being specially low on this parameter.
image2 <- load.image("CLA vs SLA.png")
plot(image2, axes=FALSE)
#I would like to suggest you to consider applying modifications to the algorithm so that the haploblocks are not only assigned on the basis of their percentage of identity, but also on their alignment density, since this factor was proved decisive to allow comparisons between chromosome-level assemblies and scaffold-level assemblies. These changes could be done in many ways and I would like to propose some new functions that would include alignment number as a factor for haploblock assignment.   
##### Suggestions ##### 
# If we have a look at the median number of alignments per bin, we get 19, 9 and 4 for bin sizes 5, 2.5 and 1 Mbp
binsize_list <- c(5.0e6, 2.5e6, 1.0e6)
names(binsize_list) <- c("bin size: 5-Mbp", "bin size: 2.5-Mbp", "bin size: 1-Mbp")
for (i in 1:3){
  print("median aln number per bin considering the real distribution")
  blocks_dataframe <- assign_blocks_mummer_edited2(aln_subset = aln_subset, aln_threshold = 18, bin_size = binsize_list[i], bin_start = 0, bin_end = max(aln_subset$re))
  aln_density <- blocks_dataframe$aln_number
  print(median(aln_density))
}
#By considering only percid and bin size 5-Mbp, we get a haploblock at 660-665 Mbp, despite a 'gap' in the section ~663-664Mbp that counts with only 2 low-percid alignments
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 5.0e06, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
#If we now use a threshold of 19, the median number of alignments in 5-Mbp bins, we get no haploblock in the region, since the number of alignments in the second bin is lower than 19
plot_aln_and_bins_edited2(aln_subset = aln_subset, aln_threshold = 19, bin_size = 5.0e06, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "target region", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "blue", color_prediction_text = "darkblue",  ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
#We can repeat the same procedure with the rest two bin sizes that you used in your work
#Bin size 2.5-Mbp, only percid
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 2.5e06, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
#Bin size 2.5-Mbp, percid & aln number
plot_aln_and_bins_edited2(aln_subset = aln_subset, aln_threshold = 9, bin_size = 2.5e06, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "target region", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "blue", color_prediction_text = "darkblue",  ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
#Bin size 1-Mbp, only percid
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 1.0e06, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "SNP haploblock RDMa", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000)
#Bin size 1-Mbp, percid & aln number
plot_aln_and_bins_edited2(aln_subset = aln_subset, aln_threshold = 4, bin_size = 1.0e06, zoom_start = zoom_start, zoom_end = zoom_end, highlighted_target = target, target_text = "target region", fill_target = "orange", color_target_text = "darkorange", fill_predictions = "blue", color_prediction_text = "darkblue",  ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 1000000, prediction_text = FALSE)
#Now, let's have a look across the whole chromosome 5B
#Only percid
plot_aln_and_bins(aln_subset = aln_subset, bin_size = 5.0e06, zoom_start = 0, zoom_end = max(aln_subset$re), highlighted_target = target, target_text = NA, fill_target = "orange", color_target_text = "darkorange", fill_predictions = "green", color_prediction_text = "darkgreen", ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 500000000, print_tables = FALSE, prediction_text = FALSE, aln_text = FALSE)
#percid & aln number
plot_aln_and_bins_edited2(aln_subset = aln_subset, aln_threshold = 19, bin_size = 5.0e06, zoom_start = 0, zoom_end = max(aln_subset$re), highlighted_target = NA, fill_predictions = "blue", color_prediction_text = "darkblue",  ymin = 98.5, cut_off = 99.99, reference_name = reference_assembly, query_name = query_assembly, dot_size = 4, x_label_gap = 500000000, print_tables = FALSE, prediction_text = FALSE, aln_text = FALSE, vline = FALSE)
#Some sections will now longer be assigned because they are not enough covered with alignments.

#Notice that I applied a restrictive criteria by using the median number of alignments. This is because comparisons between chromosome- and scaffold-level assemblies, as this one between Lancer and Paragon, have lower alignment coverage than their counterparts, so probably high thresholds should be used to allow to compared between them. Furthermore, the method could be expanded to, for example, show color ranges depending on the number of alignments across the x axis. The website could benefit from this new algorithm by, for example, giving information about how good alignment density is in each haploblock between pairwise comparisons.
#I would really like to know what you think about my work. I also would like to tell you that I would be very proud to collaborate with you in the future.

#Thanks for your attention!