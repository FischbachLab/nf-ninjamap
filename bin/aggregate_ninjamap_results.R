#!/usr/bin/env Rscript

#knitr::opts_chunk$set(echo = TRUE)
#install.packages("tidyverse")
#devtools::install_github("tidyverse/tidyverse")
library(tidyverse)
library(fs)
library(ggbeeswarm)
library(lubridate)
library(zip)
library(foreach)
library(waffle)
library(ggplot2)
source("/work/aggregate_data_util_functions.R")
pdf(NULL)

args <- commandArgs(trailingOnly = TRUE)

s3_nm_base_path  <- args[1]  # input s3 path e.g., 's3://genomics-workflow-core/Results/Ninjamap'
db_name <-  args[2]  # db name, e.g., hCom2
study_name <- args[3] # project name, e.g., Jeff-20250304
#study_raw_output_path <- args[4]

time_stamp = now() %>% format("%Y%m%d")

local_path = "aggregated_ninjamap_results"
if (! dir.exists(local_path)){
  dir.create(local_path, recursive = TRUE, mode = "0777")
}
#mkdir("aggregated_ninjamap_results")
workdir = './aggregated_ninjamap_results' # work dir
#print(paste0("Current workdir: ", workdir))


local_base = paste(workdir,db_name, sep='/')
s3_base = paste0(s3_nm_base_path,"/",db_name)
study_s3_base = paste(s3_base,study_name, sep='/')
#study_output_base = make_dir(paste(local_base,study_name, sep = '/'))

study_output_base = workdir
study_analysis_path = make_dir(paste(study_output_base,"analysis", sep = '/'))
study_raw_output_path = make_dir(paste(study_output_base,"raw_output", sep = '/'))
study_figures_path = make_dir(paste(study_output_base,"figures", sep = '/'))

# download_all_runs(runs_list,study_raw_output_path,study_s3_base)
download_from_s3(study_raw_output_path, study_s3_base)

strain_dropouts = suppressMessages(aggregate_abundance_data(parent_path = study_raw_output_path))
strain_read_stats = suppressMessages(aggregate_read_stats(parent_path = study_raw_output_path))

host_contamination_stats = suppressMessages(aggregate_host_stats(parent_path = study_raw_output_path))

zipr(paste0(study_raw_output_path,".zip"), study_raw_output_path)
dir_delete(study_raw_output_path)

sd_read_fraction_matrix = strain_dropouts %>%
  # Remove rows with NA values 
  na.omit() %>%
  select(sample_id, Strain_Name, Read_Fraction) %>%
  group_by(Strain_Name) %>%
  #ungroup %>%
  #summarise(n = n()) 
  spread(sample_id, Read_Fraction)

sd_percent_cov_matrix = strain_dropouts %>%
  na.omit() %>%
  select(sample_id, Strain_Name, Percent_Coverage) %>%
  group_by(Strain_Name) %>%
  spread(sample_id, Percent_Coverage)

sd_cov_depth_matrix = strain_dropouts %>%
  na.omit() %>%
  select(sample_id, Strain_Name, Coverage_Depth) %>%
  group_by(Strain_Name) %>%
  spread(sample_id, Coverage_Depth)

strain_dropouts %>%
  write_csv(paste0(study_analysis_path,'/',time_stamp,"_",study_name,".long.csv"),append = FALSE,col_names = TRUE)
strain_read_stats %>%
  write_csv(paste0(study_analysis_path,'/',time_stamp,"_",study_name,".read_stats.csv"),append = FALSE,col_names = TRUE)

sd_read_fraction_matrix %>%
  write_csv(paste0(study_analysis_path,'/',time_stamp,"_",study_name,".readFraction.csv"),append = FALSE,col_names = TRUE)
sd_percent_cov_matrix %>%
  write_csv(paste0(study_analysis_path,'/',time_stamp,"_",study_name,".percCoverage.csv"),append = FALSE,col_names = TRUE)
sd_cov_depth_matrix %>%
  write_csv(paste0(study_analysis_path,'/',time_stamp,"_",study_name,".covDepth.csv"),append = FALSE,col_names = TRUE)

host_contamination_stats %>%
  write_csv(paste0(study_analysis_path,'/',time_stamp,"_",study_name,".host_contaminants.csv"),append = FALSE,col_names = TRUE)



# db_meta %>%
#   anti_join(., sd_cov_depth_matrix, by = c("bin_name" = "Strain_Name")) %>%
#   select(bin_name)


fragments_stats_df = aggregate_fragment_stats(strain_read_stats, how="all")

waffle_df = strain_read_stats %>%
  process_fragment_stats() %>%
  summarize_stats_by_steps()

sum(waffle_df$values)

p <- waffle_df %>%
  waffle_plot(title = "NinjaMap Read Fate", subtitle = "Median across all samples")
# save plot to file without using ggsave
png(paste(study_figures_path,"read_fates.waffle_plot.png", sep="/"))
print(p)
dev.off()
# old eps file that MAC preview doesn't support any longer
#ggsave( paste(study_figures_path,"read_fates.waffle_plot.eps", sep="/"), device = "eps",width = 11,height = 8.5,units = "in", dpi = "retina")

foreach(sample_name=unique(strain_read_stats$sample_id)) %do% {
  p = NULL
  print(sample_name)
  # exclude specific samples
  if (sample_name != "PCR-control-1_S28" && sample_name != "PCR-control-2_S32") {
   
    p = strain_read_stats %>%
      filter(sample_id == sample_name) %>%
      process_fragment_stats() %>%
      summarize_stats_by_steps() %>%
      waffle_plot(title = "NinjaMap Read Fate", subtitle = paste0("Sample: ",sample_name))
   
    # save plot to file without using ggsave
    png(paste0(study_figures_path,"/",sample_name,".read_fates.waffle_plot.png"))
    print(p)
    dev.off()
    #ggsave(plot = p,
     #      filename = paste0(study_figures_path,"/",sample_name,".read_fates.waffle_plot.eps"),
     #      device = "eps",
     #      width = 11,
     #      height = 8.5,
     #      units = "in",
    #       dpi = "retina")
  }
}

plot_fragment_stats(fragments_stats_df)
ggsave(paste(study_figures_path,"fragment_stats.pdf", sep="/"), height = 12, width = 16, units = "in",dpi = "retina")

aggregate_fragment_stats(strain_read_stats, how="post_qc") %>%
  plot_fragment_stats()
ggsave(paste(study_figures_path,"fragment_stats.post_qc.pdf", sep="/"), height = 12, width = 16, units = "in",dpi = "retina")

strain_dropouts %>%
  mutate(formatted_name = sample_id) %>%
  brian_plot(title = "Strain abundance by Sample")+
  geom_point(aes(color=Strain_Name), show.legend = FALSE)+
  geom_line(aes(group=Strain_Name), alpha=0.1,show.legend = FALSE)
ggsave(paste(study_figures_path,"brian_plot.basic.pdf", sep="/"), height = 12, width = 16, units = "in",dpi = "retina")

# host_contamination_stats header
#sample_id,Host_Contaminant,Total_mapped_paired_reads,"Mapped_Paired_rate(%)  select(starts_with('Mapped_Paired_rate'))

# Rename the column header
names(host_contamination_stats)[4] <- "Mapped_Paired_rate"

host_contamination_stats %>%
ggplot()+
      geom_col(aes(sample_id, Total_mapped_paired_reads, fill=Host_Contaminant), position="dodge") +
      geom_line(aes(x=sample_id, y=Mapped_Paired_rate), color="blue", group = 1) +
      scale_y_continuous(sec.axis = sec_axis(~./10000,labels = scales::percent, name ="Percentage"))+
  labs(title="Host Contamination Detection", x="Samples", y= "Contamination Redas")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

host_contamination_stats %>%
ggplot()+
      geom_col(aes(sample_id, Mapped_Paired_rate, fill=Host_Contaminant), position="dodge") +
  labs(title="Host Contamination Detection- Sample Mapping Rate(Filtered Reads)", x="Samples", y= "Contamination Rate")+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 

ggsave(paste(study_figures_path,"Host_Contamination_stats.pdf", sep="/"), height = 8, width = 12, units = "in",dpi = "retina")

# Eliminate the generation of "Rplots.pdf"
dev.off()
