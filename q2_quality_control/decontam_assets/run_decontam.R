#!/usr/bin/env Rscript

options(error = function() {
  sink(stderr())
  on.exit(sink(NULL))
  traceback(3)
  if (!interactive()) {
    q(status = 1)
  }
})

library("decontam")
library("optparse")

cat(R.version$version.string, "\n")
errQuit <- function(mesg, status=1) { message("Error: ", mesg); q(status=status) }

option_list = list(
  make_option(c("--asv_table_path"), action="store", default='NULL', type='character',
              help="File path to table in .csv format "),
  make_option(c("--threshold"), action="store", default='NULL', type='character',
              help="threshold value for decontam algorithm"),
  make_option(c("--decon_method"), action="store", default='NULL', type='character',
              help="algorithm mode"),
  make_option(c("--output_track"), action="store", default='NULL', type='character',
              help="File path to tracking tsv file. If already exists, will be overwritten"),
  make_option(c("--meta_table_path"), action="store", default='NULL', type='character',
              help="File path to metadata in .tsv format"),
  make_option(c("--freq_con_column"), action="store", default='NULL', type='character',
              help="Name of column for frequency method"),
  make_option(c("--prev_control_or_exp_sample_column"), action="store", default='NULL', type='character',
              help="Name of column for prevalence method"),
  make_option(c("--prev_control_sample_indicator"), action="store", default='NULL', type='character',
              help="Indicator to identify control samples")
)
opt = parse_args(OptionParser(option_list=option_list))

inp.loc <- opt$asv_table_path
threshold <- if(opt$threshold=='NULL') NULL else as.numeric(opt$threshold)
out.track <- opt$output_track
metadata.loc<-opt$meta_table_path
decon.mode<-opt$decon_method
prev.control.col <- opt$prev_control_or_exp_sample_column
prev.id.controls<-opt$prev_control_sample_indicator
freq.con.col<-opt$freq_con_column

if(!file.exists(inp.loc)) {
  errQuit("Input ASV table does not exist.")
}else if(!file.exists(metadata.loc)) {
  errQuit("Input metadata file does not exist.")
}else{
  print("Input ASV and Metadata files found")
}

meta_data_cols <-function(asv_df, metadata_df, control.col){
  control_vec<-c()
  index<-0
  for (id in colnames(metadata_df)) {
    index=index+1
    if(id == control.col){
      control_vec<-metadata_df[,c(index)]
    }
  }
  # We need to align the metadata to the asv table
  mapped_ids <- match( rownames(asv_df),rownames(metadata_df))
  control_vec <- control_vec[mapped_ids]
  # drop sample IDs not in the asv table
  control_vec <- na.omit(control_vec)
  return(control_vec)
}

outputer<-function(decon_output, out.track){
  cat("Write output\n")

  write.table(decon_output, out.track, sep="\t",
              row.names=TRUE, col.names=NA, quote=FALSE)
  q(status=0)
}

asv_df <- read.csv(file = inp.loc, check.names=FALSE)
rownames(asv_df) <- asv_df[, 1]
asv_df <- asv_df[, -1]
numero_df <- as.matrix(sapply(asv_df, as.numeric))

metadata_df<-read.csv(file = metadata.loc, check.names=FALSE)
rownames(metadata_df) <- metadata_df[, 1]

if(decon.mode == 'prevalence'){
  control_vec <- meta_data_cols(asv_df, metadata_df, prev.control.col)
  true_false_control_vec<-grepl(prev.id.controls,control_vec)
  prev_contam <- isContaminant(numero_df, neg=true_false_control_vec, threshold=threshold, detailed=TRUE, normalize=TRUE, method='prevalence')
  outputer(prev_contam, out.track)
}else if(decon.mode == 'frequency'){
  temp_quant_vec <- meta_data_cols(asv_df, metadata_df, freq.con.col)
  quant_vec<-as.numeric(temp_quant_vec)
  freq_contam <- isContaminant(numero_df, conc=quant_vec, threshold=threshold, detailed=TRUE, normalize=TRUE, method='frequency')
  outputer(freq_contam, out.track)
}else{
  prev_control_vec <- meta_data_cols(asv_df, metadata_df, prev.control.col)
  temp_quant_vec <- meta_data_cols(asv_df, metadata_df, freq.con.col)
  quant_vec<-as.numeric(temp_quant_vec)
  true_false_control_vec<-grepl(prev.id.controls, prev_control_vec)
  comb_contam <- isContaminant(numero_df, neg=true_false_control_vec, conc=quant_vec, threshold=threshold, detailed=TRUE, normalize=TRUE, method='combined')
  outputer(comb_contam, out.track)
}
