#!/usr/bin/R
suppressPackageStartupMessages(library("argparser"))

# create command line options
parser <- arg_parser("Run CopywriteR for a Trio")
parser <- add_argument(parser, "--outdir", help="Where to store the output of the copywriteR calculation")
parser <- add_argument(parser, "--patient.id", help = "Id of the patient")
parser <- add_argument(parser, "--id.file", help="bam file of tumor ID file")
parser <- add_argument(parser, "--rel.file", help="bam file of tumor REL file")
parser <- add_argument(parser, "--cr.file", help="bam file of CR (normal() file")
parser <- add_argument(parser, "--targeted.regions", help="bed file which specifies the capture regions")
parser <- add_argument(parser, "--gene.anno", help="bed file which contains the gene annotations")
parser <- add_argument(parser, "--bin.size", help="The bin size for the read count data")


# read command line options
argv <- parse_args(parser)

argvlist <- unlist(argv)[-1:-2]
#if (any(is.na(argvlist))) {
#  print(argvlist[is.na(argvlist)])
#  stop("Please supply all required arguments!")
#} else {
  print(argvlist)
#}

suppressPackageStartupMessages(library(CopywriteR))
suppressPackageStartupMessages(library(GenomicFeatures))

sample_cpwr_dir = argv$outdir 
input_file_tumour = argv$id.file
input_file_tumour_REL = argv$rel.file
input_file_normal = argv$cr.file
capture_regions_bed = argv$targeted.regions
gene_anno_bed = argv$gene.anno
patient_id = argv$patient.id
bin_size = as.integer(argv$bin.size)


p.id = sprintf("%s_ID", patient_id)
p.rel = sprintf("%s_REL", patient_id)
p.cr = sprintf("%s_CR", patient_id)
id.vs.cr = sprintf("%s.vs.%s", p.id, p.cr)
rel.vs.cr = sprintf("%s.vs.%s", p.rel, p.cr)

bin_kb = substr(as.character(bin_size), 1,2)
reffolder = paste("hg19_", bin_kb, "kb", sep="")

print(paste("Creating reffolder", reffolder))

#sample_cpwr_dir = "/home/mpschr/Documents/projects/exon/out_triotest/PE08_copywriteR"
#input_file_tumour = '/home/mpschr/Documents/projects/exon/out_triotest/PE08_ID/PE08_ID-ready.bam'
#input_file_tumour_REL = '/home/mpschr/Documents/projects/exon/out_triotest/PE08_REL/PE08_REL-ready.bam'
#input_file_normal = '/home/mpschr/Documents/projects/exon/out_triotest/PE08_CR/PE08_CR-ready.bam'
#capture_regions_bed = '/home/mpschr/Data/agilent_sure_select_v5+utr/S04380219_Covered.bed'
#gene_anno_bed = '/home/mpschr/bin/bcbionextgen/data/genomes/Hsapiens/GRCh37/rnaseq-2014-07-14/ref-transcripts.bed'

preCopywriteR(output.folder = file.path(sample_cpwr_dir),
              bin.size = bin_size,
              ref.genome = "hg19",
              prefix = "")

sample.control = data.frame(samples=c(input_file_tumour, input_file_tumour_REL, input_file_normal), 
                            controls=c(input_file_normal, input_file_normal, input_file_normal))

bp.param <- SnowParam(workers = 3, type = "SOCK")


CopywriteR(sample.control = sample.control,
           destination.folder = file.path(sample_cpwr_dir), 
           reference.folder = file.path(sample_cpwr_dir, reffolder),
           capture.regions.file = capture_regions_bed, 
           bp.param = bp.param)
plotCNA(file.path(sample_cpwr_dir))

load(file.path(sample_cpwr_dir,'CNAprofiles','segment.Rdata'))
logdiff_df = as.data.frame(segment.CNA.object$data)

log2_df =read.csv(file.path(sample_cpwr_dir,'CNAprofiles','log2_read_counts.igv'), sep = '\t', comment.char = '#')
print(colnames(log2_df))
colnames(log2_df)[5:7] <- c(p.id, p.rel, p.cr)
print(colnames(log2_df))
log2_df[id.vs.cr] = logdiff_df[,3]
log2_df[rel.vs.cr] = logdiff_df[,4]
log2_df$segment_number = row.names(log2_df)
print(colnames(log2_df))
log2_ranges <-makeGRangesFromDataFrame(log2_df, keep.extra.columns = T)

gene_anno = read.csv(gene_anno_bed, sep = '\t', header = F)[,1:4]
colnames(gene_anno) <- c('Chromosome','Start','End','Gene')
gene_anno_ranges = makeGRangesFromDataFrame(gene_anno, keep.extra.columns = T)

merged_all = mergeByOverlaps(query = gene_anno_ranges, subject = log2_ranges)
wanted_merge_cols = c('log2_ranges.segment_number','log2_ranges.seqnames', 'log2_ranges.start',  'log2_ranges.end', 'Gene')
merged = unique(as.data.frame(merged_all)[wanted_merge_cols])
collapsed = aggregate(Gene ~ ., data = merged, paste, collapse=',')

wanted_out_cols = c('Chromosome', 'Start', 'End', 'Gene', p.id, p.rel, p.cr, id.vs.cr, rel.vs.cr)
out_df = merge.data.frame(log2_df, collapsed[c('log2_ranges.segment_number','Gene')], 
                          by.x = 'segment_number', by.y = 'log2_ranges.segment_number', all.x=TRUE)[wanted_out_cols]
out_df = out_df[order(out_df[,1], out_df[,2]),]


out_file_name = file.path(sample_cpwr_dir,'CNAprofiles', 'log2_CNA.igv')
igv_comment_line = "#track viewLimits=-3:3 graphType=heatmap color=255,0,0"
write( igv_comment_line, file=out_file_name, append=FALSE )
suppressWarnings(write.table(out_df, out_file_name, quote=FALSE, sep = '\t', row.names=FALSE, append = TRUE, na = ''))


# -------------------- segments

segmentation.df = segment.CNA.object$output
segmentation.df$ID = gsub(".ready.bam", "", gsub(c("log2."), "", segmentation.df$ID))
segmentation.df$chrom[segmentation.df$chrom == 23] = 'X'
segmentation.df$chrom[segmentation.df$chrom == 24] = 'Y'
segmentation.df$copynumber = 2^segmentation.df$seg.mean+1
segmentation.df$seg.length.kb = (segmentation.df$loc.end - segmentation.df$loc.start) / 10000
segout_file_name = file.path(sample_cpwr_dir,'CNAprofiles', 'log2_CNA.segmented.tsv')
suppressWarnings(write.table(segmentation.df, segout_file_name, quote=FALSE, sep = '\t', row.names=FALSE,  na = ''))



  
