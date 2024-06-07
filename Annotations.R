## prepare for regular ES update
snp_vcf <- list.files("/storage/lupski/dat/INPUT_hg38/", pattern = "*[0-9|A-Z]_snp.vcf.gz$", recursive = TRUE)
snp_vcf.df <- data.frame(sample = sapply(strsplit(snp_vcf,"/"), '[[', 1), 
                         snp_vcf = sapply(strsplit(snp_vcf,"/"), '[[', 2))

indel_vcf <- list.files("/storage/lupski/dat/INPUT_hg38/", pattern = "*[0-9|A-Z]_indel.vcf.gz$", recursive = T)
indel_vcf.df <- data.frame(sample = sapply(strsplit(indel_vcf,"/"), '[[', 1), 
                           indel_vcf = sapply(strsplit(indel_vcf, "/"), '[[', 2))

sample_vcf <- merge(snp_vcf.df, indel_vcf.df, by = "sample", all = T)%>%
  distinct(snp_vcf, indel_vcf, sample, .keep_all = T) %>%
  filter(!sample %in% c("BAB9637","BAB9638","BAB9639","BH12473-1"), !snp_vcf %like% "BH16500", !snp_vcf %like% "tmp") %>%
  mutate(root_path = "/storage/lupski/dat/INPUT_hg38/")

# for WES data other than HGSC, like from Genedx, BG, etc. (Haowei call these files as "reprocessed files"). On May 17 2024, Haowei shared excel sheet for these reprocessed samples
snp_vf_wgl <- list.files("/storage/lupski/dat/WGL2CMG_hg38/", pattern = "*[0-9|A-Z]_snp.vcf.gz$", recursive = TRUE)
indel_vf_wgl <- list.files("/storage/lupski/dat/WGL2CMG_hg38/", pattern = "*[0-9|A-Z]_indel.vcf.gz$", recursive = TRUE)
snp_vcf_wgl.df <- data.frame(sample = sapply(strsplit(snp_vf_wgl, "/"), '[[',1), 
                             snp_vcf = sapply(strsplit(snp_vf_wgl, "/"), '[[',2))
indel_vcf_wgl.df <- data.frame(sample = sapply(strsplit(indel_vf_wgl,"/"), '[[',1), 
                               indel_vcf = sapply(strsplit(indel_vf_wgl,"/"), '[[',2))

sample_vcf_wgl <- merge(snp_vcf_wgl.df ,indel_vcf_wgl.df ,by = "sample",all = T) %>%
  distinct(snp_vcf, indel_vcf, sample, .keep_all = T) %>%
  filter(!sample %in% c("BAB9637","BAB9638","BAB9639","BH12473-1"), !snp_vcf%like%"BH16500", !snp_vcf%like%"tmp") %>%
  mutate(root_path = "/storage/lupski/dat/WGL2CMG_hg38/")
sample_vcf_wes <- rbindlist(list(sample_vcf, sample_vcf_wgl))
#### "sample_vcf_wes" is the input file for annotation


# # Function to do the entire VA update. Need to update. It's incomplete
# RUN_update_hg38 <- function(meta){
#   variant_key_path <- "~/VA/hg38_variant_key/"
#   # somalier is the program for relatedness calculation
#   somalier_extracted_path <- "~/VA/hg38_somalier_extracted/"
#   somalier_path <- "/storage/lupski/app/bin/somalier" #v0.2.17
#   merge_query <- function(sample,out_path,snp_path,indel_path,bcftool_path="/usr/local/bin/bcftools"){
#     cmd <- glue::glue('{bcftool_path} concat -a {snp_path} {indel_path} | {bcftool_path} query -i ',
#                       "'FILTER=", '"PASS"', "' -f", " '[%CHROM-%POS-%REF-%ALT]\\t%FORMAT[\\t%GT\\t%VR\\t%RR\\t%SAMPLE]\\n'",
#                       ' > {out_path}',
#                       bcftool_path=bcftool_path,
#                       snp_path=snp_path,
#                       indel_path=indel_path,
#                       out_path=out_path)
#     system(cmd)
#   }
  
  # ## for WGS
  # sample_vcf_wgs <- meta %>% filter(APPLICATION=="WHOLE GENOME SEQUENCING",!is.na(snp_vcf),!is.na(indel_vcf))
  # lapply(1:nrow(sample_vcf_wgs), function(i){
  #   print(paste0("WGS:",i))
  #   raw_path <- "/storage/lupski/dat/WGS_hg38/"
  #   process_path <- "/storage/lupski/home/va/data/raw_data_hg38/"
  #   raw_dir <- paste0(raw_path,sample_vcf_wgs$sample[i])
  #   process_dir <- paste0(process_path,sample_vcf_wgs$sample[i])
  #   dir.create(process_dir)
  #   sample=as.character(sample_vcf_wgs$sample[i])
  #   snp_path=paste0(raw_path,sample,"/",as.character(sample_vcf_wgs$snp_vcf[i]))
  #   indel_path=paste0(raw_path,sample,"/",as.character(sample_vcf_wgs$indel_vcf[i]))
  #   out_path=paste0(process_dir,"/xatlas.pass.variant.tsv")
  #   key_path=paste0(variant_key_path,sample_vcf_wgs$sample[i],"_wgs_variant_key.tsv")
  #   if(file.exists(gsub("tsv","parquet",out_path))){
  #     if(!file.exists(key_path)){
  #       arrow::read_parquet(gsub("tsv","parquet",out_path),as_data_frame = F)%>%
  #         filter(V5+V6>20,grepl("chr[0-9]+-[1-9]|X",V1))%>%
  #         mutate(V5=as.numeric(V5),V6=as.numeric(V6),V8=round(V5/(V5+V6),digits = 3))%>%
  #         mutate(zyg=ifelse(V8>0.1&V8<0.9,1,ifelse(V8>0.9,2,0)))%>%
  #         dplyr::select(V1,zyg)%>%
  #         dplyr::collect()%>%write.table(.,file = key_path,quote = F,row.names = F,col.names = F,sep = " ")
  #     }
  #   } else {
  #     print(paste("## Merge vcf for",sample))
  #     system(paste("tabix -p vcf",snp_path))
  #     system(paste("tabix -p vcf",indel_path))
  #     merge_query(sample,out_path,snp_path,indel_path,bcftool_path="/usr/local/bin/bcftools")
  #     data.table::fread(out_path,stringsAsFactors = F)%>%arrow::write_parquet(.,gsub("tsv","parquet",out_path))
  #     arrow::read_parquet(gsub("tsv","parquet",out_path),as_data_frame = F)%>%
  #       filter(V5+V6>20,grepl("chr[0-9]+-[1-9]|X",V1))%>%
  #       mutate(V5=as.numeric(V5),V6=as.numeric(V6),V8=round(V5/(V5+V6),digits = 3))%>%
  #       mutate(zyg=ifelse(V8>0.1&V8<0.9,1,ifelse(V8>0.9,2,0)))%>%
  #       dplyr::select(V1,zyg)%>%
  #       dplyr::collect()%>%write.table(.,file = key_path,quote = F,row.names = F,col.names = F)
  #     }
  #   })
  # wgs_key_path <- list.files("hg38_variant_key",pattern = "*wgs_variant_key.tsv",full.names = T)
  # cmd <- paste0("./countKeys ", paste(wgs_key_path,collapse = " ")," > total_wgs_key.tsv")
  # system(cmd,intern = F)
  #sample_vcf_wes <- meta %>% filter(APPLICATION=="WHOLE EXOME SEQUENCING",!is.na(snp_vcf),!is.na(indel_vcf))
  
  ## for WES
  
  merge_query <- function(sample, out_path, snp_path, indel_path, bcftool_path = "/usr/local/bin/bcftools"){
    cmd <- glue::glue('{bcftool_path} concat -a {snp_path} {indel_path} | {bcftool_path} query -i ',
                      "'FILTER=", '"PASS"', "' -f", " '[%CHROM-%POS-%REF-%ALT]\\t%FORMAT[\\t%GT\\t%VR\\t%RR\\t%SAMPLE]\\n'",
                      ' > {out_path}',
                      bcftool_path = bcftool_path,
                      snp_path = snp_path,
                      indel_path = indel_path,
                      out_path = out_path)
    system(cmd)
  }
  
  # from vcf files, extract variant information
  lapply(1:nrow(sample_vcf_wes), function(i){
    print(paste0("WES:", i))
    # define variables 
    sample = as.character(sample_vcf_wes$sample[i])
    
    ## PATHS
    # define existing input and output paths and create new directories based on sample name
    # input
    raw_path <- sample_vcf_wes$root_path[i]
    raw_dir <- paste0(raw_path, sample_vcf_wes$sample[i])
    snp_path = paste0(raw_path, sample, "/", as.character(sample_vcf_wes$snp_vcf[i]))
    indel_path = paste0(raw_path,sample, "/", as.character(sample_vcf_wes$indel_vcf[i]))
    
    # main output
    process_path <- "/storage/lupski/home/va/data/raw_data_hg38/"
    process_dir <- paste0(process_path, sample_vcf_wes$sample[i])
    dir.create(process_dir)
    out_path = paste0(process_dir, "/xatlas.pass.variant.tsv")
    
    # other outputs
    # somalier is the program for relatedness calculation
    somalier_path <- "/storage/lupski/app/bin/somalier" #v0.2.17
    somalier_extracted_path <- "~/VA/hg38_somalier_extracted/"
    sample_extract_path = paste0(somalier_extracted_path, sample, ".somalier")
    # variant key
    variant_key_path <- "~/VA/hg38_variant_key/"
    key_path = paste0(variant_key_path, sample_vcf_wes$sample[i], "_wes_variant_key.tsv")
    index_key_path = paste0(variant_key_path, sample_vcf_wes$sample[i], "_wes_key_indexed.tsv")
    
    # parquet is used to handle files across different platforms like R, SQL, python, etc.
    # check existence of files
    # gsub because tsv file is stored in the out_path 
    parquet_exist = file.exists(gsub("tsv", "parquet", out_path))
    #parquet_exist=FALSE
    key_path_exist = file.exists(key_path)
    index_key_path_exist = file.exists(index_key_path)
    #index_key_path_exist=FALSE
    if(parquet_exist){
      if(!key_path_exist){
        arrow::read_parquet(gsub("tsv", "parquet", out_path), as_data_frame = F) %>%
          # I think creating a key for each variant (basically numbering the variant as they occur in order per sample) 
          # And selecting only those variants with total  ref+alt read > 20 (I checked for sample BH15932-1 and found that this filtering step removed almost half of the variants!)
          filter(V5 + V6 > 20, grepl("chr[0-9]+-[1-9]|X", V1)) %>%
          # calculate VR/TR ratio in V8 column
          mutate(V5 = as.numeric(V5), V6 = as.numeric(V6), V8 = round(V5/(V5+V6), digits = 3)) %>%
          # Heterozygotes VR/TR ratio is between 0.1-0.9
          mutate(zyg = ifelse(V8 > 0.1 & V8 < 0.9, 1, ifelse(V8>0.9, 2, 0))) %>%
          dplyr::select(V1, zyg) %>%
          dplyr::collect() %>% 
          # write key in /storage/lupski/home/va/VA/hg38_variant_key
          write.table(.,file = key_path,quote = F,row.names = F,col.names = F)
      }
    } else if(!parquet_exist){
      print(paste("## Merge vcf for", sample))
      system(paste("tabix -p vcf", snp_path))
      system(paste("tabix -p vcf", indel_path))
      merge_query(sample, out_path, snp_path, indel_path, bcftool_path = "/usr/local/bin/bcftools")
      data.table::fread(out_path, stringsAsFactors = F) %>% arrow::write_parquet(., gsub("tsv", "parquet", out_path))
      arrow::read_parquet(gsub("tsv", "parquet", out_path), as_data_frame = F)%>%
        filter(V5 + V6 > 20, grepl("chr[0-9]+-[1-9]|X",V1))%>%
        mutate(V5 = as.numeric(V5), V6 = as.numeric(V6), V8 = round(V5/(V5+V6), digits = 3)) %>%
        mutate(zyg = ifelse(V8 > 0.1 & V8 < 0.9, 1, ifelse(V8 > 0.9, 2, 0)))%>%
        dplyr::select(V1, zyg) %>%
        dplyr::collect() %>%
        write.table(., file = key_path,quote = F,row.names = F,col.names = F)
    } else {
      print(paste("write out key for:", sample))
      arrow::read_parquet(gsub("tsv", "parquet", out_path), as_data_frame = F) %>%
        filter(V5 + V6 > 20, grepl("chr[0-9]+-[1-9]|X", V1)) %>%
        dplyr::select(V1) %>%
        dplyr::collect() %>% 
        write.table(., file = key_path)
    }
    if(!index_key_path_exist){
      print(paste("index variant key for", sample))
      # indexing keys
      # rsvr enc: Encode input variants specified as CHROM, POS, REF ALT as RSVR IDs.
      cmd <- paste0("awk -F'[-:]' '{print $1, $2, $3, $4}' ", 
                    key_path, 
                    " | sed 's/^chr//' | /storage/lupski/app/rsvr-master/rsvr enc -s -p -d ' ' > ", 
                    "hg38_variant_key/", sample, "_wes_key_indexed.tsv")
      system(cmd)
    }
    
    
    ## deprecated
    # if(!file.exists(sample_extract_path)){
    #   cmd=paste0(somalier_path," extract -d hg38_somalier_extracted/ --sites sites.hg38.vcf.gz -f GRCh38_1000Genomes/GRCh38_full_analysis_set_plus_decoy_hla.fa ",snp_path)
    #   system(cmd,intern = F)
    # }
  })
  wes_key_path <- list.files("hg38_variant_key", pattern = "*_wes_variant_key.tsv", full.names = T)
  write.table(file = "wes_key_path.list", wes_key_path, quote = F, row.names = F, col.names = F)
  cmd2 <- paste0("./countKeys_v2 -f wes_key_path.list > total_wes_key.tsv")
  system(cmd2, intern = F)
  # get keys with more than 20 alleles
  cmd3 <- "awk -F':' '$2 > 20 { print $1 }' total_wes_key.tsv > hg38_internal_prox0.001_key"
  system(cmd3, intern = F)
  # remove common variant from gnomad and esp6500
  cmd4 <- "./subtractKeys_v2 total_wes_key.tsv hg38_esp6500siv2_0.001_key hg38_gnomad211_exome_0.001_key hg38_internal_prox0.001_key > total_wes_key_subtracted.tsv"
  system(cmd4, intern = F)
  # index the rare variant database
  #total_wes_key_subtracted <- fread("total_wes_key_subtracted.tsv", sep = ":", header = F)
  # Filter rows where variant key has "chr" in it. The ones without "chr" are from hg19 and the ones with "chr" are from hg38
  total_wes_key_subtracted <- fread(cmd = "awk -F'[-:]' '{print $1, $2, $3, $4, $5}' total_wes_key_subtracted.tsv", header = F)[V1 %in% paste0("chr",c(1:22,"X")),]
  total_wes_key_subtracted[, tempcol := gsub("chr", "", V1)]
  
  # Next five lines rearrange the data.table by the increasing order of chr and position
  total_wes_key_subtracted[, col1_transformed := as.numeric(replace(tempcol, tempcol == "X", 23))]
  setorder(total_wes_key_subtracted, col1_transformed, V2)
  total_wes_key_subtracted[, col1_transformed := NULL]
  total_wes_key_subtracted[, tempcol := NULL]
  
  write.table(total_wes_key_subtracted,file = "total_wes_key_subtracted.sorted.tsv", quote = F, sep = "\t",row.names = F,col.names = F)
  cmd5 <-"awk '{print $1, $2, $3, $4}' total_wes_key_subtracted.sorted.tsv|sed 's/^chr//'|/storage/lupski/app/rsvr-master/rsvr enc -s -p -d ' '> total_wes_key_subtracted_indexed.tsv"
  system(cmd5, intern = F)
  
  ## append variant_key table to the SQL database
  # all keys
  include_keys <- fread("total_wes_key_subtracted_indexed.tsv", header = F, select = "V1")
  # Function to append
  # Later this function will be called on wes_key_path which has files for indexed keys for each sample
  index_files <- function(filenames) {
    conn <- dbConnect(SQLite(), dbname = "va.sqlite")
    # Get the list of samples already in the database
    existing_samples <- dbGetQuery(conn, "SELECT DISTINCT SAMPLE FROM hg38_keys")
    # Load keys to be included
    for (filename in filenames) {
      sample <- sub("\\_wes_key_indexed.tsv", "", basename(filename))  # Strip text from the filename
      
      # Check if sample in sqldb has already been indexed
      if(sample %in% existing_samples$sample) {
        message(paste("Sample", sample, "has already been indexed. Skipping."))
        next
      }
      
      # Read file into a data frame and store as hg38_keys sql table
      data <- fread(filename, header = F)[, c(1, 6)] %>% distinct(V1, .keep_all = T)
      setnames(data, c("V1", "V6"), c("RSVR_ID", "GT"))
      data <- filter(data, RSVR_ID %in% include_keys$V1)
      data$SAMPLE <- sample
      dbWriteTable(conn, "hg38_keys", data, append = TRUE)
    }
    dbDisconnect(conn)
  }
  
  # call the above function
  wes_key_path <- list.files("hg38_variant_key", pattern = "*_wes_key_indexed.tsv", full.names = T)
  ## create rvid table
  tryCatch({
    index_files(wes_key_path)
  }, error = logError)
  
  ## perform variant count and store in hg38_key_counts sql table
  count_key_files <- function() {
    conn <- dbConnect(SQLite(), dbname = "va.sqlite")
    
    # Drop the table if it already exists
    dbExecute(conn, "DROP TABLE IF EXISTS hg38_key_counts")
    
    # Perform the query and save results in a new table
    dbExecute(conn, "
    CREATE TABLE hg38_key_counts AS
    SELECT RSVR_ID, COUNT(DISTINCT SAMPLE) as sample_counts
    FROM hg38_keys
    WHERE GT > 0
    GROUP BY RSVR_ID
  ")
    dbDisconnect(conn)
  }
  tryCatch({
    count_key_files()
  }, error = logError)
  
  ## perform sample concat
  conn <- dbConnect(SQLite(), dbname = "va.sqlite")
  cmg <- tbl(conn, "hg38_keys") %>% 
    collect()
  nsample <- length(unique(cmg$SAMPLE))*2
  concat_tbl <- cmg %>%
    group_by(RSVR_ID) %>%
    # n_distinct is equivalent of length(unique(x))
    summarise(samples = ifelse(n_distinct(SAMPLE) <= 10, toString(unique(SAMPLE)), '>10'),
              cmg_freq = sum(GT)/nsample)
  #dbWriteTable(conn, "hg38_keys_concat", data.frame(RSVR_ID = integer64(), samples = character(),cmg_freq=double()), overwrite = TRUE)
  dbWriteTable(conn, "hg38_keys_concat", concat_tbl, overwrite = TRUE)
  
  ## write VARIANT table
  indexed_keys <- fread("total_wes_key_subtracted_indexed.tsv", stringsAsFactors = F, header = F)
  names(indexed_keys) <- c("RSVR_ID","CHROM","POS","REF","ALT")
  
  indexed_keys <- indexed_keys[CHROM %in% c(1:22, "X", "Y"),] %>% distinct(RSVR_ID, .keep_all = T) %>% arrange(RSVR_ID)
  dbWriteTable(conn, "VARIANT", indexed_keys, row.names = FALSE, overwrite = TRUE)
  unanno_keys <- as.data.table(indexed_keys)
  unanno_keys[, strand := "+"]
  unanno_keys[, chr := paste0("chr", CHROM)]
  write.table(unanno_keys[,c("chr", "POS", "strand", "REF", "ALT", "RSVR_ID")], file = "unanno_key.tsv", quote = F, sep = "\t", row.names = F, col.names = T)
  
  ## using opencravat for the annotation
  # First, scp the unanno_key.tsv into your local directory
  # scp va@10.66.4.211:/storage/lupski/home/va/VA/unanno_key.tsv .
  # Second, go to OpenCravat -> upload unanno_key.tsv file -> select following options:
  # AlphaMissense, CADD Exome, ClinGen Allele Ancestry, ClinVar, DIDA, GnomAD3, LoFtool, OMIM, PhyloP, PolyPhen-2,
  # SIFT, REVEL, SwissProt-Domains, SpliceAI

  alphamissense <- c("RSVR_ID","alphamissense.am_pathogenicity","alphamissense.am_class","alphamissense.transcript_id","alphamissense.uniprot_id")
  cadd<- c("RSVR_ID", "cadd.score", "cadd.phred")
  clingen <- c("RSVR_ID","clingen_allele_registry.allele_registry_id","clingen_allele_registry.disease","clingen_allele_registry.mode_of_inheritance",
               "clingen_allele_registry.assertion","clingen_allele_registry.evidence_codes","clingen_allele_registry.summary_of_interpretation")
  clinvar <- c("RSVR_ID","clinvar.sig","clinvar.disease_refs","clinvar.disease_names","clinvar.rev_stat","clinvar.id","clinvar.sig_conf","clinvar_acmg.ps1_id","clinvar_acmg.pm5_id")
  dida <- c("RSVR_ID","dida.name","dida.effect","dida.relation","dida.fam","dida.funct","dida.dist","dida.pub","dida.all")
  loftool <- c("RSVR_ID","loftool.loftool_score")
  omim <- c("RSVR_ID","omim.omim_id")
  mapping <- c("RSVR_ID", "coding", "hugo", "transcript", "so", "exonno", "cchange", "achange", "all_mappings")
  phylop <- c("RSVR_ID","phylop.phylop100_vert","phylop.phylop100_vert_r","phylop.phylop30_mamm","phylop.phylop30_mamm_r","phylop.phylop17_primate","phylop.phylop17_primate_r")
  polyphen <- c("RSVR_ID","polyphen2.hdiv_pred","polyphen2.hvar_pred","polyphen2.hdiv_rank","polyphen2.hvar_rank")
  revel <- c("RSVR_ID","revel.transcript","revel.score","revel.rankscore","revel.all")
  sift <- c("RSVR_ID","sift.transcript","sift.prediction","sift.confidence","sift.score","sift.rankscore","sift.med","sift.seqs","sift.all","sift.multsite")
  spliceai <- c("RSVR_ID","spliceai.ds_al","spliceai.ds_ag","spliceai.ds_dl","spliceai.ds_dg","spliceai.dp_al","spliceai.dp_ag","spliceai.dp_dl","spliceai.dp_dg")
  swissprot <- c("RSVR_ID","swissprot_domains.uniprotkb","swissprot_domains.domain","swissprot_domains.intramem","swissprot_domains.motif","swissprot_domains.peptide","swissprot_domains.repeat","swissprot_domains.topo","swissprot_domains.transmem","swissprot_domains.pubmed","swissprot_domains.all")
  gnomad3 <- c("RSVR_ID","gnomad3.af","gnomad3.af_afr","gnomad3.af_asj","gnomad3.af_eas","gnomad3.af_fin","gnomad3.af_lat","gnomad3.af_nfe","gnomad3.af_oth","gnomad3.af_sas")
  
  
  # read the output from OpenCravat 
  anno_key <- fread("grep -v '#' 20240607_unanno_key.tsv.variant.tsv", header = T)
  setnames(anno_key, "samples", "RSVR_ID")
  
  ## seperate the tables
  # anno_key[complete.cases(anno_key[,..cadd]),..cadd]
  # anno_key[complete.cases(anno_key[,..mapping]),..mapping]
  anno_key_filtered <- anno_key[gnomad3.af<0.01 | is.na(gnomad3.af), ]
  
  ## additional annotation
  pidd_2021 <- readLines("/storage/lupski/home/va/workspace/variant_analyzer/r_scripts/gene_lists/PIDD_2021.txt")
  pidd_2022 <- readLines("/storage/lupski/home/va/workspace/variant_analyzer/r_scripts/gene_lists/PIDD_2022.txt")
  piddPredicted_2021 <- readLines("/storage/lupski/home/va/workspace/variant_analyzer/r_scripts/gene_lists/PIDD_predicted_2021.txt")
  piddPredicted_2022 <- readLines("/storage/lupski/home/va/workspace/variant_analyzer/r_scripts/gene_lists/PIDD_predicted_2022.txt")
  pidd_2022_extend <- readLines("/storage/lupski/home/va/workspace/variant_analyzer/r_scripts/gene_lists/PIDD_2022_extended.txt")
  cmt_genes <- readLines("/storage/lupski/home/va/workspace/variant_analyzer/r_scripts/gene_lists/cmt_claudia_genes.txt")
  tbm_candidate_genes <- readLines("/storage/lupski/home/va/workspace/variant_analyzer/r_scripts/gene_lists/tbm_candidate_genes.txt")
  tbm_extended_genes <- readLines("/storage/lupski/home/va/workspace/variant_analyzer/r_scripts/gene_lists/tbm_extended_genes.txt")
  pLI = read.table('db/fordist_cleaned_exac_r03_march16_z_pli_rec_null_data.txt', header = T) %>%
    dplyr::select(gene, pLI, pRec, pNull)
  mim2gene <- read.csv('~/VA/db/mim2gene.txt', sep='\t', stringsAsFactors = FALSE, header = TRUE, skip = 4)
  genemap2 <- read.csv('~/VA/db/genemap2.txt', sep='\t', stringsAsFactors = FALSE, header = TRUE, skip = 3)
  OMIM_gene_inh <- mim2gene %>% 
    mutate(MIM.Number = X..MIM.Number, MIM.type = MIM.Entry.Type..see.FAQ.1.3.at.https...omim.org.help.faq.) %>% 
    merge(., genemap2, by='MIM.Number') %>% 
    select(MIM.Number, MIM.type, Approved.Symbol, Phenotypes, Ensembl.Gene.ID) %>%
    mutate(AD = ifelse(Phenotypes %like% "Autosomal dominant | autosomal dominant", "AD", ""),
           AR = ifelse(Phenotypes %like% "Autosomal recessive | autosomal recessive", "AR", ""))%>%
    mutate(Raw_gene_name = Approved.Symbol,
           Phenotypes_OMIM = paste0(MIM.Number, ";", Phenotypes)) %>%
    dplyr::select(Raw_gene_name, Phenotypes_OMIM, AD, AR)
  anno_key_filtered$pidd_2021 <- anno_key_filtered$hugo %in% pidd_2021
  anno_key_filtered$pidd_2022 <- anno_key_filtered$hugo %in% pidd_2022
  anno_key_filtered$piddPredicted_2021 <- anno_key_filtered$hugo%in% piddPredicted_2021
  anno_key_filtered$piddPredicted_2022 <- anno_key_filtered$hugo%in% piddPredicted_2022
  anno_key_filtered$pidd_2022_extend <- anno_key_filtered$hugo%in% pidd_2022_extend
  
  anno_key_filtered$CMT <- anno_key_filtered$hugo %in% cmt_genes
  anno_key_filtered$TBM <- anno_key_filtered$hugo %in% tbm_candidate_genes
  anno_key_filtered$TBM_extended <- anno_key_filtered$hugo %in% tbm_extended_genes
  
  ## annotate pedgree
  ## 'cmg_pedgree_tbl' stores the pedgree information 
  cmg_pedgree <- tbl(conn, "cmg_pedgree_tbl") %>% 
    collect()
  existing_samples <- dbGetQuery(conn, "SELECT DISTINCT SAMPLE FROM hg38_keys")
  cmg <- tbl(conn, "hg38_keys") %>% 
    collect()
  
  miss_pedlist <- existing_samples$SAMPLE[!existing_samples$SAMPLE %in% unique(unlist(cmg_pedgree))]
  tmplist <- lapply(miss_pedlist, function(sample){
    print(sample)
    if(!sample %in% cmg_pedgree$Parent1 &!sample %in% cmg_pedgree$Parent2 &!sample %in% cmg_pedgree$Proband){
      uniq_var <- fread(paste0("hg38_variant_key/",sample,"_wes_key_indexed.tsv"), header = F)[,c(1,6)] %>% 
        distinct(V1, .keep_all = T)
      setnames(uniq_var, c("V1","V6"), c("RSVR_ID","GT"))
      uniq_var <- filter(uniq_var, 
                         RSVR_ID %in% include_keys$V1)
      # based on hg38_keys
      cmg_samples <-  cmg %>%
        filter(RSVR_ID %in% local(uniq_var$RSVR_ID)) %>%
        select(SAMPLE) %>%
        collect()
      rltbl <- sort(table(unlist(cmg_samples$SAMPLE)), decreasing = T)
      if (length(rltbl) < 5) {
        return(NULL)
      }
      parent1 <- parent2 <- ""
      if(sample %like% "_XATLAS"){
        rltbl <- sort(c(rltbl[!duplicated((gsub("_XATLAS", "", names(rltbl))), fromLast = TRUE) & 
                                !duplicated((gsub("_XATLAS", "", names(rltbl))))], 
                        rltbl[duplicated((gsub("_XATLAS", "", names(rltbl))))]), decreasing = T)
        pt <- sample
        if (length(rltbl) < 3) {
          return(NULL)
        }
        if (rltbl[2] > 50 & rltbl[2] / rltbl[1] > 0.2) parent1 <- names(rltbl)[2]
        if (rltbl[3] > 50 & rltbl[3] / rltbl[1] > 0.2) parent2 <- names(rltbl)[3]
      }else{
        pt <- names(rltbl)[1]
        if (rltbl[2] > 50 & rltbl[2] / rltbl[1] > 0.2) parent1 <- names(rltbl)[2]
        if (rltbl[3] > 50 & rltbl[3] / rltbl[1] > 0.2) parent2 <- names(rltbl)[3]
      }
      return(data.table(Proband=pt, Parent1=parent1, Parent2=parent2))
    }
  })
  ped_add <- rbindlist(tmplist)
  ped_add$Affected <- ""
  ped_add$Gender <- ""
  cmg_pedgree2 <- rbindlist(list(cmg_pedgree,ped_add)) %>%
    distinct(Proband,.keep_all = T) %>%
    filter(!is.na(Proband))
  setDT(cmg_pedgree2)
  # note: some cases are in hg37(BH15492,BH15495,BH15496,BH15505,BH15506,BH15572,BH15578,BH15579)
  dbWriteTable(conn, "cmg_pedgree_tbl", cmg_pedgree2, overwrite=T, row.names=F)
  cmg_pedgree_trios <- tbl(conn, "cmg_pedgree_tbl") %>%
    filter(Parent1 != "", Parent2 != "") %>%
    collect() %>%
    filter(!Proband %like% "-2|-3")%>%
    filter(!Parent1%like%"-1_XATLAS" & !Parent2 %like% "-1_XATLAS") %>%
    filter(gsub("_XATLAS", "", Parent1) != Parent2 & gsub("_XATLAS","", Parent2) != Parent1)
  setDT(cmg_pedgree_trios)
  IDs <- unique(cmg$SAMPLE)
  
  lapply(1:length(IDs), function(i){
    print (paste0("write out reannotate samples: ", i))
    sample = IDs[i]
    process_path <- "/storage/lupski/home/va/data/raw_data_hg38/"
    process_dir <- paste0(process_path, sample)
    parquet_path = paste0(process_dir, "/xatlas.pass.variant.parquet")
    xls31_anno = paste0(process_dir, "/xatlas.pass.filtered.variant.csv")
    #ll_xls_v2 = paste0(process_dir, "/xatlas.pass.filtered.anno.variant.csv")
    ll_xls_v3 = paste0(process_dir, "/xatlas.pass.filtered_0.01.anno.variant.csv")
    VEP_input = paste0(process_dir, "/xatlas.pass.filtered.anno.variant.vcf")
    RVSD_anno = paste0(process_dir, "/RVSD_anno.csv")
    if ( !file.exists(parquet_path)){
      print ("parquet file does not exists")
      print (sample)		
    }
    if(!file.exists(RVSD_anno)){
      variants_indexed <- fread(paste0("hg38_variant_key/", sample, "_wes_key_indexed.tsv"))
      rare_rsvr <- indexed_keys$RSVR_ID[indexed_keys$RSVR_ID %in% variants_indexed$V1]
      variants_indexed <- variants_indexed[V1 %in% rare_rsvr]
      setnames(variants_indexed,"V1","RSVR_ID")
      keys <-  variants_indexed %>%
        mutate(key = paste0("chr", V2, "-", V3, "-", V4, "-", V5)) %>%
        select(key, RSVR_ID)
      variants <- arrow::read_parquet(parquet_path,as_data_frame = F) %>%
        filter(grepl("chr[0-9]+-[1-9]|X", V1)) %>%
        collect() %>%
        merge(., keys,by.x = "V1", by.y="key", all.y=T) %>%
        as.data.table()
      variants$project_name <- "N/A"
      variants$InhFrom <- "N/A"
      variants$Compound_het <- "N/A"
      variants$Potential_Parent1 <- "N/A"
      variants$Potential_Parent2 <- "N/A"
      variants_p1 <- cmg[which(cmg$SAMPLE == p1), ]
      variants_p2 <- cmg[which(cmg$SAMPLE == p2), ]
      if(sample %in% cmg_pedgree_trios$Proband & nrow(variants_p1) != 0 & nrow(variants_p2) != 0){
        p1 <- cmg_pedgree_trios[Proband == sample, ]$Parent1
        p2 <- cmg_pedgree_trios[Proband == sample, ]$Parent2
        variants <- variants %>%
          mutate(V5 = as.numeric(V5), V6 = as.numeric(V6)) %>%
          #filter(V5+V6 > 20, V5 > 2) %>%
          mutate(InhFrom = ifelse(!RSVR_ID %in% variants_p1$RSVR_ID & !RSVR_ID %in% variants_p2$RSVR_ID, "DeNovo",
                                ifelse(RSVR_ID %in% variants_p1$RSVR_ID & RSVR_ID %in% variants_p2$RSVR_ID, paste0(p1, ";", p2),
                                       ifelse(RSVR_ID %in% variants_p1$RSVR_ID & !RSVR_ID %in% variants_p2$RSVR_ID, p1, p2))))%>%
          mutate(Zygosity_V2 = ifelse(V5/(V6+V5) > 0.35 & V5/(V6+V5) <0.65, "Het",
                                    ifelse(V5/(V6+V5) > 0.9,"Hom", "Het-mosaicism"))) %>%
          merge(., anno_key[, c("RSVR_ID", "hugo")], by= "RSVR_ID") %>%
          group_by(hugo) %>%
          summarise(
            Identifier = V7,
            key = V1,
            RSVR_ID = RSVR_ID,
            Zygosity = V4,
            Zygosity_V2 = Zygosity_V2,
            vR = V5,
            tR = V5+V6,
            Gene.refGene = hugo,
            InhFrom = InhFrom,
            Compound_het = ifelse(n() > 1, "TRUE", "FALSE")) %>% 
          ungroup()
        variants$Potential_Parent1 <- p1
        variants$Potential_Parent2 <- p2
        
      } 
      else {
        variants <- variants%>%
          mutate(V5 = as.numeric(V5),V6 = as.numeric(V6)) %>%
          filter(V5+V6>20, V5>2) %>%
          mutate(Zygosity_V2 = ifelse(V5/(V6+V5) > 0.35 & V5/(V6+V5) < 0.65,"Het",
                                    ifelse(V5/(V6+V5) > 0.9,"Hom", "Het-mosaicism"))) %>%
          merge(., anno_key[, c("RSVR_ID","hugo")], by = "RSVR_ID") %>%
          group_by(hugo) %>%
          summarise(
            Identifier = V7,
            key = V1,
            RSVR_ID = RSVR_ID,
            Zygosity = V4,
            Zygosity_V2 = Zygosity_V2,
            vR = V5,
            tR = V5+V6,
            Gene.refGene = hugo,
            InhFrom = InhFrom,
            Compound_het = ifelse(n()>1, "TRUE", "FALSE")) %>%
          ungroup()
        variants$Potential_Parent1 <- "NA"
        variants$Potential_Parent2 <- "NA"
      }
      info <- rep('.', nrow(variants))
      identifier <- rep(".", nrow(variants))
      QUAL <- rep(".", nrow(variants))
      SAMPLE <- rep(".", nrow(variants))
      add_variant <- variants$key
      varPosList <- strsplit(add_variant, ":|_|>|-")
      varPosList <- varPosList[sapply(varPosList,length) == 4]
      merged <- cbind(sapply(varPosList,"[[",1),
                      sapply(varPosList,"[[",2),
                      identifier, 
                      sapply(varPosList,"[[",3),
                      sapply(varPosList,"[[",4), QUAL, info, SAMPLE)
      chrOrder <- c(paste("chr", 1:22, sep = ""), "chrX", "chrY", "chrM")
      merged <- data.table::as.data.table(merged)
      merged$V1 <- factor(merged$V1, levels = chrOrder)
      merged <- merged %>%
        mutate(V2 = as.numeric(V2)) %>%
        arrange(V1, V2)
      write.table(merged,file=VEP_input,quote=FALSE,sep='\t',row.names = FALSE,col.names=FALSE)
      system(paste0("bgzip -f ", VEP_input))
      write.csv(variants,file=xls31_anno,row.names = F)
    }
    if(TRUE){
      print(paste("write out reannotate file for ", sample))
      variants <- fread(xls31_anno,stringsAsFactors = F)%>%
        merge(., concat_tbl, by = "RSVR_ID", all.x=T)%>%
        merge(., anno_key_filtered, by = "RSVR_ID", all.x=T)%>%
        merge(., pLI, by.x="Gene.refGene", by.y = "gene", all.x=T)%>%
        select(-c("hugo.x", "hugo.y", "uid", "numsample", "note_variant","tags"))%>%
        distinct(RSVR_ID, .keep_all = T)
      # write.csv(variants,file = ll_xls_v2,row.names = F,na = ".")
      write.csv(variants, file = ll_xls_v3, row.names = F,na = ".")
    }
    
  })
  rm(list = ls())
  gc()
  
  # EDN of the RSID based annotaion 
  ##-----   
