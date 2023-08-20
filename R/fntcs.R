#' load genes of somtic TOP-ART criterion as character vector#'
#' @return A character vector.
#' @export
load_somatic_topart_genes <- function(){
  c('C1orf86','MTOR','MAD2L2','RPA2','SFPQ','PLK3','RAD54L','USP1','GADD45A',
    'CDC7','PRMT6','PARP1','SMC6','GEN1','DNMT3A','FANCL','ERCC3','MCM6',
    'PSMD14','NABP1','SUMO1','INO80D','BARD1','SMARCAL1','STK36','FANCD2',
    'TOP2B', 'WDR48','BAP1','POLQ','MCM2','TOPBP1','ATR','RNF168','POLN','HELQ',
    'FAM175A','SMARCA5','MND1','PAPD7','ERCC8','CDK7','POLK','MSH3','XRCC4',
    'RAD50','HUS1B','MDC1','FANCE','POLH','MCM3','MMS22L','REV3L','HDAC2',
    'MCM9','SHPRH','AP5Z1','RPA3','POLM','HUS1','RFC2','NAMPT','CDK5','XRCC2',
    'ESCO2','WRN','GINS4','POLB','SPIDR','PRKDC','NBN','RAD54B','RRM2B','RAD21',
    'NSMCE2','TONSL','RECQL4','SMARCA2','FANCG','SMC5','FANCC','CDC14B',
    'RAD23B','INIP','SWI5','ABL1','IPMK','PTEN','HELLS','SFR1','SMC3','RRM1',
    'FANCF','SSRP1','DDB1','FEN1','KAT5','MUS81','POLD3','MRE11A','ATM','H2AFX',
    'CHEK1','RAD52','RAD51AP1','NABP2','TIMELESS','NAP1L1','UBE2N','UNG',
    'GTF2H3','BRCA2','LIG4','APEX1','REC8','G2E3','FANCM','MNAT1','ZFYVE26',
    'RAD51B','MLH3','YY1','XRCC3','FAN1','RAD51','INO80','TP53BP1','DUT',
    'MORF4L1','FANCI','BLM','EME2','MEIOB','DNASE1L2','SLX4','USP7','ERCC4',
    'SMG1','PALB2','NSMCE1','SLX1B','SLX1A','USP10','GINS2','FANCA','RPA1',
    'ZSWIM7','TOP3A','ATAD5','LIG3','RAD51D','CDK12','PSMC3IP','BRCA1','EME1',
    'RAD51C','BRIP1','ESCO1','RBBP8','MUM1','XAB2','SWSAP1','C19orf40','ERCC1',
    'LIG1','CCDC155','PNKP','RAD21L1','AP5S1','MCM8','SPO11','CDC45','CHEK2',
    'POLR2F','DMC1','XRCC6','MAPK12','FANCB','UBA1','TEX11','UBE2A','BRCC3') %>%
    return()
}
#' load genes of germline TOP-ART criterion as character vector
#' @export
#' @return A character vector.
load_germline_topart_genes <- function(){
  c('FANCL','ERCC3','BARD1','FANCD2','BAP1','ATR','FAM175A','RAD50','FANCE',
    'XRCC2','WRN','NBN','RECQL4','FANCG','FANCC','PTEN','FANCF','MRE11A','ATM',
    'BRCA2','MLH3','RAD51','FANCI','BLM','SLX4','ERCC4','PALB2','FANCA',
    'RAD51D','BRCA1','RAD51C','BRIP1','CHEK2','FANCB') %>%
    return()
}
#' run YAPSA analysis of mutational signatures
#'
#' @param FILE_REFGEN indexed fasta file of reference genome .fa 
#' @param SEQ_METHOD character either "WES" or "WGS".
#' @param TBL_SOM_SNV_YAPSA tibble/dataframe containing relevant SNVs for YAPSA
#' analysis.
#' @return A tibble containing the results for each called signature.
#' @export
#' @import YAPSA
#' @importFrom dplyr as_tibble
run_YAPSA_signature_analysis <- function(FILE_REFGEN, SEQ_METHOD, TBL_SOM_SNV_YAPSA){
  #require(YAPSA)
  # load data for signatures and cutoffs from the package
  data(sigs)
  data(cutoffs)
  # create the SNV mutational cataloge
  mutation_catalogue_list <- create_mutation_catalogue_from_df(
    this_df = TBL_SOM_SNV_YAPSA,
    this_seqnames.field = "CHROM",
    this_refGenome = Rsamtools::FaFile(FILE_REFGEN),
    this_wordLength = 3,
    this_rownames = rownames(AlexCosmicValid_sig_df),
    this_verbose = 0)
  mutation_catalogue_df <- as_tibble(mutation_catalogue_list$matrix)
  # perform signature anaylysis
  if(SEQ_METHOD == "WES"){
    sample_col <- "Valid_norm"
    cor_list <- targetCapture_cor_factors[[targetCapture]]
    corrected_catalogue_df <- normalizeMotifs_otherRownames(mutation_catalogue_df,
                                                            cor_list$rel_cor)
    CosmicValid_LCDlist <- LCD_complex_cutoff(
      in_mutation_catalogue_df = corrected_catalogue_df,
      in_signatures_df = AlexCosmicValid_sig_df,
      in_cutoff_vector = cutoffCosmicValid_rel_df[6, ],
      in_filename = NULL,
      in_method = "abs",
      in_sig_ind = AlexCosmicValid_sigInd_df)
  } else {
    sample_col <- "Valid_abs"
    corrected_catalogue_df <- mutation_catalogue_df
    CosmicValid_LCDlist <- LCD_complex_cutoff(
      in_mutation_catalogue_df = corrected_catalogue_df,
      in_signatures_df = AlexCosmicValid_sig_df,
      in_cutoff_vector = cutoffCosmicValid_abs_df[6, ],
      in_filename = NULL,
      in_method = "abs",
      in_sig_ind = AlexCosmicValid_sigInd_df)
  }
  ## create confidence intervals
  tbl_mutsig_conf <- variateExp(
    in_catalogue_df = corrected_catalogue_df,
    in_sig_df = AlexCosmicValid_sig_df[, rownames(CosmicValid_LCDlist$exposures), drop = FALSE],
    in_exposures_df = CosmicValid_LCDlist$exposures,
    in_sigLevel = 0.025, in_delta = 0.4) %>%
    mutate(sample=sample_col)
  return(tbl_mutsig_conf)
}
#' calculate signature score
#' @import dplyr
#' @export
estimate_signature_score <- function(TBL_YAPSA=NULL, 
                                     SEQ_METHOD,
                                     run_yapsa=F,
                                     FILE_REFGEN=NULL,
                                     TBL_SOM_SNV_YAPSA=NULL,
                                     TOTAL_SNVS_YAPSA=NULL){
  if(run_yapsa==T){
    if(nrow(TBL_SOM_SNV_YAPSA)==0){
      ## if no somatic small variants are provided... assuming no variants present
      ## signature score --> 0
      meta_info <- list(AC3_exp=0, AC3_conf=NA, 
                        signature_score=0,
                        total_snvs_yapsa=0,
                        P1_message="no input snvs and no YAPSA result provided -> assuming there are no somatic SNVs -> 0 pts") %>%
        return()
    } else {
      TBL_YAPSA <- run_YAPSA_signature_analysis(FILE_REFGEN, SEQ_METHOD, TBL_SOM_SNV_YAPSA)
    }
  }
  #rel_norm_type <- ifelse(SEQ_METHOD=='WGS', 'Valid_abs', 'Valid_norm')
  df_signatures <- TBL_YAPSA %>%
    filter(#norm_type==rel_norm_type,
      !exposure==0,
      sig %in% paste0("AC", c(1:30))) %>% 
    rowwise() %>%
    mutate(conf=ifelse(between(0, lower, upper), "low", "high")) %>%
    select(sig, conf, exposure)
  if(is.null(TOTAL_SNVS_YAPSA)){
    TOTAL_SNVS_YAPSA <- df_signatures %>% 
      pull(exposure) %>%
      sum()
  }
  AC3_data <- df_signatures[which(df_signatures$sig=="AC3"),]
  if(nrow(AC3_data)==0){
    signature_score <- 0
    message <- "AC3 not detected"
  } else if(TOTAL_SNVS_YAPSA<25){
    signature_score <- 0
    message <- "number of SNVs not sufficient (below 25)"
  } else if(AC3_data$conf=="low"){
    signature_score <- 1
    message <- "AC3 detected -> zero in confidence interval"
  } else {
    signature_score <- 2
    message <- "AC3 detected -> zero not in confidence interval"
  }
  if(signature_score==0){
    meta_info <- list(AC3_exp=0, AC3_conf=NA, 
                      signature_score=0)
  } else {
    meta_info <- AC3_data %>% 
      select(AC3_exp=exposure, AC3_conf=conf) %>%
      as.list() %>%
      append(list(signature_score=signature_score))
  }
  return(meta_info %>% 
           append(list(P1_message=paste("mutational signature analysis based on",
                                        TOTAL_SNVS_YAPSA,"SNVs ->",
                                        message, "->", signature_score, "pts"),
                       total_snvs_yapsa=TOTAL_SNVS_YAPSA
           )))
}
#' @return A List.
#' @export
#' @importFrom dplyr case_when
estimate_rearrangement_score <- function(HRD, LST){
  rearrangement_score <- case_when(
    HRD+LST>20 ~ 2,
    HRD+LST>10 ~ 1,
    TRUE ~ 0
  )
  message = paste("HRD:", HRD, "LST:", LST, "->", rearrangement_score, "pts")
  return(list(rearrangement_score=rearrangement_score, P2_message=message, HRD=HRD, LST=LST))
}
#' @return A List.
#' @import GenomicRanges
#' @import dplyr
#' @export
estimate_germline_score <- function(GR_GERM_SMALL_VARS, 
                                    GERM_TA_GENES){
  if(is.null(GR_GERM_SMALL_VARS)){
    germline_score <- 0
    gr_germ <- NULL
    G2_message=NULL
  } else {
    if(length(GenomicRanges::ranges(GR_GERM_SMALL_VARS))==0){  
      germline_score <- 0
      gr_germ <- NULL
      G2_message=NULL
    } else {
      gr_germ_raw <- GR_GERM_SMALL_VARS %>%
        .[(GenomicRanges::elementMetadata(.)[,"ACMG_class"] %in% c("4","5"))] %>%
        .[(GenomicRanges::elementMetadata(.)[,"gene"] %in% GERM_TA_GENES)]
      if(length(GenomicRanges::ranges(gr_germ_raw))==0){
        germline_score <- 0
        gr_germ <- NULL
        G2_message <- NULL
      } else {
        germline_score <- 1
        G2_message <- gr_germ_raw %>% 
          as_tibble() %>%
          select(relevant_genes=gene,
                 ACMG_class)
        gr_germ <- gr_germ_raw
      }  
    }
  }  
  return(list(germline_score=germline_score,
              gr_germ=gr_germ,
              G2_message=G2_message
  ))
}
#' @return A List.
#' @import dplyr
#' @export
estimate_somatic_score <- function(GR_CNV, 
                                   GR_SOM_SMALL_VARS, 
                                   PURITY, 
                                   SEX, 
                                   gr_germ, 
                                   DEBUG, 
                                   SOM_TA_GENES,
                                   GR_SOM_TA_GENES,
                                   FILE_BAM,
                                   FILE_RNA_BAM,
                                   PLOIDY,
                                   FILE_VCF,
                                   showReadDetail){
  # gr_germ <- g2_criterion$gr_germ
  # purity=PURITY
  # ploidy = PLOIDY
  # sex=SEX
  # somCna=GR_CNV
  # somSmallVars=GR_SOM_SMALL_VARS
  # germSmallVars=gr_germ
  # geneModel=GR_SOM_TA_GENES
  # bamDna=FILE_BAM
  # bamRna=FILE_RNA_BAM
  # includeHomoDel=TRUE
  # includeIncompleteDel=FALSE
  # showReadDetail=FALSE
  # byTcn=FALSE
  # vcf=raw_vcf
  # assumeSomCnaGaps=FALSE
  # colnameTcn=NULL
  # colnameCnaType=NULL
  # printLog = T
  # distCutOff <- 5000
  detailed_genewise_output <- capture.output(
    pred_zyg <- predict_zygosity(
      purity=PURITY, 
      ploidy = PLOIDY,
      sex=SEX,
      somCna=GR_CNV, 
      somSmallVars=GR_SOM_SMALL_VARS, 
      germSmallVars=gr_germ, 
      geneModel=GR_SOM_TA_GENES,
      bamDna=FILE_BAM,
      bamRna=FILE_RNA_BAM,
      includeHomoDel=TRUE,
      includeIncompleteDel=FALSE,
      showReadDetail=showReadDetail,
      byTcn=FALSE,
      vcf=FILE_VCF,
      assumeSomCnaGaps=TRUE,
      colnameTcn=NULL,
      colnameCnaType=NULL,
      printLog = TRUE
      
    )
  )
  ## define somatic score
  if(is.null(pred_zyg$eval_per_gene)){
    ## if no gene has any affections
    somatic_score=0
    corr_eval_per_gene <- NULL
  } else {
    
    var_origin <- pred_zyg$eval_per_variant %>%
      group_by(gene) %>%
      summarize(is_som=paste(origin, collapse=", ") %>% str_match("somatic") %>% as.character())
    
    corr_eval_per_gene <- pred_zyg$eval_per_gene %>%
      left_join(var_origin, by="gene") %>%
      mutate(score=case_when(
        status == "all_copies_affected" ~ 2,
        is.na(is_som) ~ 0,
        TRUE ~ 1
      )) %>% select(-is_som)
    somatic_score <- max(corr_eval_per_gene$score)
  }
  
  
  #   if(str_detect(paste(pred_zyg$eval_per_gene$status, collapse = " "), 
  #                      "all_copies_affected")){
  #   ## if any gene has all affected
  #   somatic_score=2
  # } else {
  #   ## now germline score 1 need to be filtered
  #   def <- pred_zyg$eval_per_variant %>%
  #     select(gene, origin) %>%
  #     left_join(pred_zyg$eval_per_gene %>% select(gene, status), by="gene") %>%
  #     ## filter out genes that are not in final eval... (lost in tumor)
  #     filter(!is.na(status)) %>%
  #     filter(origin=="somatic")
  #   if(nrow(def)==0){
  #     somatic_score <- 0
  #   } else {
  #     somatic_score <- 1
  #   }
  # }
  return(list(somatic_score=somatic_score, 
              G1_tbl_eval_per_variant=pred_zyg$eval_per_variant,
              G1_tbl_eval_per_gene=corr_eval_per_gene,
              G1_tbl_phasing_info=pred_zyg$phasing_info,
              uncovered_variants=pred_zyg$uncovered_input,
              ext_snp_phasing=pred_zyg$ext_snp_phasing,
              read_detail=pred_zyg$readpair_info,
              detailed_genewise_output=detailed_genewise_output)
  )
} 


# G1_tbl_eval_per_variant=g1_criterion$G1_tbl_eval_per_variant
# G1_tbl_eval_per_gene=g1_criterion$G1_tbl_eval_per_gene
# G1_phasing_info=g1_criterion$G1_tbl_phasing_info
# G2_message=g2_criterion$G2_message
#' @keywords internal
#' @importFrom knitr kable
#' @import dplyr
#' @import stringr
create_output_message <- function(overview, G1_tbl_eval_per_variant, G1_tbl_eval_per_gene, 
                                  G2_message, G1_phasing_info){
  header <- overview[1:5] %>% 
    gather() %>% 
    mutate_all(.funs = insert_tabs) %>%
    mutate(mes= paste(key, value, sep="  ")) %>% 
    pull(mes) %>% 
    paste(collapse="\n") 
  P1_out=paste0("signature score: ", overview$signature_score,
                "\nscoring:\n  ", overview$P1_message) %>% str_replace_all(" ->", "\n  ->")
  P2_out=paste0("rearrangement score: ", overview$rearrangement_score,
                "\nscoring:\n  ", overview$P2_message) %>% str_replace_all(" ->", "\n  ->")
  G2_out1 <- paste0("germline score: ", overview$germline_score) 
  G1_out1 <- paste0("somatic score: ", overview$somatic_score,"\nevaluation per variant:\n  ")
  if(!is.null(G2_message)){
    G2_out2 <- paste(knitr::kable(head(G2_message, n=nrow(G2_message))),collapse = "\n  ") %>% paste0("\nscoring:\n  ", .)
  } else {
    G2_out2 <- NULL
  }
  if(!is.null(G1_tbl_eval_per_gene)){
    if(!is.null(G1_phasing_info)){
      G1_phasing_info <- G1_phasing_info %>% dplyr::rename(wt_cp=left_wt_cp, no_ovlp=no_overlap, none_rw=none_raw, classes=class_comb) %>%
        relocate(1:13, DNA_rds, RNA_rds, classes, info)
      #  if(nrow(G1_phasing_info)!=0){
      phasing_list <-lapply(unique(G1_phasing_info$gene), function(GENE){
        tbl_gene <- G1_phasing_info %>% filter(gene==GENE) %>% #select(-gene) %>%
          apply(.,1,paste, collapse="\t", simplify=F) %>% 
          unlist() %>%
          c(paste(names(G1_phasing_info) %>% .[which(.!="name")],collapse="\t"),.) %>%
          paste(collapse="\n\t") %>%
          ## change things in final message
          str_replace_all("wt-copies left", "wt-cp") %>%
          str_replace_all("muts on diff reads", "at diff")%>% 
          str_replace_all("unclear: phasing -> ", "")%>%
          str_replace_all("unclear: ", "")%>%
          silver()
        return(c(gene=GENE, phasing_info=tbl_gene))
      }) %>% 
        bind_rows()  
    } else {
      phasing_list <- tibble(gene=NA, phasing_info=NA)
    }
    sorting <- G1_tbl_eval_per_gene %>%
      arrange(desc(score)) %>%
      pull(gene) %>% c(.,G1_tbl_eval_per_variant$gene) %>% unique()
    G1_out3 <- knitr::kable(head(G1_tbl_eval_per_gene %>% mutate(gene=factor(gene, levels=sorting)) %>% 
                                   arrange(gene), n=nrow(G1_tbl_eval_per_gene))) %>% 
      as.list() %>% unlist() %>%
      tibble(gene_info=.) %>%
      mutate(gene=str_split(gene_info, "\\|") %>% map_chr(.,2) %>% str_replace_all(" ","")) %>%
      left_join(phasing_list) %>%
      select(-gene) %>%
      mutate(
        gene_info=ifelse(str_detect(gene_info,"unclear"),
                         paste0("\033[31m", gene_info, "\033[39m"),
                         gene_info),
        mes=ifelse(!is.na(phasing_info),
                   paste(gene_info, phasing_info, sep="\n\t"),
                   paste0(gene_info,""))) %>%
      pull(mes) %>%
      paste0(collapse="\n  ") %>% 
      paste0("\n\nevaluation per gene:\n  ",.)
  } else {
    G1_out3 <- NULL
    sorting <- NULL
  }
  if(!is.null(G1_tbl_eval_per_variant)){
    if(is.null(sorting)){
      sorting <- G1_tbl_eval_per_variant$gene
    }
    G1_out2 <- paste(knitr::kable(head(G1_tbl_eval_per_variant %>% 
                                         mutate(gene=factor(gene, levels=sorting)) %>%
                                         arrange(gene) %>%
                                         select(gene, class, chr, pos, af=af, tcn, cna_type, aff_cp, pre_info) %>%
                                         mutate(class=str_replace(class,"nonsynonymous", "nonsyn") %>% 
                                                  str_replace("frameshift", "fs") %>%
                                                  str_replace("deletion", "del") %>%
                                                  str_replace("insertion", "ins"),
                                                pre_info=str_replace(pre_info, "LOH detected", "LOH") %>% 
                                                  str_replace("somatic-variant", "som") %>%
                                                  str_replace("germline-variant", "germ") %>%
                                                  str_replace("left wt-copies", "wt-cp") %>%
                                                  str_replace("all copies affected", "all aff") %>%
                                                  str_replace("not all affected", "not all aff") %>%
                                                  str_replace("variant lost in tumor", "lost in tumor")
                                         ) %>%
                                         mutate_at(.vars=c("af", "tcn", "aff_cp"), .funs=as.numeric) %>%
                                         mutate_at(.vars=c("af", "tcn", "aff_cp"), .funs=round, digits=2), 
                                       n=nrow(G1_tbl_eval_per_variant))),collapse = "\n  ") %>%
      paste0("\nevaluation per variant:\n  ",.)
  } else {
    G1_out2 <- NULL
  }
  G1_out1 <- paste0("somatic score: ", overview$somatic_score)
  output_message_raw <- list(
    list(header, P1_out, P2_out,G2_out1) %>% compact() %>%
      unlist() %>%
      paste(collapse="\n\n"),   
    G2_out2, "\n\n",
    G1_out1,
    G1_out2,
    G1_out3
  ) %>% compact() %>% unlist() %>% paste(collapse="") %>%
    paste0(.,"\n")
  
  return(output_message_raw)
}
#' @export
insert_tabs <- function(COL){
  mx <- max(nchar(COL)) 
  lapply(COL, function(inp){
    paste0(inp, paste0(rep(" ", mx-nchar(inp)), collapse = ""))
  }) %>% c(recursive=T)
}
#' @keywords internal
create_gene_overview <- function(G1_tbl_eval_per_gene, stf){
  if(!is.null(G1_tbl_eval_per_gene)){
    G1_tbl_eval_per_gene %>% filter(score==stf) %>% pull(gene) %>% paste(collapse=", ") %>%
      na_if("") %>%
      return()
  } else {
    return(NA)
  }
}
#' @keywords internal
#' @import GenomicRanges
#' @import dplyr
check_input_data <- function(SAMPLE_ID=NULL,
                             SEQ_METHOD,
                             HRD,
                             LST,
                             PURITY,
                             SEX,
                             FILE_BAM,
                             FILE_RNA_BAM,
                             GR_CNV,               # meta cols required: tcn, cna_type
                             GR_SOM_SMALL_VARS=NULL,    # meta cols required: gene, af, ref, alt
                             GR_GERM_SMALL_VARS=NULL,   # meta cols required: gene, af, ref, alt, ACMG_class
                             GR_GENCODE_EXON, # meta cols required: gene
                             
                             TBL_YAPSA=NULL,
                             FILE_REFGEN=NULL,
                             TOTAL_SNVS_YAPSA=NULL,
                             GR_SOM_SNV_YAPSA=NULL,
                             
                             PRINT_OUTPUT=F,
                             STORE_OUTPUT=F,
                             OUTPUT_DIR=NULL,
                             MAN_SOM_GENES=NULL,
                             MAN_GERM_GENES=NULL,
                             DEBUG=F,
                             FILTER_SOM_SMALL_VARS){
  
  complete_input <- c("SEQ_METHOD", "HRD", "LST", "PURITY", "SEX", "FILE_BAM",
                      "GR_CNV", "GR_GENCODE_EXON", "TBL_YAPSA") %>%
    as.character() %>%
    lapply(function(x){exists(as.character(x))}) %>%
    tibble(layer=names(.), exists=.) %>%
    filter(exists==F)
  
  if(nrow(complete_input)!=0){
    message(paste("following inputs are required but missing:", paste(complete_input$layer, collapse = "; "), "\n  aborting..."))
    return(list(proceed=F))
  } else if(!file.exists(FILE_BAM)){
    ## muss ja net immer da sein... also könnte man auch so amchen dasses nur warnt
    message("Input FILE_BAM seems to not exists\n  aborting...")
    return(list(proceed=F))
    # } else if(!file.exists(FILE_RNA_BAM)){  
    #  hier auch noch sagen was passiert falls nicht da 
    #} else if(!i){
  } else {
    ## from here... all basically required inputs are given... now check if they are correct
    if(!SEQ_METHOD %in% c("WES", "WGS")){
      message("input SEQ_METHOD must be a character: either \'WES\' or \'WGS\'\n  aborting...")
      return(list(proceed=F))
    }
    
    if(is.na(as.numeric(HRD))){
      message(paste("input HRD must be numeric or a character that can be converted to numeric;\n  ", HRD, 
                    "can not be converted to numeric\n  aborting..."))
      return(list(proceed=F))
    }
    if(is.na(as.numeric(LST))){
      message(paste("input LST must be numeric or a character that can be converted to numeric;\n ", LST, 
                    "can not be converted to numeric\n  aborting..."))
      return(list(proceed=F))
    }
    allowed_sex <- c("male", "m", "female", "f") %>% c(.,toupper(.))
    if(!SEX %in% allowed_sex){
      message(paste("input SEX must be one of", paste(allowed_sex, collapse = "\', \'") %>% paste0("\'", ., "\'"), "\n  aborting..."))
      return(list(proceed=F))
    }
    if(is.na(as.numeric(PURITY))){
      message(paste("input PURITY must be numeric or a character that can be converted to numeric;\n  ", PURITY, 
                    "can not be converted to numeric\n  aborting..."))
      return(list(proceed=F))
    } else if(!between(PURITY, 0, 1)){
      message(paste("input PURITY must be a numerci value between 0 and 1;\n  abortin..."))
      return(list(proceed=F))
    }
    # if(!class(GR_CNV)[1]=="GRanges"){
    #   message(paste("input GR_CNV must be a GRanges object; given input appears to be:", paste(unlist(class(GR_CNV)), collapse = ";")))
    #   return(list(proceed=F))
    # } else if(!("tcn" %in% names(GenomicRanges::elementMetadata(GR_CNV))&"cna_type" %in% names(GenomicRanges::elementMetadata(GR_CNV)))){
    #   message("input GR_CNV requires the following metadata columns: \'tcn\' and \'cna_type\'\n  aborting...")
    #   return(list(proceed=F))
    # } else {
    #   if(sum(is.na(GenomicRanges::elementMetadata(GR_CNV)[,"tcn"]))>0){
    #     message(paste("warning: tcn column of input GR_CNV contains", 
    #                   sum(is.na(GenomicRanges::elementMetadata(GR_CNV)[,"tcn"])),
    #                   "NA values;\n  they will be taken as neutral TCN (tcn=2)"))
    #     GenomicRanges::elementMetadata(GR_CNV)[,"tcn"][which(is.na(GenomicRanges::elementMetadata(GR_CNV)[,"tcn"]))] <- 2
    #   }
    #   if(sum(is.na(GenomicRanges::elementMetadata(GR_CNV)[,"cna_type"]))>0){
    #     message(paste("warning: cna_type column of input GR_CNV contains", 
    #                   sum(is.na(GenomicRanges::elementMetadata(GR_CNV)[,"cna_type"])),
    #                   "NA values;\n  they will be taken as hetero-zygous"))
    #     GenomicRanges::elementMetadata(GR_CNV)[,"cna_type"][which(is.na(GenomicRanges::elementMetadata(GR_CNV)[,"cna_type"]))] <- "no CNA type found"
    #   }
    #   #new_GR_CNV <- insert_missing_cnv_regions(GR_CNV, GR_GENCODE_EXON, SEX)
    #   new_GR_CNV <- GR_CNV
    # }
    ### check somatic small variant input
    # if(!is.null(GR_SOM_SMALL_VARS)){
    #   if(!class(GR_SOM_SMALL_VARS)[1]=="GRanges"){
    #     message(paste("input GR_SOM_SMALL_VARS must be a GRanges object; given input appears to be:", 
    #                   paste(unlist(class(GR_SOM_SMALL_VARS)), collapse = ";")))
    #     return(list(proceed=F))
    #   } else if(!(("gene" %in% names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS))|
    #                "GENE" %in% names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)))&
    #               ("af"   %in% names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS))|
    #                "AF"   %in% names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)))&
    #               ("ref"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS))|
    #                "REF"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)))&
    #               ("alt"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS))|
    #                "ALT"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS))))){
    #     message("input GR_SOM_SMALL_VARS requires the following metadata columns: \'gene\'/\'GENE\', \'ref\'/\'REF\', \'alt\'/\'ALT\' and \'af\'/\'AF\'\n  aborting...")
    #     return(list(proceed=F))
    #   } else {
    #     col_gene <- str_match(names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)), "gene|GENE") %>%
    #       stri_remove_na()
    #     col_af <- str_match(names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)), "af|AF") %>%
    #       stri_remove_na()
    #     col_ref <- str_match(names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)), "ref|REF") %>%
    #       stri_remove_na()
    #     col_alt <- str_match(names(GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)), "alt|ALT") %>%
    #       stri_remove_na()
    #     GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)[,"gene"] <- GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)[,col_gene]
    #     GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)[,"maf"] <- GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)[,col_af]
    #     GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)[,"ref"] <- GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)[,col_ref]
    #     GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)[,"alt"] <- GenomicRanges::elementMetadata(GR_SOM_SMALL_VARS)[,col_alt]
    #   }
    # }
    # ### check germline small variant input
    # if(!is.null(GR_GERM_SMALL_VARS)){
    #   if(!class(GR_GERM_SMALL_VARS)[1]=="GRanges"){
    #     message(paste("input GR_GERM_SMALL_VARS must be a GRanges object; given input appears to be:", 
    #                   paste(unlist(class(GR_GERM_SMALL_VARS)), collapse = ";")))
    #     return(list(proceed=F))
    #   } else if(!(("gene" %in% names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS))|
    #                "GENE" %in% names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)))&
    #               "ACMG_class" %in% names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS))&
    #               ("af"   %in% names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS))|
    #                "AF"   %in% names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)))&
    #               ("ref"  %in% names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS))|
    #                "REF"  %in% names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)))&
    #               ("alt"  %in% names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS))|
    #                "ALT"  %in% names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS))))){
    #     message("input GR_GERM_SMALL_VARS requires the following metadata columns: \'gene\', \'ref\', \'alt\', \'ACMG_class\'  and \'af\'\n  aborting...")
    #     return(list(proceed=F))
    #   } else {
    #     col_gene <- str_match(names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)), "gene|GENE") %>%
    #       stri_remove_na()
    #     col_af <- str_match(names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)), "af|AF") %>%
    #       stri_remove_na()
    #     col_ref <- str_match(names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)), "ref|REF") %>%
    #       stri_remove_na()
    #     col_alt <- str_match(names(GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)), "alt|ALT") %>%
    #       stri_remove_na()
    #     num_ACMG_class <- GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)[,"ACMG_class"] %>%
    #       str_replace_all("Likely Pathogenic|likely pathogenic|LP", "4") %>%
    #       str_replace_all("Likely Benign|likely benign|LB", "2") %>%
    #       str_replace_all("Pathogenic|pathogenic|P", "5") %>%
    #       str_replace_all("Benign|benign|B", "1") %>%
    #       str_replace_all("Uncertain Significance|uncertain significance|US|VUS", "3")
    #     GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)[,"gene"] <- GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)[,col_gene]
    #     GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)[,"maf"] <- GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)[,col_af]
    #     GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)[,"ref"] <- GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)[,col_ref]
    #     GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)[,"alt"] <- GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)[,col_alt]
    #     GenomicRanges::elementMetadata(GR_GERM_SMALL_VARS)[,"ACMG_class"] <- num_ACMG_class
    #   }
    # }
    ################
    ## YAPSA check
    ################
    if(!is.null(TBL_YAPSA)){
      ## if a yapsa table is given
      required_cols <- c("sig", "exposure", "lower", "upper")
      tbl_check_presence <- tibble(col=required_cols) %>%
        rowwise() %>%
        mutate(exists=ifelse(col %in% names(TBL_YAPSA),
                             T, F)) %>%
        filter(exists==F)
      
      if(nrow(tbl_check_presence)!=0){
        message(paste("input TBL_YAPSA must contain the following columns:", paste(required_cols, collapse = "; "),
                      "\n the following are missing:", paste(tbl_check_presence$col, collapse = ", ")))
        
        ## check if enough alternative info is given to perform signature analysis
        run_yapsa=T
      } else {
        run_yapsa=F      
        if(is.null(TOTAL_SNVS_YAPSA)){
          message("input TOTAL_SNV_YAPSA is missing... taking sum of exposure column from TBL_YAPSA")
        }
        TBL_SOM_SNV_YAPSA <- NULL
      }
    } else {
      run_yapsa=T
    }
    if(run_yapsa==T){
      ## no yapsa table given
      if(is.null(FILE_REFGEN)){
        message("input FILE_REFGEN (indexed fasta file) is required to run mutational signature analysis")
        return(list(proceed=F)) 
      } else if(file.exists(FILE_REFGEN)){
        message(paste("input FILE_REFGEN (indexed fasta file) is required to run mutational signature analysis but seems to not exist:",
                      FILE_REFGEN))
        return(list(proceed=F)) 
      }
      if(is.null(GR_SOM_SNV_YAPSA)){
        message("input GR_SOM_SNV_YAPSA not provided")
        if(is.null(GR_SOM_SMALL_VARS)){
          message("\nand no somatic small variants provided (GR_SOM_SMALL_VARS) !! \n... assuming there are no somatic variants in the sample\n !! SIGNATURE SCORE WILL BE ZERO !!")
          TBL_SOM_SNV_YAPSA <- tibble()
        } else {
          message("extracting SNVs from input GR_SOM_SMALL_VARS")
          TBL_SOM_SNV_YAPSA <- as_tibble(GR_SOM_SMALL_VARS) %>%
            filter(width==1) %>%
            select(CHROM=seqnames, POS=start, REF=ref, ALT=alt)
        }
        
      } else {
        if(("ref"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA))|
            "REF"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA)))&
           ("alt"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA))|
            "ALT"  %in% names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA)))){
          col_ref <- str_match(names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA)), "ref|REF") %>%
            stri_remove_na()
          col_alt <- str_match(names(GenomicRanges::elementMetadata(GR_SOM_SNV_YAPSA)), "alt|ALT") %>%
            stri_remove_na()        
          TBL_SOM_SNV_YAPSA <- as_tibble(GR_SOM_SNV_YAPSA) %>%
            filter(width==1) %>%
            select(CHROM=seqnames, POS=start, REF=all_of(col_ref), ALT=all_of(col_alt))
        } else {
          message("input GR_SOM_SNV_YAPSA requires the following meta data columns: REF/ref; ALT/alt")
          return(list(proceed=F))
        }
        
        
      }
    }
    ##### ref genome
    if(!class(GR_GENCODE_EXON)[1]=="GRanges"){
      message(paste("input GR_GENCODE_EXON must be a GRanges object; given input appears to be:", 
                    paste(unlist(class(GR_GENCODE_EXON)), collapse = ";")))
      return(list(proceed=F))
    } else if(!("gene" %in% names(GenomicRanges::elementMetadata(GR_GENCODE_EXON)))){
      message("input GR_GENCODE_EXON requires the following metadata columns: \'gene\'\n  aborting...")
      return(list(proceed=F))
    }     
    ## if the input needs to be filtered for functional variants (exonic)
    if(FILTER_SOM_SMALL_VARS==T&!is.null(GR_SOM_SMALL_VARS)){
      GR_SOM_SMALL_VARS <- IRanges::subsetByOverlaps(GR_SOM_SMALL_VARS, 
                                                     GR_GENCODE_EXON)
    }
    ## warnings  
    if(is.null(SAMPLE_ID)&STORE_OUTPUT==T){
      cat("No sample identification; storing output with default sample ID\n  set SAMPLE_ID=\'your_sample_id\'")
    }
    if(is.null(GR_SOM_SMALL_VARS)){
      message("no somatic small variants provided (GR_SOM_SMALL_VARS=NULL); \n  assuming there are no relevant variants present")
    }
    if(is.null(GR_GERM_SMALL_VARS)){
      message("no germline variants provided (GR_GERM_SMALL_VARS=NULL); \n  assuming there are no relevant variants present")
    }
    if(is.null(MAN_SOM_GENES)){
      SOM_TA_GENES <- load_somatic_topart_genes()
    } else {
      #cat("Using manual list of genes for somatic score")
      SOM_TA_GENES <- MAN_SOM_GENES
    }
    if(is.null(MAN_GERM_GENES)){
      GERM_TA_GENES <- load_germline_topart_genes()
    } else {
      GERM_TA_GENES <- MAN_GERM_GENES
    }
    return(
      list(
        proceed=T,
        SOM_TA_GENES=        SOM_TA_GENES,
        GERM_TA_GENES=       GERM_TA_GENES,
        GR_CNV=              GR_CNV,
        GR_GERM_SMALL_VARS=  GR_GERM_SMALL_VARS,
        GR_SOM_SMALL_VARS=   GR_SOM_SMALL_VARS,
        run_yapsa=           run_yapsa,
        TBL_SOM_SNV_YAPSA=   TBL_SOM_SNV_YAPSA
      )
    )
  }  
}
#' calcuation of TOP-ART score
#' @import GenomicRanges
#' @import dplyr
#' @import readr
#' @export
calculate_TOP_ART_score <- function(SAMPLE_ID=NULL,
                                    SEQ_METHOD,
                                    HRD,
                                    LST,
                                    PURITY,
                                    PLOIDY,
                                    SEX,
                                    
                                    ## files
                                    FILE_BAM,
                                    FILE_RNA_BAM,
                                    FILE_VCF=NULL,
                                    
                                    ## granges
                                    GR_CNV,                   # meta cols required: tcn, cna_type
                                    GR_SOM_SMALL_VARS=NULL,   # meta cols required: gene, af, ref, alt
                                    GR_GERM_SMALL_VARS=NULL,  # meta cols required: gene, af, ref, alt, ACMG_class
                                    GR_GENCODE_EXON,     # meta cols required: gene
                                    GR_FULL_VARS=NULL,
                                    
                                    ## yapsa options
                                    TBL_YAPSA=NULL,
                                    FILE_REFGEN=NULL,
                                    TOTAL_SNVS_YAPSA=NULL,
                                    GR_SOM_SNV_YAPSA=NULL,
                                    
                                    ## options
                                    PRINT_OUTPUT=F,
                                    STORE_OUTPUT=F,
                                    OUTPUT_DIR=NULL,
                                    MAN_SOM_GENES=NULL,
                                    MAN_GERM_GENES=NULL,
                                    DEBUG=F,
                                    FILTER_SOM_SMALL_VARS=T,
                                    showReadDetail=F
){
  ## for manual runs
  # DEBUG=F
  # FILTER_SOM_SMALL_VARS=F
  # STORE_OUTPUT=F
  # MAN_SOM_GENES=NULL
  # MAN_GERM_GENES=NULL
  
  adapted_input <- check_input_data(SAMPLE_ID,
                                    SEQ_METHOD,
                                    HRD,
                                    LST,
                                    PURITY,
                                    SEX,
                                    FILE_BAM,
                                    FILE_RNA_BAM,
                                    GR_CNV,              
                                    GR_SOM_SMALL_VARS,   
                                    GR_GERM_SMALL_VARS,
                                    GR_GENCODE_EXON, 
                                    
                                    TBL_YAPSA,
                                    FILE_REFGEN,
                                    TOTAL_SNVS_YAPSA,
                                    GR_SOM_SNV_YAPSA,
                                    
                                    PRINT_OUTPUT,
                                    STORE_OUTPUT,
                                    OUTPUT_DIR,
                                    MAN_SOM_GENES,
                                    MAN_GERM_GENES,
                                    DEBUG,
                                    FILTER_SOM_SMALL_VARS)
  if(adapted_input$proceed==T){
    SOM_TA_GENES=adapted_input$SOM_TA_GENES
    GERM_TA_GENES=adapted_input$GERM_TA_GENES
    
    ## i think i dont need this because its done inside zygosity predictor but i keep it in any case
    
    ## dont remove thios otherwise it will not work for germine score
    GR_CNV=adapted_input$GR_CNV
    GR_GERM_SMALL_VARS=adapted_input$GR_GERM_SMALL_VARS
    GR_SOM_SMALL_VARS=adapted_input$GR_SOM_SMALL_VARS
    
    
    GR_SOM_TA_GENES <- GR_GENCODE_EXON %>%
      .[(GenomicRanges::elementMetadata(.)[,"gene"] %in% SOM_TA_GENES)]
    TBL_SOM_SNV_YAPSA <- adapted_input$TBL_SOM_SNV_YAPSA
    ## score calculation
    p1_criterion <- estimate_signature_score(TBL_YAPSA, 
                                             SEQ_METHOD,
                                             adapted_input$run_yapsa,
                                             FILE_REFGEN,
                                             TBL_SOM_SNV_YAPSA,
                                             TOTAL_SNVS_YAPSA)
    p2_criterion <- estimate_rearrangement_score(as.numeric(HRD), 
                                                 as.numeric(LST))
    g2_criterion <- estimate_germline_score(GR_GERM_SMALL_VARS, 
                                            GERM_TA_GENES)
    g1_criterion <- estimate_somatic_score(GR_CNV, 
                                           GR_SOM_SMALL_VARS, 
                                           PURITY, 
                                           SEX, 
                                           g2_criterion$gr_germ, 
                                           DEBUG,
                                           SOM_TA_GENES,
                                           GR_SOM_TA_GENES,
                                           FILE_BAM,
                                           FILE_RNA_BAM,
                                           PLOIDY,
                                           FILE_VCF,
                                           showReadDetail)
    total_score <- p1_criterion$signature_score+
      p2_criterion$rearrangement_score+
      g1_criterion$somatic_score+
      g2_criterion$germline_score
    ## output creation
    overview <- list(
      sample_id=ifelse(!is.null(SAMPLE_ID),
                       SAMPLE_ID, "default"),
      seq_method=SEQ_METHOD,
      sex=SEX,
      purity=PURITY,
      top_art_score=total_score,
      germline_score=g2_criterion$germline_score,
      germline_variants=g2_criterion$G2_message$relevant_genes %>% paste(collapse=", "),
      somatic_score=g1_criterion$somatic_score,
      genes_all_cp_aff=create_gene_overview(g1_criterion$G1_tbl_eval_per_gene, 2),
      genes_wt_cp_left=create_gene_overview(g1_criterion$G1_tbl_eval_per_gene, 1)
    ) %>%
      append(p1_criterion) %>%
      append(p2_criterion) %>%
      compact() %>%
      as_tibble()
    output_message <- create_output_message(overview, 
                                            g1_criterion$G1_tbl_eval_per_variant, 
                                            g1_criterion$G1_tbl_eval_per_gene, 
                                            g2_criterion$G2_message,
                                            g1_criterion$G1_tbl_phasing_info #%>%
                                            # dplyr::rename(wt_cp=left_wt_cp,
                                            #             DNA=DNA_rds,
                                            #            RNA=RNA_rds)
    )
    ### data export
    if(STORE_OUTPUT==T){
      if(is.null(OUTPUT_DIR)){
        message("variable OUTPUT_DIR required to store output")
      } else {
        write_tsv(overview, file=file.path(OUTPUT_DIR, "combined_criteria.tsv"))
        write_lines(output_message, 
                    file=file.path(OUTPUT_DIR, "output_message.txt"))
        if(!is.null(g1_criterion$G1_tbl_eval_per_variant)){
          write_tsv(g1_criterion$G1_tbl_eval_per_variant, 
                    file=file.path(OUTPUT_DIR, "G1_evaluation_per_variant.tsv"))      
        }
        if(!is.null(g1_criterion$G1_tbl_eval_per_gene)){
          write_tsv(g1_criterion$G1_tbl_eval_per_gene, 
                    file=file.path(OUTPUT_DIR, "G1_evaluation_per_gene.tsv"))     
        }
        if(!is.null(g1_criterion$G1_tbl_phasing_info)){
          write_tsv(g1_criterion$G1_tbl_phasing_info, 
                    file=file.path(OUTPUT_DIR, "G1_phasing_info.tsv")) 
        }
      }
    }
    if(PRINT_OUTPUT==T){
      cat(output_message)
    }
    return(list(total_score=total_score, 
                output_message=output_message,
                overview=overview,
                G1_eval_per_variant=g1_criterion$G1_tbl_eval_per_variant, 
                G1_eval_per_gene=g1_criterion$G1_tbl_eval_per_gene, 
                G2_message=g2_criterion$G2_message,
                G1_phasing_info=g1_criterion$G1_tbl_phasing_info,
                uncovered_input=g1_criterion$uncovered_variants,
                ext_snp_phasing=g1_criterion$ext_snp_phasing,
                read_detail=g1_criterion$read_detail,
                detailed_genewise_output=g1_criterion$detailed_genewise_output)
    )
  } else {
    cat("cannot calculate TOP-ART score")
    return(NULL)
  }
}