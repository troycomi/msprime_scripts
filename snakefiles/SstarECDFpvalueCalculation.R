
suppressMessages(library(ggplot2))
suppressMessages(library(dplyr))
suppressMessages(library(data.table))
suppressMessages(library(tidyr))
suppressMessages(library(scales))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(stringr))
suppressMessages(library(optparse))

suppressMessages(require(bit64))
##########################

###########
###########
print('FUNCTION: Generate ECDFs')
generate.ecdf.region_ind.fn <- function(null.dt){
 print(' GENERATE ECDFS')
 #max_snps_ecdf <<- 0
 for( p in sort(unique(null.dt$pop))){
  assign(paste0('max_snps_ecdf','.',p), 0, inherits = TRUE)
  for( i in sort(unique(as.numeric(null.dt$n_region_ind_snps)))){
   if(i>0){
   print(i)
   nam <<- paste0('null.f.region_ind.', i, '.',p,'.ecdf')
    if(nrow(filter(null.dt, s_star>0, n_region_ind_snps==i, pop==p))>0){
     print(nam)
     assign(nam, ecdf(filter(null.dt, s_star>0, n_region_ind_snps==i, pop==p)$s_star), inherits = TRUE)
     assign(paste0('max_snps_ecdf','.',p), i, inherits = TRUE)
    }
   }
  }
 }
}
#############
#############
print('FUNCTION: Calculate S*-pvalue from ecdf')
estimate.pval.ecdf.region_ind.fn <- function(X){
  s_star <- as.numeric(X[["s_star"]])
  n_snps <- as.numeric(X[["n_region_ind_snps"]])
  pop <- X[["pop"]]
  max_snps <- eval(as.name(paste0('max_snps_ecdf.',pop)))
  if (n_snps==0){
    X[["sstarpval_region_ind_snps"]] <- NA
    }
  else if (n_snps<=max_snps) {
          if( exists(paste0("null.f.region_ind.",n_snps,".",pop,".ecdf")) ){
            ecdf.fn <- match.fun(paste0("null.f.region_ind.",n_snps,".",pop,".ecdf"))
            s_star_pval <- 1-ecdf.fn(s_star)
            X[["sstarpval_region_ind_snps"]] <- round(x = s_star_pval, digits = 4)
            }
  #        else if { }
    }
  else if (n_snps>max_snps) {
    ecdf.fn <- match.fun(paste0("null.f.region_ind.",max_snps,".",pop,".ecdf"))
    s_star_pval <- 1-ecdf.fn(s_star)
    X[["sstarpval_region_ind_snps"]] <- round(x = s_star_pval, digits = 4)
    }
  return(X[c("chrom","winstart","winend","ind_id","pop","s_star","n_region_ind_snps","sstarpval_region_ind_snps",
            "num_s_star_snps","hap_1_s_start","hap_1_s_end","hap_2_s_start","hap_2_s_end",
            "n_s_star_snps_hap1","n_s_star_snps_hap2","s_star_haps")])
}
#############
#############
print('FUNCTION: Write OUTPUT TABLES, FILTERED FOR S* AND MATCH PVALUES')
write.filtered.bed.fn <- function(dt, outputdir, mdl, admix, chrom, spval, matchpval){
 print(' Writing filtered .bed file')
 dat_1 <- dt %>%
     filter(s_star>0) %>%
     filter(sstarpval_region_ind_snps<=spval) %>%
     filter(match_pvalue<=matchpval) %>%
     filter(haplotype==1) %>%
     select(msp_ID, winstart, winend) %>%
     #select(msp_ID, hap_1_s_start, hap_1_s_end) %>%
     setnames(c('msp_ID','start','end')) %>%
     as.data.table()

 dat_2 <- dt %>%
     filter(s_star>0) %>%
     filter(sstarpval_region_ind_snps<=spval) %>%
     filter(match_pvalue<=matchpval) %>%
     filter(haplotype==2) %>%
     select(msp_ID, winstart, winend) %>%
     #select(msp_ID, hap_2_s_start, hap_2_s_end) %>%
     setnames(c('msp_ID','start','end')) %>%
     as.data.table()

 dat.bed <- rbind(dat_1,dat_2)

 options(scipen=10)
 dat.bed %<>% mutate(start=str_trim(as.character(start),side = "both")) %>% mutate(end=str_trim(as.character(end), side = "both"))
 write.table(x = dat.bed,
           file = paste0(outputdir,'/',mdl,'_',chrom,'_',admix,'.sstar_sig_',spval,'.match_sig_N_',matchpval,'.isc_0','.bed'),
           quote = FALSE,
           sep = '\t',
           row.names = FALSE,
           col.names = TRUE)
 options(scipen=0)
}
#############
#############
print('FUNCTION: Write OUTPUT TABLES, ALL S* AND MATCH PVALUES')
write.all.bed.fn <- function(dt, outputdir, mdl, admix, chrom){
 print(' Writing unfiltered .bed file')
 dat.bed <- rbind(dt %>%
     			filter(s_star>0) %>%
     			filter(haplotype==1) %>%
     			#select(msp_ID, winstart, winend, sstarpval_region_ind_snps) %>%
     			select(msp_ID, winstart, winend, sstarpval_region_ind_snps, match_pvalue) %>%
     			#select(msp_ID, hap_1_s_start, hap_1_s_end, sstarpval_region_ind_snps, match_pvalue) %>%
     			setnames(c('msp_ID','start','end', 'sstarpval_region_ind_snps', 'match_pvalue')) %>%
                #setnames(c('msp_ID','start','end', 'sstarpval_region_ind_snps')) %>%
     			as.data.table(),
		  dt %>%
     			filter(s_star>0) %>%
     			filter(haplotype==2) %>%
     			#select(msp_ID, winstart, winend, sstarpval_region_ind_snps) %>%
     			select(msp_ID, winstart, winend, sstarpval_region_ind_snps,match_pvalue) %>%
     			#select(msp_ID, hap_2_s_start, hap_2_s_end, sstarpval_region_ind_snps,match_pvalue) %>%
     			setnames(c('msp_ID','start','end', 'sstarpval_region_ind_snps', 'match_pvalue')) %>%
                #setnames(c('msp_ID','start','end', 'sstarpval_region_ind_snps')) %>%
     			as.data.table()
		)

 options(scipen=10)
 dat.bed %<>% mutate(start=str_trim(as.character(start),side = "both")) %>% mutate(end=str_trim(as.character(end), side = "both"))
 write.table(x = dat.bed,
           file = paste0(outputdir,'/',mdl,'_',chrom,'_',admix,'.sstar_sig_','ALL','.match_sig_N_MH_','ALL','.isc_0','.bed'),
           quote = FALSE,
           sep = '\t',
           row.names = FALSE,
           col.names = TRUE)
 options(scipen=0)
}


#####################################

commandline.arguments.fn <- function(){
print(' Import command line arguments')
option_list = list(
    make_option(c("--inputdir"), action="store", default=NA, type='character', help="Directory containing input files"),
    make_option(c("--outputdir"), action="store", default=NA, type='character', help="Set output directory for bedfiles"),
    make_option(c("--mdl"), action="store", default=NA, type='character', help="Model type, e.g. Tenn_nonAfr"),
    make_option(c("--null_dir"), action="store", default='null', type='character', help="Null directory name"),
    make_option(c("--null_tag"), action="store", default='n1_0.0_n2_0.0', type='character', help="Null model tag, e.g. n1_0.0_n2_0.0"),
    make_option(c("--ecdf"), action="store", default=NA, type='character', help='Specify stored ECDF RData set to use'),
    make_option(c("--admix_dir"), action="store", default=NA, type='character', help="Admix directory name"),
    make_option(c("--admix_tag"), action="store", default=NA, type='character', help="Admix model tag, e.g. n1_0.02_n2_0.0"),
    make_option(c("--max_chrm_admix"), action="store", default=NA, type='numeric', help="Number of chromosomes to test from admix data"),
    make_option(c("--max_chrm_null"), action="store", default=NA, type='numeric', help="Number of chromosomes to test from null data"),
    make_option(c("--sstarpval"), action="store", default=0.01, type='numeric', help="Sstar pvalue cutoff for significance"),
    make_option(c("--matchpval"), action="store", default=0.05, type='numeric', help="Match pvalue cutoff for significance"),
    make_option(c("--nofilter"), action="store_true", default=TRUE, help="Print complete output, w/o filtering [default]"),
    make_option(c("--filter"), action="store_false", dest="nofilter", help="Print the filtered output"),
    make_option(c("--ecdf_only"), action="store_true", dest="ecdf_only", default=FALSE, help="Generate the ecdfs and then exit")
)

opt <<- parse_args(OptionParser(option_list=option_list))
}

# Read and assign the commandline arguments to variables
commandline.arguments.fn()
#####################################
#####################################

inputdir <- opt$inputdir
mdl <- opt$mdl
admix <- opt$null_tag
dir <- opt$null_dir
maxchrm <- opt$max_chrm_null

if ( opt$ecdf_only==TRUE ){ print("WARNING: Only Generating ecdfs") }

if ( file.exists(opt$ecdf) ){
    print( ' use specified ecdf file ')
    ecdf_data <- opt$ecdf
} else {
    ecdf_data <- paste0(inputdir,dir,'/SstarECDF_maxchrm_',maxchrm,'.RData.gz')
}

if (file.exists(ecdf_data)){
    print('LOAD ECDF DATA')
    print( ecdf_data )
    load(file = ecdf_data, verbose=TRUE)
    #print(paste0(' max_snps_ecdf: ', max_snps_ecdf))
} else {
    print('LOAD NULL DATA')
    sim_chrms <- fread(paste0('cat ',inputdir,dir,"/*.chr_list"))
    null.dt <- data.table(NULL)
    for( i in seq(1,as.numeric(maxchrm),by = 1)){
    #for( i in seq(1,nrow(sim_chrms),by = 1)){
        c <- sim_chrms[i][[1]]
        print(paste0(' Loading NULL chromosome number: ',c))
        infile <- paste0(inputdir,dir,'/RegionFiles/', c,".windowcalc.gz")
        dat <- fread(paste0('zcat ', infile), header=TRUE, select=c('s_star', 'n_region_ind_snps', 'pop'))
        dat <- filter(dat, s_star>0)
        null.dt <- rbind(null.dt, dat)
        remove(dat)
        gc()
    }
    ####################
    print(' Run generate.ecdf.fn null data')

    generate.ecdf.region_ind.fn(null.dt = null.dt)
    print(' ecdf generated')
    print(' remove null data')
    remove(null.dt)

     print(' save.image')
     save.image(file = ecdf_data, compress=TRUE, safe=TRUE)
     if ( opt$ecdf_only==TRUE ){
         print("Only Generating ecdfs")
         quit()
     }
}

## Redefine the commandline options here incase there are some conflicting ones already in the RData file we just loaded
commandline.arguments.fn()
####################

print('LOAD ADMIX DATA')

 # inputdir <- '/Genomics/akeylab/abwolf/SimulatedDemographic/Sstar/test/multi_sample/'
 # mdl <- 'Tenn_nonAfr'
 # admix <- 'n1_0.05_n2_0.0'
 # dir <- admix
 # maxchrm <- 5

inputdir <- opt$inputdir
mdl <- opt$mdl
admix <- opt$admix_tag
dir <- opt$admix_dir
maxchrm <- opt$max_chrm_admix

 sim_chrms <- fread(paste0('cat ',inputdir,dir,"/*.chr_list"))
 for( i in seq(1,as.numeric(maxchrm),by=1) ){
    out <- data.table(NULL)
    c <- sim_chrms[i][[1]]
    print(paste0(' Loading ADMIX chromosome number: ',c))
    infile <- paste0(inputdir,dir,'/RegionFiles/', mdl, "_",c,'_',admix,".windowcalc_out.gz")
    if( length(readLines(infile)) == 0 ) { next }

    dat <- fread(paste0('zcat ',infile), header=TRUE,
                select=c('chrom','winstart','winend','n_region_ind_snps',
                        'ind_id','pop','s_star','num_s_star_snps',
                        'hap_1_s_start','hap_1_s_end','hap_2_s_start','hap_2_s_end',
                        'n_s_star_snps_hap1','n_s_star_snps_hap2','s_star_haps'),
                na.strings=c("NA", "None",'.'))
                
    print(' Run estimate.pval.ecdf.fn')
    #out <- as.data.table(t(apply(X = dat,MARGIN = 1,FUN = estimate.pval.ecdf.region_ind.fn, max_snps=max_snps_ecdf)))
    out <- as.data.table(t(apply(X = dat,MARGIN = 1,FUN = estimate.pval.ecdf.region_ind.fn)))
    setnames(out,c("chrom","winstart","winend","ind_id","pop","s_star","n_region_ind_snps","sstarpval_region_ind_snps",
              "num_s_star_snps","hap_1_s_start","hap_1_s_end","hap_2_s_start","hap_2_s_end",
              "n_s_star_snps_hap1","n_s_star_snps_hap2","s_star_haps"))
    print(' estimate.pval.fn complete')

    out[,chrom:=as.numeric(chrom)]
    out[,winstart:=as.numeric(winstart)]
    out[,winend:=as.numeric(winend)]
    out[,ind_id:=as.character(ind_id)]
    out[,pop:=as.character(pop)]
    out[,s_star:=as.numeric(s_star)]
    out[,sstarpval_region_ind_snps:=as.numeric(sstarpval_region_ind_snps)]

    print(' Assign S* haplotype')

    req.snp.frac <- 0.8
    out$s_star_hap_1 <- (as.numeric(out$n_s_star_snps_hap1) / as.numeric(out$num_s_star_snps)) >= req.snp.frac
    out$s_star_hap_2 <- (as.numeric(out$n_s_star_snps_hap2) / as.numeric(out$num_s_star_snps)) >= req.snp.frac

    print(' Define haplotype sets')
    hap_1 = out[s_star_hap_1==TRUE & s_star_hap_2==FALSE]
    hap_2 = out[s_star_hap_1==FALSE & s_star_hap_2==TRUE]

    hap_1[,haplotype:=1]
    hap_2[,haplotype:=2]

    print(' rbind')
    out = rbind(
      hap_1,
      hap_2
    )

    out[,msp_ID:=paste0(ind_id,':',haplotype,'_',chrom)]
    out = out  %>% arrange(msp_ID, winstart, winend) %>% as.data.table()

    #######################

    print(paste0(' Loading MATCHPVAL chromosome number: ',c))
    infile <- paste0(inputdir,dir,'/match_pvalues/null-*/','pvalue_table_', mdl, "_",c,'_',admix,"_*.tsv.gz")
    admix.match_pvals <- fread(paste0('zcat ',infile), header=TRUE,
                                select=c('chr','start','end','haplotype','pvalue'),
                                col.names=c('chrom','winstart','winend','ID','match_pvalue'))
    admix.match_pvals[,msp_ID:=paste0(ID,'_',chrom)]
    admix.match_pvals = admix.match_pvals %>% arrange(msp_ID, winstart, winend) %>% as.data.table()


    print(' Assign matchpvals')
    out <- left_join(out, admix.match_pvals, by = c('msp_ID', 'winstart', 'winend')) %>% as.data.table()

    remove(admix.match_pvals)

    ######################

    print(' WRITE OUTPUT TABLES')
    print(paste0(' nofilter flag: ',opt$nofilter))

    if( opt$nofilter==FALSE ){
    write.filtered.bed.fn(dt = out,
    		outputdir = as.character(opt$outputdir),
    		mdl = opt$mdl,
            chrom = c,
    		admix = opt$admix_tag,
    		spval = opt$sstarpval,
    		matchpval = opt$matchpval)
    } else {
    write.all.bed.fn(dt = out,
    		outputdir = as.character(opt$outputdir),
    		mdl = opt$mdl,
            chrom = c,
    		admix = opt$admix_tag)
    }
 }

print(' FIN')
