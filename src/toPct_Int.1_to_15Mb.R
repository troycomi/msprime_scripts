args = commandArgs(TRUE)

sample.1_to_15Mb = args[1]
pop.size = as.numeric(args[2])
mdl = args[3]
thrhld = as.numeric(args[4])
n_chr = as.numeric(args[5])
# Vernot 2015: thrhld = 0.000316
# Sankararaman 2014: thrhld = 0.000113

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))

toPct_Int.1_to_15Mb.fn = function(sampled.data.1_to_15Mb, pop.size, Mdl, n_chr){
	 Pct_Intr.format = sampled.data.1_to_15Mb %>% 
		setnames(c('sim.tsk','hapstrt','hapend', 'sample_ID', 'Winsize', 'n1', 'n2')) %>%
		mutate(Mdl=Mdl) %>%
		select(-(sample_ID)) %>%
		mutate(len_bp=(hapend-hapstrt)) %>% 
		group_by(Mdl, Winsize, sim.tsk) %>%		# group by the Window size and the chr# 
		summarise(sum_len_bp=sum(as.numeric(len_bp))) %>% 		# caclulate total introgressed based for each chr for a given window size
		select(Mdl, Winsize,sim.tsk,sum_len_bp) %>% 
		mutate(Pct_Int=sum_len_bp/pop.size/(Winsize*1000000))
	
	Pct_Intr.format <- as.data.table(Pct_Intr.format)	
	setkey(Pct_Intr.format, sim.tsk, Winsize, Mdl)
  	Pct_Intr.format = Pct_Intr.format[CJ(rep(1:n_chr), unique(Winsize), unique(Mdl))]
  	Pct_Intr.format[is.na(Pct_Int), Pct_Int:=0]

	Pct_Intr.format$n1 <- sampled.data.1_to_15Mb$n1
	Pct_Intr.format$n2 <- sampled.data.1_to_15Mb$n2
	Pct_Intr.format$admix <- paste(Pct_Intr.format$n1, Pct_Intr.format$n2, sep="_")

	
    	return(Pct_Intr.format)
}


prop_blw_thrhld.fn = function(data.1_to_15Mb.PctInt, THRHLD, Mdl){
	## First, test that if you filter the data by the THRHLD, some windows are actually below thrhld
	# If some windows in this data set are actually below the thrhld
	if(nrow(data.1_to_15Mb.PctInt %>% filter(Pct_Int<THRHLD))>0){
		data.1_to_15Mb.PctInt.prop_blw_thrhld = data.1_to_15Mb.PctInt %>%
		filter(Pct_Int<THRHLD) %>%		# Find all chr for all winsizes where the pct_int is below the thrhld (i.e. the chr is depleted)
		mutate(tally=1) %>%			# add a flag
		group_by(Mdl, Winsize) %>%		# group Window sizes together
		summarise(sum_tally=sum(as.numeric(tally))) %>%	# count up how many chr for a given window size are depleted
		mutate(prop_blw_thrhld=sum_tally/n_chr)	# turn that count into to proportion of the total number of simulated chromosomes
		
		data.1_to_15Mb.PctInt.prop_blw_thrhld = as.data.table(data.1_to_15Mb.PctInt.prop_blw_thrhld)

		setkey(data.1_to_15Mb.PctInt.prop_blw_thrhld, Winsize, sum_tally, prop_blw_thrhld, Mdl)
   		# Fill in any rows where the prop_blw_thrhld = 0 
		for(i in seq(5,10)){
		      if( nrow(filter(data.1_to_15Mb.PctInt.prop_blw_thrhld, Winsize==i))!=1 ){
		        data.1_to_15Mb.PctInt.prop_blw_thrhld = rbind(data.1_to_15Mb.PctInt.prop_blw_thrhld,
		                                                      data.table(Mdl=mdl,
		                                                                 Winsize=i,
		                                                                 sum_tally=0,
		                                                                 prop_blw_thrhld=0))
		      }
		}
		# Add the admix column
		data.1_to_15Mb.PctInt.prop_blw_thrhld$n1 <- data.1_to_15Mb.PctInt$n1
		data.1_to_15Mb.PctInt.prop_blw_thrhld$n2 <- data.1_to_15Mb.PctInt$n2
		data.1_to_15Mb.PctInt.prop_blw_thrhld$admix <- data.1_to_15Mb.PctInt$admix

		return(data.1_to_15Mb.PctInt.prop_blw_thrhld)
	}else{
		## If none of the windows in this data set are below the thrhld (i.e. none of the chr at any window size are depleted)
		# Write an empty data table where for each window size, report 0 chr are below thrhld 
		data.empty = data.table()
		for(i in seq(5,10)){
			data.empty = rbind(data.empty, data.table(Mdl=Mdl, 
								  Winsize=i, 
								  sum_tally=0, 
								  prop_blw_thrhld=0, 
								  n1=unique(data.1_to_15Mb.PctInt$n1), 
								  n2=unique(data.1_to_15Mb.PctInt$n2),
								  admix=unique(data.1_to_15Mb.PctInt$admix)))
		}
		return(data.empty)
	}
}



dt.sample.1_to_15Mb = fread(sample.1_to_15Mb, header=FALSE, verbose=FALSE, showProgress=FALSE)

dt.sample.1_to_15Mb.PctInt = toPct_Int.1_to_15Mb.fn(dt.sample.1_to_15Mb, pop.size, mdl, n_chr)

write('to PctInt complete', stderr())

#write.table(dt.sample.1_to_15Mb.PctInt, file="", quote=FALSE, sep='\t', row.names = FALSE, col.names=TRUE)

dt.sample.1_to_15Mb.PctInt.prop_blw_thrhld = prop_blw_thrhld.fn(dt.sample.1_to_15Mb.PctInt, thrhld, mdl)

write.table(dt.sample.1_to_15Mb.PctInt.prop_blw_thrhld, file="", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)

write('prop_blw_thrhld complete', stderr())
