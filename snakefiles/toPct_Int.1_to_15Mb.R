args = commandArgs(TRUE)

sample.1_to_15Mb = args[1]
pop.size = as.numeric(args[2])
mdl = args[3]
thrhld = as.numeric(args[4])
chr_list = args[5]
# Vernot 2015: thrhld = 0.000316
# Sankararaman 2014: thrhld = 0.000113

suppressMessages(library("data.table"))
suppressMessages(library("dplyr"))

toPct_Int.1_to_15Mb.fn = function(sampled.data.1_to_15Mb, pop.size, Mdl, chr_list){
	 Pct_Intr.format = sampled.data.1_to_15Mb %>%
		mutate(Mdl=Mdl) %>%
		mutate(len_bp=(hapend-hapstrt)) %>%
		group_by(Mdl, Winsize, sim.tsk) %>%		# group by the Window size and the chr#
		summarise(sum_len_bp=sum(as.numeric(len_bp))) %>% 		# caclulate total introgressed based for each chr for a given window size
		select(Mdl, Winsize,sim.tsk,sum_len_bp) %>%
		mutate(Pct_Int=sum_len_bp/pop.size/(Winsize*1000000)) %>%
	  as.data.table()

	setkey(Pct_Intr.format, sim.tsk, Winsize, Mdl)
	Pct_Intr.format <- merge(dt.chr_list, Pct_Intr.format, all=TRUE, by="sim.tsk") %>%
	                   mutate(Pct_Int = ifelse(test = is.na(Pct_Int), yes = 0, no = Pct_Int)) %>%
	                   mutate(n1 = unique(sampled.data.1_to_15Mb$n1)) %>%
	                   mutate(n2 = unique(sampled.data.1_to_15Mb$n2)) %>%
	                   mutate(admix = paste(n1, n2, sep="_")) %>%
	                   as.data.table()

    	return(Pct_Intr.format)
}


prop_blw_thrhld.fn = function(data.1_to_15Mb.PctInt, THRHLD, Mdl){
	## First, test that if you filter the data by the THRHLD, some windows are actually below thrhld
	# If some windows in this data set are actually below the thrhld
	if(nrow(data.1_to_15Mb.PctInt %>% filter(Pct_Int<THRHLD))>0){
		data.1_to_15Mb.PctInt.prop_blw_thrhld = data.1_to_15Mb.PctInt %>%
		filter(!is.na(Mdl)) %>%
		filter(Pct_Int<THRHLD) %>%		# Find all chr for all winsizes where the pct_int is below the thrhld (i.e. the chr is depleted)
		mutate(tally=1) %>%			# add a flag
		group_by(Mdl, Winsize) %>%		# group Window sizes together
		summarise(sum_tally=sum(as.numeric(tally))) %>%	# count up how many chr for a given window size are depleted
		mutate(prop_blw_thrhld=sum_tally/n_chr)	%>% # turn that count into to proportion of the total number of simulated chromosomes
    as.data.table()


		setkey(data.1_to_15Mb.PctInt.prop_blw_thrhld, Winsize, sum_tally, prop_blw_thrhld, Mdl)
   		# Fill in any rows where the prop_blw_thrhld = 0
		#for(i in seq(5,10)){
		for(i in seq(1,10)){
		      if( nrow(filter(data.1_to_15Mb.PctInt.prop_blw_thrhld, Winsize==i))!=1 ){
		        data.1_to_15Mb.PctInt.prop_blw_thrhld = rbind(data.1_to_15Mb.PctInt.prop_blw_thrhld,
		                                                      data.table(Mdl=mdl,
		                                                                 Winsize=i,
		                                                                 sum_tally=0,
		                                                                 prop_blw_thrhld=0))
		      }
		}
		# Add the admix column
		data.1_to_15Mb.PctInt.prop_blw_thrhld <- data.1_to_15Mb.PctInt.prop_blw_thrhld %>%
		                                          mutate(n1 = unique(data.1_to_15Mb.PctInt$n1)) %>%
		                                          mutate(n2 = unique(data.1_to_15Mb.PctInt$n2)) %>%
		                                          mutate(admix = unique(data.1_to_15Mb.PctInt$admix))

		return(data.1_to_15Mb.PctInt.prop_blw_thrhld)
	}else{
		## If none of the windows in this data set are below the thrhld (i.e. none of the chr at any window size are depleted)
		# Write an empty data table where for each window size, report 0 chr are below thrhld
		data.empty = data.table()
		#for(i in seq(5,10)){
		for(i in seq(1,10)){
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

#######################
#######################

dt.sample.1_to_15Mb = fread(paste0('zcat ',sample.1_to_15Mb),
                            header=FALSE, verbose=FALSE, showProgress=FALSE,
                            col.names = c('sim.tsk','hapstrt','hapend', 'Winsize', 'n1', 'n2') )
dt.chr_list = fread(chr_list,
                    header=FALSE, verbose=FALSE, showProgress = FALSE,
                    col.names = c('sim.tsk'))

n_chr = nrow(dt.chr_list)

dt.sample.1_to_15Mb.PctInt = toPct_Int.1_to_15Mb.fn(dt.sample.1_to_15Mb, pop.size, mdl, dt.chr_list)

write('to PctInt complete', stderr())

#write.table(dt.sample.1_to_15Mb.PctInt, file="", quote=FALSE, sep='\t', row.names = FALSE, col.names=TRUE)

dt.sample.1_to_15Mb.PctInt.prop_blw_thrhld = prop_blw_thrhld.fn(dt.sample.1_to_15Mb.PctInt, thrhld, mdl)

write.table(dt.sample.1_to_15Mb.PctInt.prop_blw_thrhld, file="", quote=FALSE, sep='\t', row.names = FALSE, col.names=FALSE)

write('prop_blw_thrhld complete', stderr())
