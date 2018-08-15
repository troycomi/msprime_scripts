#Converting 15Mb simulated bedfile of introgressed haplotyes to YMb coordinates. 
# File format
#   chr         hap.start   hap.end     len    Ne  ID
# sim_0.task_1	2048116	    4132969   2000000 
from __future__ import print_function
import sys
import copy
import gzip

max_len = int(sys.argv[3])  #specify the max simulated size

out_bed = gzip.open(sys.argv[2], 'wb') #create a new file to write the split haplotype locations to


#max_chrom_len = 1000000 
max_chrom_len = 5000000  #start at 5Mb window size, increase by 1Mb until reaching 10Mb size windows

while max_chrom_len <= max_len:
	open_bed = gzip.open(sys.argv[1], 'rb') #open the original introgressed haplotype bed file
	for line in open_bed:
		line_stripped = line.strip() #take a line from the original bed file, strip the end
    		line_list = line_stripped.split('\t') #split it into a list
    		if int(line_list[1]) < max_chrom_len: #if the haplotype start position is less that YMb
        		if int(line_list[2]) <= max_chrom_len: #and if the haplotype end position is equal or less than YMb
        			line_list.append(str(max_chrom_len/1000000)) #add the window size to the end of the haplotype list
				line_join = '\t'.join(line_list) #join the list to a string. tab-delimited
        			out_bed.write(line_join + '\n') #and write the string to the new bed file
        
        		elif int(line_list[2]) > max_chrom_len: #if the haplotype end is beyond YMb,
        			line_list[2] = str(max_chrom_len) #to the original line, set a hard end at YMb
        			line_list.append(str(max_chrom_len/1000000))
				line_join = '\t'.join(line_list)
        			out_bed.write(line_join + '\n')
	max_chrom_len += 1000000
				
            
open_bed.close()
out_bed.close()

print('fin',file=sys.stderr)
