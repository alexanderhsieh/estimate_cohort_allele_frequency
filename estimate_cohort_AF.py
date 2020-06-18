## generates a table of <chr:pos:ref:alt> <number of samples carrying variant allele>
## estimated on the basis of genotypes

import sys
from optparse import OptionParser
import subprocess
import os
import gzip
import io

####################################################################################################
## handle arguments
####################################################################################################
parser = OptionParser()
parser.add_option('-i', '--input', dest='input_gvcf',help='input merged gvcf')
#parser.add_option('-m', '--min_vaf', dest='min_var',help='minimum VAF to be considered a carrier')
parser.add_option('-n', '--cohort_size', dest='cohort_size',help='number of samples in cohort')
parser.add_option('-o', '--output', dest='output_file',help='output tab-separated variants file')
(options, args) = parser.parse_args()

## check all arguments present
#if (options.input_gvcf == None or options.min_vaf == None or options.output_file):
if (options.input_gvcf == None or options.cohort_size == None or options.output_file == None):
	print('\n' + '## ERROR: missing arguments' + '\n')
	parser.print_help()
	print('\n')
	sys.exit()


input_gvcf = options.input_gvcf
#min_vaf = options.min_vaf
cohort_size = options.cohort_size
output_file = options.output_file


####################################################################################################
## parse gvcf
####################################################################################################
outf = open(output_file, 'w')
header = '\t'.join(['var_id', 'n_carriers', 'cohort_allele_frequency'])
outf.write(header + '\n')

carriers = {}
#with gzip.open(input_gvcf, 'r') as f:
with open(input_gvcf, 'r') as f:
	for line in f:
		tmp = line.strip().split('\t')
		if not line.startswith('##'):
			if tmp[0] == '#CHROM':
				idx = {col:index for index, col in enumerate(tmp)}
			else:
				chr, pos, ref, alt = tmp[idx['#CHROM']], tmp[idx['POS']], tmp[idx['REF']], tmp[idx['ALT']]
				
				format = tmp[idx['FORMAT']]
				
				key = '.:.:.:.'

				n = 0 ## counter for number of individuals carrying this variant with vaf > min_vaf

				## only process variant lines
				if not alt == '<NON_REF>':
					## need to deal with biallelic and multiallelic sites separately
					list_alt = alt.split(',')
					
					# biallelic case
					if len(list_alt) == 2:

						# get shortlist of sample columns where sample is NON ./. and shows evidence for alt allele
						shortlist = [s for s in tmp[idx['FORMAT']+1:] if not any(substring in s for substring in ['0/0', './.', '0|0', '.|.']) ]
						n = len(shortlist)

						freq = float(n)/float(cohort_size)

						alt_allele = alt.strip(',<NON_REF>').strip('<NON_REF>,')

						key = ':'.join([chr, pos, ref, alt_allele])

						## write out result
						carriers[key] = n
						outf.write(key + '\t' + str(n) + '\t' + str(freq) + '\n')

					# multiallelic case
					elif len(list_alt) > 2:
						
						for j in range(len(list_alt)):
							
							n = 0
							
							if not list_alt[j] == '<NON_REF>':
								
								if j == 0: # first alt allele
									shortlist = [s for s in tmp[idx['FORMAT']+1:] if any(substring in s for substring in ['/1', '|1', '1/', '1|']) ]
									
									alt_allele = list_alt[j]
									key = ':'.join([chr, pos, ref, alt_allele])
									n = len(shortlist)

									## write out result
									if n > 0:
										carriers[key] = n
										freq = float(n)/float(cohort_size)
										outf.write(key + '\t' + str(n) + '\t' + str(freq) + '\n')

								elif j > 0: # 
									shortlist = [s for s in tmp[idx['FORMAT']+1:] if any(substring in s for substring in ['/2', '|2', '2/', '2|', '/3', '|3', '3/', '3|', '/4', '|4', '4/', '4|']) ]


									alt_allele = list_alt[j]
									key = ':'.join([chr, pos, ref, alt_allele])
									n = len(shortlist)
								
									## write out result
									if n > 0:
										carriers[key] = n
										freq = float(n)/float(cohort_size)
										outf.write(key + '\t' + str(n) + '\t' + str(freq) + '\n')
						
						
	
					## testing
					'''
					if key == 'chr1:873548:C:T':
						print(key + '\t' + str(n))
						break
					'''


outf.close()

