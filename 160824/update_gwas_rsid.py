# change the rsid 
# for each line in gwas: 
# 	if markername starts with rs:
#		output line
#	else: 
#		parse line
#		get chr and position
# 		get rsid from vcf file using tabix 
#		example tabix ALL.chr3.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz 3:51616403-51616403 | cut -f1-3
#		if rsid is not empty:
#			update the rsid of gwas 
#			output the gwas line
#		else: 
#			report error 
#			# do not output gwas line
import sys, subprocess
n=0
for line in sys.stdin:
	n=n+1
	if n%100000==0:
		sys.stderr.write(str(n)+' lines processed.\n')
	if line.startswith('rs') or line.startswith('markername'):
		sys.stdout.write(line)
	else:
		split_line=line.strip().split()
		chrom=split_line[1]
		pos=split_line[2]
		cmd="tabix /srv/persistent/bliu2/shared/1000genomes/phase1v3/ALL.chr{chrom}.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.vcf.gz {chrom}:{pos}-{pos} | cut -f3-3".format(chrom=chrom,pos=pos)
		result=subprocess.check_output(cmd,shell=True)
		rsid=result.strip().split('\n')[-1]

		if rsid!=".":
			split_line[0]=rsid.strip()
			sys.stdout.write("\t".join(split_line)+'\n')
		else:
			sys.stderr.write('no rsid for {chrom}:{pos}\n'.format(chrom=chrom,pos=pos))
