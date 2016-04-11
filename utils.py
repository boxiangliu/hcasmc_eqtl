import os,sys,shutil
import subprocess as sp
from time import strftime
from multiprocessing import Pool

def setwd(wd):
	os.chdir(wd)
	print "[ cwd ] " + os.getcwd()

def reportDone():
	name = sys.argv[0] + '.done'
	cmd = ['touch', name]
	sp.call(cmd)
	print "[ done ] created %s"%name

def call(cmd,dry_run=False):
	print "[ cmd ] %s"%cmd
	if not dry_run: sp.call(cmd, shell=True)

def readSampleList(sample_list):
	with open(sample_list, 'r') as f: 
		sample_list = [x for x in f.readlines()]
	sample_list = [x.strip() for x in sample_list] # remove white space.
	print "[ readSampleList ] samples are:"
	for sample in sample_list: print sample
	print "%i samples"%len(sample_list)
	return sample_list

def merge(bams, merged):
	bams = " ".join(bams)
	cmd = 'samtools merge %s %s'%(merged, bams)
	call(cmd)
	print "[ merge ] %s finished."%merged

def sort(bam):
	cmd = "samtools sort %s %s"%(bam, bam.replace('bam','sorted'))
	call(cmd)
	print "[ sort ] %s finished."%bam

def index(bam):
	cmd = "samtools index %s"%bam
	call(cmd)
	print "[ index ] %s finished."%bam

def compressVcf(sample_dir, sample='raw_variants.g.vcf'):
	cmd = "bgzip -c %s/%s > %s/%s"%(sample_dir, sample, sample_dir, sample.replace('vcf','vcf.gz'))
	call(cmd)
	print "[ compressVcf ] finished."

def indexVcf(sample_dir, sample='raw_variants.g.vcf.gz'):
	cmd = "tabix -p vcf %s/%s"%(sample_dir, sample)
	call(cmd)
	print "[ indexVcf] finished."

def printDateTime(pipe = 'stderr'):
	if pipe == 'stderr': 
		sys.stderr.write("-"*40 + strftime("%Y-%m-%d %H:%M:%S") + "-"*40)
	else: 
		sys.stdout.write("-"*40 + strftime("%Y-%m-%d %H:%M:%S") + "-"*40)

