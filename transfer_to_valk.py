#!/usr/bin/env python
from utils import *
from paramiko import SSHClient, AutoAddPolicy
from scp import SCPClient


def transfer(samples, files, remote_path, hostname, username, password):
	# setup SSH and SCP: 
	ssh = SSHClient()
	ssh.load_system_host_keys()
	ssh.set_missing_host_key_policy(AutoAddPolicy())
	ssh.connect(hostname = hostname, username = username, password = password)
	scp = SCPClient(ssh.get_transport())

	for sample in samples:
		print "transferring %s"%sample
		for f in files:
			scp.put(files = ["%s/%s"%(sample,f)], remote_path = "%s/%s"%(remote_path,sample))

	ssh.close()

# def remove(rtn):
# 	for line in rtn.split('\n'):
# 		line = line.strip()
# 		if line == "": continue
# 		direct = line.split('/')[0]
# 		print "deleting %s..."%direct
# 		shutil.rmtree(direct)

if __name__ == "__main__":
	wd = sys.argv[1]
	sample_list = sys.argv[2]
	setwd(wd)
	samples = readSampleList(sample_list)
	files = ['raw_variants.g.vcf.gz','raw_variants.g.vcf.gz.tbi','realignment_targets.list','recal_data.table','recal_reads.bai','recal_reads.bam']
	remote_path = '/home/diskstation/wgs/WGS_HCASMC_BWA_mapping'
	transfer(samples, files, remote_path, "valkyr.stanford.edu", 'bosh','bosh')
	reportDone()