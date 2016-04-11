import os, subprocess, shutil
from paramiko import SSHClient, AutoAddPolicy
from scp import SCPClient

def get_finished_jobs():
    os.chdir('/srv/gsfs0/projects/montgomery/bliu2/HCASMC_eQTL/data')
    cmd = 'ls */*.vcf.gz'
    rtn = subprocess.check_output(cmd, shell = True)
    print "%i jobs in total"%len(rtn.split('\n'))
    print rtn
    return rtn

def open_ssh_channel(hostname, username = 'bliu2', password = 'Lbx_911011'):
    # setup SSH and SCP: 
    ssh = SSHClient()
    ssh.load_system_host_keys()
    ssh.set_missing_host_key_policy(AutoAddPolicy())
    ssh.connect(hostname = hostname, username = username, password = password)
    scp = SCPClient(ssh.get_transport())
    return scp

def transfer(src, dst): 
    ''' src is a list of files to be transfered''' 
    scp.put(files = src, remote_path = dst)

# def transfer(rtn):

#     for line in rtn.split('\n'):
#         line = line.strip()
#         if line == "": continue
#         direct = line.split('/')[0]
#         print "transferring%s"%direct
#         scp.put(files = ["%s/%s"%(direct,'raw_variants.g.vcf.gz'),"%s/%s"%(direct,'raw_variants.g.vcf.gz.tbi'),"%s/%s"%(direct,'realignment_targets.list'),"%s/%s"%(direct,'recal_data.table'),"%s/%s"%(direct,'recal_reads.bai'),"%s/%s"%(direct,'recal_reads.bam')], remote_path = '/mnt/data/WGS_HCASMC/%s/'%direct)

def close_ssh_channel(ssh):
    ssh.close()



setwd(argv.sys[1])
sample_dirs = readSampleList(sys.argv[2])
ssh = open_ssh_channel('scg3.stanford.edu')
for sample_dir in sample_dirs: 
    src = "/mnt/data/WGS_HCASMC/%s/raw_variants.g.vcf.gz"%(sample_dir)
    dst = "/srv/gsfs0/projects/montgomery/bliu2/HCASMC_eQTL/data/%s/"%(sample_dir)
    if os.path.exists(dst): os.mkdir(dst)
    transfer(src, dst)
rtn = get_finished_jobs()
transfer(rtn)

                         