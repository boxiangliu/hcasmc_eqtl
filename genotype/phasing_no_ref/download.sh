# download the beagle reference file: 
mkdir /srv/persistent/bliu2/tools/beagle/reference
cd /srv/persistent/bliu2/tools/beagle/reference
for i in $(seq 1 22) X; do
echo $i
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr$i.1kg.phase3.v5a.bref
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr$i.1kg.phase3.v5a.vcf.gz
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/individual_chromosomes/chr$i.1kg.phase3.v5a.vcf.gz.tbi
done

# download the genetic map: 
wget http://bochet.gcc.biostat.washington.edu/beagle/genetic_maps/plink.GRCh37.map.zip

# download panel file:
wget http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase3_v5a/integrated_call_samples_v3.20130502.ALL.panel 
