# Install ChIPseq pipeline: 
cd /usrs/bliu2
git clone https://github.com/kundajelab/TF_chipseq_pipeline
cd TF_chipseq_pipeline
./install_dependencies.sh 


# Install ChromHMM: 
wget --directory-prefix=/srv/persistent/bliu2/tools http://compbio.mit.edu/ChromHMM/ChromHMM.zip
cd /srv/persistent/bliu2/tools/
unzip ChromHMM.zip 
rm ChromHMM.zip