cd /srv/persistent/bliu2/tools/msCentipede
python call_binding.py --task learn test/Ctcf_chr10_motifs.txt.gz test/Gm12878_Rep1.bam test/Gm12878_Rep2.bam
python call_binding.py --task infer test/Ctcf_chr10_motifs.txt.gz test/Gm12878_Rep1.bam test/Gm12878_Rep2.bam

python plot_accessibility_profile.py test/Ctcf_chr10_motifs.txt.gz
python plot_accessibility_profile.py test/Ctcf_chr10_motifs_msCentipede_binding_posterior.txt.gz
