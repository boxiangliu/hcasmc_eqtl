#!/usr/bin/bash
scripts=/srv/persistent/bliu2/HCASMC_eQTL/scripts
# Bash script to perform nominal pass for all tissues at every subsample size
for tissue in `cat /srv/persistent/bliu2/HCASMC_eQTL/data/gtex/gtex.v6p.eqtl.tissues.txt`; do
	bash $scripts/160816/nominal.pass.subsamples.v6p.core.sh ${tissue}
done
