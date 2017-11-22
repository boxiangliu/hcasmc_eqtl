Author: Mike Gloudemans

Existing HCASMC GWAS, eQTL, and sQTL files were formatted and prepared for colocalization analysis with the following command:

`bash hcasmc_data_prep.sh`

I then ran the FINEMAP colocalization pipeline on these files. 

In this folder, I have placed the colocalization pipeline config file `settings_used.config` along with the `git log` at
the time this pipeline was run (`git_status.txt`). If it is ever necessary to rerun the pipeline in exact the same way, we could do so by
checking out the appropriate version from the git repository, and running the old pipeline using the included config file.
