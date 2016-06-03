#!/bin/bash

#source /etc/profile

#module load r/3.0.1

TraitsFileName=$1
MaxFactorsN=15
MaxIterations=10000
BoundTol=0.001
VarTol=0.00001
e_pa=0.1
e_pb=10
a_pa=0.001
a_pb=0.1
OutDir=$2
tissue=hcasmc

echo "Rscript $scripts/get_peer_correction.Extended.R $TraitsFileName $MaxFactorsN $MaxIterations $BoundTol $VarTol $e_pa $e_pb $a_pa $a_pb $OutDir $tissue"
# R -f /users/joed3/goats/code/PEER/get_peer_correction.Extended.R --vanilla --slave --args $TraitsFileName $MaxFactorsN $MaxIterations $BoundTol $VarTol $e_pa $e_pb $a_pa $a_pb $OutDir $tissue
Rscript /srv/persistent/bliu2/HCASMC_eQTL/scripts/160527/get_peer_correction.Extended.R $TraitsFileName $MaxFactorsN $MaxIterations $BoundTol $VarTol $e_pa $e_pb $a_pa $a_pb $OutDir $tissue

echo "GetPeerExtended.sh DONE"


