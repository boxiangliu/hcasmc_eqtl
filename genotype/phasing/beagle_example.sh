#!/bin/bash
wd=/srv/persistent/bliu2/tools/beagle
java=/srv/persistent/bliu2/tools/jre1.8.0_91/bin/java
mkdir $wd
cd $wd

if [ ! -f beagle.03May16.862.jar ]; then
  echo
  echo "Downloading beagle.03May16.862.jar"
  wget http://faculty.washington.edu/browning/beagle/beagle.03May16.862.jar
fi

if [ ! -f bref.03May16.862.jar ]; then
  echo
  echo "Downloading bref.03May16.862.jar"
  wget http://faculty.washington.edu/browning/beagle/bref.03May16.862.jar
fi

echo

if [ ! -f test.03May16.862.vcf.gz ]; then
    echo
    echo "*** Downloading some 1000 Genomes Project data to file: test.03May16.862.vcf.gz ***"
    wget http://faculty.washington.edu/browning/beagle/test.03May16.862.vcf.gz
fi

echo
echo "*** Creating test files: ref.03May16.862.vcf.gz target.03May16.862.vcf.gz ***"
echo
zcat test.03May16.862.vcf.gz | cut -f1-190 | tr '/' '|' | gzip > ref.03May16.862.vcf.gz
zcat test.03May16.862.vcf.gz | cut -f1-9,191-200 | gzip > target.03May16.862.vcf.gz


echo
echo "*** Running test analysis with \"gt=\" argument ***"
echo
$java -jar beagle.03May16.862.jar gt=test.03May16.862.vcf.gz out=out.gt

echo
echo "*** Running test analysis with \"gl=\" argument ***"
echo
$java -jar beagle.03May16.862.jar gl=test.03May16.862.vcf.gz out=out.gl

echo
echo "*** Running test analysis with \"ref=\" argument ***"
echo
$java -jar beagle.03May16.862.jar ref=ref.03May16.862.vcf.gz gt=target.03May16.862.vcf.gz out=out.ref

echo
echo "*** Making \"bref\" file ***"
echo
$java -jar bref.03May16.862.jar ref.03May16.862.vcf.gz

echo
echo "*** Running test analysis with \"bref=\" argument ***"
echo
$java -jar beagle.03May16.862.jar ref=ref.03May16.862.bref gt=target.03May16.862.vcf.gz out=out.bref


echo 
echo "*** Truncate the reference panel ***"
echo 
zcat ref.03May16.862.vcf.gz | head -500 | gzip > ref.truncated.03May16.862.vcf.gz 


echo 
echo "*** Running test analysis with truncated reference ***"
echo 
$java -jar beagle.03May16.862.jar ref=ref.truncated.03May16.862.vcf.gz gt=target.03May16.862.vcf.gz out=out.ref.truncated


echo 
echo "*** modify target vcf ***"
echo 
zcat target.03May16.862.vcf.gz | awk 'BEGIN { OFS = "\t"} { if ($3=="rs138720731") {t=$4; $4=$5; $5=t; print;} else {print} }' | awk 'BEGIN { OFS = "\t"} { if ($3=="rs73387790") {$4="C"; $5="T"; print;} else {print} }' | gzip > target.mod.03May16.862.vcf.gz


echo 
echo "*** Test conform-gt ***"
echo 
wget https://faculty.washington.edu/browning/conform-gt/conform-gt.03May16.570.jar
$java -jar conform-gt.03May16.570.jar ref=ref.truncated.03May16.862.vcf.gz gt=target.mod.03May16.862.vcf.gz chrom=22 match=POS out=out.conform



echo 
