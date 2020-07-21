#!/usr/bin/env bash
#ln -s ~/Documents/GitHub/sigminer.wrapper/sigflow.R ~/.local/bin/sigflow

## Extraction PART:

#1
sigflow extract -i tcga_laml.maf.gz -o test_results/test_maf -m MAF -r 10 -T 4
if [ $? -ne 0 ]; then
    t1="failed at 'test_maf' in mode 'MAF'"
else
    t1="succeed in test #1"
fi

#2
sigflow extract -i tcga_laml.maf.gz -o test_results/test_maf_SBS -m SBS -r 10 -T 4
if [ $? -ne 0 ]; then
    t2="failed at 'test_maf' in mode 'SBS'"
else
    t2="succeed in test #2"
fi

#3
sigflow extract -i tcga_laml.maf.gz -o test_results/test_maf_ID -m ID -r 10 -T 4
if [ $? -ne 0 ]; then
    t3="failed at 'test_maf' in mode 'ID'"
else
    t3="succeed in test #3"
fi

#4
sigflow extract -i vcf/ -o test_results/test_vcf -m MAF -r 10 -T 4
if [ $? -ne 0 ]; then
    t4="failed at 'test_vcf' in mode 'MAF'"
else
    t4="succeed in test #4"
fi

$5
sigflow extract -i example_cn.tsv -o test_results/test_cn -m CN -r 10 -T 4
if [ $? -ne 0 ]; then
    t5="failed at 'test_cn' in mode 'CN'"
else
    t5="succeed in test #5"
fi


## Output test results
echo ===========================================
echo "             Test Results               "
echo $t1
echo $t2
echo $t3
echo $t4
echo $t5
