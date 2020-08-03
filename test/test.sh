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

#5
sigflow extract -i example_cn.tsv -o test_results/test_cn -m CN -r 10 -T 4
if [ $? -ne 0 ]; then
    t5="failed at 'test_cn' in mode 'CN'"
else
    t5="succeed in test #5"
fi

#6
sigflow extract -i tcga_laml.maf.gz -o test_results/test_maf_manual -m SBS -r 10 -T 4 --manual -g hg19
if [ $? -ne 0 ]; then
    t6_1="failed at 'test_manual' step 1 in mode 'SBS'"
else
    t6_1="succeed in test #6 step 1"
fi

sigflow extract -i tcga_laml.maf.gz -o test_results/test_maf_manual -m SBS -r 10 -T 4 --manual -N 3
if [ $? -ne 0 ]; then
    t6_2="failed at 'test_manual' step 2 in mode 'SBS'"
else
    t6_2="succeed in test #6 step 2"
fi

## Fitting PART:

#7
sigflow fit -i tcga_laml.maf.gz -o test_results/test_fitting -m MAF
if [ $? -ne 0 ]; then
    t7="failed at 'test_fitting' in mode 'MAF'"
else
    t7="succeed in test #7"
fi

#8
sigflow bt -i tcga_laml.maf.gz -o test_results/test_bt -m MAF -r 10
if [ $? -ne 0 ]; then
    t8="failed at 'test_bt' in mode 'MAF'"
else
    t8="succeed in test #8"
fi

## SigProfiler PART
sigflow extract -i tcga_laml.maf.gz -o test_results/test_sigprofiler -m MAF -r 10 -T 4 --sigprofiler
if [ $? -ne 0 ]; then
    t9="failed at 'test_sigprofiler' in mode 'MAF'"
else
    t9="succeed in test #9"
fi

## Output test results
echo ===========================================
echo "             Test Results               "
echo $t1
echo $t2
echo $t3
echo $t4
echo $t5
echo "$t6_1; $t6_2"
echo $t7
echo $t8
echo $t9
echo ===========================================
