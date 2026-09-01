07-DML-adult-comparison
================
Kathleen Durkin
2025-05-16

In `06.1-differential-methylation-split-lifestage` I identified loci in
the larval samples that are differentially methylated based on parental
OA treatment. I now want to see whether any of those DMLs were also
present in the adults

Note: I’m first checking the adult DML tables, but I may also want to
look a the methylomes of both the larvae and the larvae’s parents in
more detail, since not all adult oysters whose gametes were sequenced in
Venkatamaran et al. 2024 were used for fertilization crosses.

Download adult DMLs (both female (eggs) and male (sperm))

``` bash
curl -L https://github.com/sr320/ceabigr/raw/refs/heads/main/data/female_dml.bed -o ../data/adult_female_dml.bed
curl -L https://github.com/sr320/ceabigr/raw/refs/heads/main/data/male_dml.bed -o ../data/adult_male_dml.bed
```

    ##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
    ##                                  Dload  Upload   Total   Spent    Left  Speed
    ##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
    ## 100  5018  100  5018    0     0  10699      0 --:--:-- --:--:-- --:--:-- 10699
    ##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
    ##                                  Dload  Upload   Total   Spent    Left  Speed
    ##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0
    ## 100  159k  100  159k    0     0   341k      0 --:--:-- --:--:-- --:--:--  341k

``` bash
echo "Larval DMLs that are also differentially methylated in female parents"
/home/shared/bedtools2/bin/bedtools intersect -a ../output/06.1-differential-methylation-split-lifestage/larvae_diff25_DML.bed -b ../data/adult_female_dml.bed

echo ""

echo "Larval DMLs that are also differentially methylated in male parents"
/home/shared/bedtools2/bin/bedtools intersect -a ../output/06.1-differential-methylation-split-lifestage/larvae_diff25_DML.bed -b ../data/adult_male_dml.bed
```

    ## Larval DMLs that are also differentially methylated in female parents
    ## 
    ## Larval DMLs that are also differentially methylated in male parents
    ## NC_035784.1  55195530    55195531    *   48.7264348554671
    ## NC_035787.1  28915240    28915241    *   -36.2175362175362
