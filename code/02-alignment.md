02-alignment
================
Kathleen Durkin
2025-04-29

I’ll be performing alignment of WGBS reads to the *C. virginica* genome
using `bismark`. Sam has already generated scripts for this purpose as
part of the [`ceasmallr`](https://github.com/sr320/ceasmallr/tree/main_)
project, but I’ll be rerunning alignment for the `FISH 546`
bioinformatics class (SP 2025).

Download trimmed reads

``` bash
# Run wget to retrieve FastQs and MD5 files
# Note: the --no-clobber command will skip re-downloading any files that are already present in the output directory
wget \
--directory-prefix ../data/trimmed-WGBS-reads \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 4 \
--no-host-directories \
--no-parent \
--quiet \
--no-clobber \
--accept="*.fq.gz,checksums.md5,multiqc_report.html,fastq_pairs.txt" https://gannet.fish.washington.edu/gitrepos/ceasmallr/output/00.00-trimming-fastp/
```

We should have 64 trimmed reads (32 samples, R1 and R2 for each sample)

``` bash
echo "How many trimmed reads files were downloaded?"
ls ../data/trimmed-WGBS-reads/*.fq.gz | wc -l
```

    ## How many trimmed reads files were downloaded?
    ## 64

Download C.virginica genome

``` bash
wget \
--directory-prefix ../data/ \
--recursive \
--no-check-certificate \
--continue \
--cut-dirs 3 \
--no-host-directories \
--no-parent \
--quiet \
--no-clobber \
 https://gannet.fish.washington.edu/gitrepos/ceasmallr/data/Cvirginica_v300/
```

Download scripts files

``` bash
curl -L -O https://raw.githubusercontent.com/sr320/ceasmallr/refs/heads/main/code/02.01-bismark-bowtie2-alignment-SLURM-array.sh

curl -L -O https://raw.githubusercontent.com/sr320/ceasmallr/refs/heads/main/code/02.02-bismark-SLURM-job.sh
```

    ##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
    ##                                  Dload  Upload   Total   Spent    Left  Speed
    ##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0  0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0100  3220  100  3220    0     0  12365      0 --:--:-- --:--:-- --:--:-- 12337
    ##   % Total    % Received % Xferd  Average Speed   Time    Time     Time  Current
    ##                                  Dload  Upload   Total   Spent    Left  Speed
    ##   0     0    0     0    0     0      0      0 --:--:-- --:--:-- --:--:--     0100  1258  100  1258    0     0   5299      0 --:--:-- --:--:-- --:--:--  5308

Modify scripts as needed.

Changes I made:

.sh file:

- changed repo directory and output directory

.job file:

- changed the notification email and working directory from Sam’s info
  to mine

- changed final line, which provides file path to the .sh script

Ran job script by running the below line from Klone terminal:

    sbatch --array=0-$(($$(wc -l < ../data/trimmed-WGBS-reads/fastq_pairs.txt) - 1)) 02.02-bismark-SLURM-job.sh

Finished running after ~26 hours, check some of the outputs. I’m a
little concerned because the email I recieved from Hyak notifying me of
completion said
`Slurm Array Summary Job_id=25714300_* (25714300) Name=bismark_job_array Failed, Mixed, ExitCode [0-2]`

``` bash
cd ../output/02.01-bismark-bowtie2-alignment-SLURM-array/

echo "Number of error files (should have 32):"
ls *.err | wc -l

echo "Number of .out files (should have 32):"
ls *.out | wc -l

echo "Number of bismark summary files (should have 32):"
ls *bismark_summary.txt | wc -l

find . -type f -name "*.out" -size -100c
```

    ## Number of error files (should have 32):
    ## 33
    ## Number of .out files (should have 32):
    ## 33
    ## Number of bismark summary files (should have 32):
    ## 28
    ## ./bismark_job_25714300_0.out
    ## ./bismark_job_25714300_12.out
    ## ./bismark_job_25714300_13.out
    ## ./bismark_job_25714300_9.out
    ## ./bismark_job_25714300_10.out

There are 33 error and .out files because I mistakenly removed a +1 in
the job command, so it started at the wrong index. That means the
jobid_0 files are basically duds. So we actually have 32 of each, as
expected.

There are only 28 bismark summary files, instead of the expected 32.
Excluding the dud jobid_0.out, there are 4 .out files that are very
small (\<100 B): ./bismark_job_25714300_12.out
./bismark_job_25714300_13.out ./bismark_job_25714300_9.out
./bismark_job_25714300_10.out

Each of these files just contains the following:

    Contents of pair:


    Files  and  have already been processed. Exiting.

So… for each of these sub jobs, the read pair is interpreted as having
already been processed. Let’s look at all the -bismark_summary.txt files
(I think there should only be one of these for each pair of processed
reads)

These files seem to be named in the format:

SampleID-jobSubID-bismark_summary.txt

``` bash
cd ../output/02.01-bismark-bowtie2-alignment-SLURM-array/
  
ls *bismark_summary.txt | sort
```

    ## CF01-CM02-Larvae-8-bismark_summary.txt
    ## CF02-CM02-Zygote-3-bismark_summary.txt
    ## CF03-CM03-Zygote-1-bismark_summary.txt
    ## CF03-CM03-Zygote-5-bismark_summary.txt
    ## CF03-CM04-Larvae-14-bismark_summary.txt
    ## CF03-CM04-Larvae-15-bismark_summary.txt
    ## CF03-CM04-Larvae-4-bismark_summary.txt
    ## CF03-CM05-Larvae-16-bismark_summary.txt
    ## CF03-CM05-Larvae-6-bismark_summary.txt
    ## CF06-CM01-Zygote-19-bismark_summary.txt
    ## CF06-CM01-Zygote-20-bismark_summary.txt
    ## CF06-CM02-Larvae-11-bismark_summary.txt
    ## CF06-CM02-Larvae-22-bismark_summary.txt
    ## CF07-CM02-Zygote-24-bismark_summary.txt
    ## CF07-CM02-Zygote-28-bismark_summary.txt
    ## CF08-CM03-Zygote-23-bismark_summary.txt
    ## CF08-CM03-Zygote-26-bismark_summary.txt
    ## CF08-CM03-Zygote-2-bismark_summary.txt
    ## CF08-CM04-Larvae-7-bismark_summary.txt
    ## EF01-EM01-Zygote-31-bismark_summary.txt
    ## EF01-EM01-Zygote-32-bismark_summary.txt
    ## EF02-EM02-Zygote-17-bismark_summary.txt
    ## EF03-EM03-Zygote-18-bismark_summary.txt
    ## EF04-EM04-Zygote-21-bismark_summary.txt
    ## EF05-EM06-Larvae-25-bismark_summary.txt
    ## EF06-EM02-Larvae-27-bismark_summary.txt
    ## EF07-EM01-Zygote-29-bismark_summary.txt
    ## EF07-EM03-Larvae-30-bismark_summary.txt

Yikes. Ok it looks like, not only are many reads not run appropriately
through Bismark, many of the reads were processed twice
(e.g. `CF03-CM03-Zygote-1-bismark_summary.txt`,
`CF03-CM03-Zygote-5-bismark_summary.txt`). That means only ~half of the
reads were processed, not 28/32.

It also looks like

``` bash
cd ../output/02.01-bismark-bowtie2-alignment-SLURM-array/

echo "How many output .bam files?"
ls *.bam | wc -l

echo "How large are the .bam files?"
ls -lh *.bam
```

    ## How many output .bam files?
    ## 18
    ## How large are the .bam files?
    ## -rw-r--r-- 1 kdurkin1 nogroup 9.0G May  3 00:52 CF01-CM02-Larvae_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup 4.7G May  3 03:43 CF02-CM02-Zygote_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup 5.1G May  3 00:04 CF03-CM03-Zygote_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup    0 May  3 07:56 CF03-CM04-Larvae_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup    0 May  3 04:56 CF03-CM05-Larvae_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup    0 May  3 03:38 CF06-CM01-Zygote_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup    0 May  2 23:14 CF06-CM02-Larvae_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup    0 May  3 04:38 CF07-CM02-Zygote_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup    0 May  3 00:49 CF08-CM03-Zygote_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup 105M May  2 20:05 CF08-CM04-Larvae_R1_001.fastp-trim.REPAIRED_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup    0 May  3 01:30 EF01-EM01-Zygote_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup 5.8G May  2 23:38 EF02-EM02-Zygote_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup 3.8G May  2 23:38 EF03-EM03-Zygote_R1_001.fastp-trim.REPAIRED_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup 2.1G May  3 00:04 EF04-EM04-Zygote_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup  16G May  3 06:32 EF05-EM06-Larvae_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup 8.1G May  3 00:16 EF06-EM02-Larvae_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup 4.0G May  3 03:52 EF07-EM01-Zygote_R1_001.fastp-trim_bismark_bt2_pe.bam
    ## -rw-r--r-- 1 kdurkin1 nogroup 8.5G May  3 04:22 EF07-EM03-Larvae_R1_001.fastp-trim_bismark_bt2_pe.bam

There are only 18 output bam files, and 6 of those are empty (size 0)
