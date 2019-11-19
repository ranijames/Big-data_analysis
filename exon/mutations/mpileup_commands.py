__author__ = 'ARJ'


MPILEUP_DELETIONS_CMD = """
    /home/mpschr/bin/samtools/bin/samtools mpileup -A -m 3 -p -B -d 100000 -Q 0 -t DP -t AD -t ADF -t ADR
    -f {GENOME}
    --positions {BED}
    {TUMOR_SAMPLES}
    {NORMAL_SAMPLE} >
    {OUTPUT_MPILEUP}
"""


#    /home/mpschr/bin/samtools/bin/samtools mpileup -u -A -m 3 -p -B -d 100000 -Q 0 -t DP -t AD -t ADF -t ADR

MPILEUP_CMD = """/home/mpschr/bin/samtools/bin/samtools mpileup -L 1000000 -u -A -B -d 100000 -Q 0 -t DP -t AD -t ADF -t ADR --VCF \
    -f {GENOME} \
    --positions {BED} \
    {TUMOR_SAMPLES} \
    {NORMAL_SAMPLE} | \
    /home/mpschr/bin/samtools/bin/bcftools norm -m-any -f /home/mpschr/bin/bcbionextgen/data/genomes/Hsapiens/GRCh37/seq/GRCh37.fa -  | \
     grep -v '<\*>' | \
    /home/mpschr/bin/samtools/bin/bcftools call -c | /home/mpschr/bin/samtools/bin/bcftools view -m 2 > {OUTPUT_VCF}
"""
