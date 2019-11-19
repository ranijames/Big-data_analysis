import argparse

__author__ = 'ARJ'

CMD_ZCAT_TEMPLATE = 'zcat {}'
CMD_VCF_FILTER = 'vcffilter -f "! DB & MQ > 30 & QD > 20"'
CMD_GREP_CODING = 'grep -i "\|coding\|ENST"'
CMD_MUT_TYPE_SUMMARY = 'grep -oP "TYPE=\w+" | sort | uniq -c'

CMD_GEMINI = 'gemini query --header -q "select chrom,start,end,gene,ref,alt,aa_change,rs_ids,cosmic_ids,aaf_adj_exac_all,type,qual,gt_depths,qual_depth,allele_count,aaf from variants where not in_dbsnp and not in_exac or aaf_adj_exac_all < 0.000001 order by type" {}'
CMD_GEMINI2 = 'gemini query --header -q "select chrom,start,end,gene,ref,alt,aa_change,rs_ids,cosmic_ids,aaf_adj_exac_all,type,qual,gt_depths,qual_depth,allele_count,aaf from variants where (is_exonic or is_coding) and qual > 10 and (not in_dbsnp and not in_exac or aaf_adj_exac_all < 0.000001) order by qual DESC" {}'


def pipe_commands(commands):
    assert type(commands) is list
    return ' | '.join(commands)

class VCFFilteringTask:

    def run(self, input_files, output_file):

        for ifile in input_files:
            readfile = CMD_ZCAT_TEMPLATE.format(ifile)
            cmds = [readfile, CMD_VCF_FILTER, CMD_GREP_CODING, CMD_MUT_TYPE_SUMMARY]
            print(pipe_commands(cmds))



def cmdline():
    # Parse arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('-i', action='append', dest='input_files', default=[])
    parser.add_argument('-o', dest='output_files', action="append")
    options = parser.parse_args()

    # Run the command line
    task = VCFFilteringTask()
    task.run(
        options.input_files,
        options.output_files
    )

if __name__ == "__main__": # detects if called from system
    print('Command line')
    cmdline()
else:
    print('No command line')


