#include <unistd.h>
#include <stdio.h>
#include "vcfpaste.h"


int usage()
{
    fprintf(stderr, "\n");
    fprintf(stderr, "About : Merge VCF GT PL.\n");
    fprintf(stderr, "Version : 0.0.1, build with htslib version : %s\n", hts_version());
    fprintf(stderr, "Usage : vcfpaste -s <site_vcf> -1 <prior_vcf> -2 <behind_vcf> -o <out_vcf>\n");
    fprintf(stderr, "   -s <file>            sites vcf (refpanel sites) [request]\n");
    fprintf(stderr, "   -1 <file>            prior_vcf (snps vcf.gz), 将该文件中genotype信息合并至site_vcf， 若和behind_vcf有相同位点，优先使用该文件 [request]\n");
    fprintf(stderr, "   -2 <file>            behind_vcf (indels vcf.gz), 将该文件中genotype信息合并至site_vcf [request]\n");
    fprintf(stderr, "   -o <file>            write output to a file [request]\n");
    fprintf(stderr, "   -h                   print help info\n");
    fprintf(stderr, "\n");
    return 1;
}

int main(int argc, char * argv[]) {
    int ch;
    char *site_vcf;
    char *prior_vcf;
    char *second_vcf;
    char *out_vcf;


    while ((ch = getopt(argc, argv, "s:1:2:o:h")) != -1) {
        switch (ch) {
            case 's':
                site_vcf = optarg;
                break;
            case '1':
                prior_vcf = optarg;
                break;
            case '2':
                second_vcf = optarg;
                break;
            case 'o':
                out_vcf = optarg;
                break;
            case 'h':
                return usage();
            case '?':
                printf("Unknown option: %c\n", (char) optopt);
                break;
        }
    }
    if (optind < 2){
        return usage();
    }
    scan_vcf(site_vcf, prior_vcf, second_vcf, out_vcf);
    return 0;
}