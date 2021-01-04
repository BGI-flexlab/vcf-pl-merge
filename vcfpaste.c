//
// Created by 黄志博 on 2020/12/31.
//
#include "vcfpaste.h"

void scan_vcf(char *fmain, char *f1, char *f2, char *out) {
    bcf_srs_t * sr =  bcf_sr_init();
//    sr->collapse = COLLAPSE_NONE;
    sr->require_index = 1;
    if(!(bcf_sr_add_reader (sr, fmain))) fprintf(stderr,"Problem opening index file for %s]\n", fmain);
    if(!(bcf_sr_add_reader (sr, f1))) fprintf(stderr, "Problem opening index file for %s]\n", f1);
    if(!(bcf_sr_add_reader (sr, f2))) fprintf(stderr, "Problem opening index file for %s]\n", f2);

    htsFile * fp = hts_open(out, "wz");
    bcf_hdr_t * hdr = bcf_sr_get_header(sr, 1);
    bcf_hdr_remove(hdr, BCF_HL_INFO, NULL);
    bcf_hdr_write(fp, hdr);
    bcf1_t *rec = bcf_init1();

    int * genotypes = (int*)malloc(2*sizeof(int));
    int n_variants = 0;
    int n_variants_prior = 0;
    int n_variants_behind = 0;
    int nset;
    int ngl_main, ngl_arr_main = 0, *gl_arr_main = NULL;
    int ngt, *gt_arr = NULL, ngt_arr = 0;
    bcf1_t * line_main, * line_gt, * line_gt2;
    while ((nset = bcf_sr_next_line(sr))) {
        bcf_clear1(rec);
        if (!bcf_sr_has_line(sr, 0)) continue;
        n_variants ++;

        line_main =  bcf_sr_get_line(sr, 0);

        if (bcf_sr_has_line(sr, 1)){
            line_gt =  bcf_sr_get_line(sr, 1);
            ngt = bcf_get_genotypes(sr->readers[1].header, line_gt, &gt_arr, &ngt_arr);
            ngl_main = bcf_get_format_int32(sr->readers[1].header, line_gt, "PL", &gl_arr_main, &ngl_arr_main);
            n_variants_prior ++;
        } else if (bcf_sr_has_line(sr, 2)){
            line_gt =  bcf_sr_get_line(sr, 2);
            ngt = bcf_get_genotypes(sr->readers[2].header, line_gt, &gt_arr, &ngt_arr);
            ngl_main = bcf_get_format_int32(sr->readers[2].header, line_gt, "PL", &gl_arr_main, &ngl_arr_main);
            n_variants_behind ++;
        } else {
            ngt = -1;
            ngl_main = bcf_get_format_int32(sr->readers[2].header, line_gt, "PL", &gl_arr_main, &ngl_arr_main);
        }

        rec->rid = line_gt->rid;
        rec->pos = line_gt->pos;
        bcf_update_alleles(hdr, rec, line_main->d.allele, line_main->n_allele);

        if(ngt == 2){
            bcf_update_genotypes(hdr, rec, gt_arr, 2);
            bcf_update_format_int32(hdr, rec, "PL", gl_arr_main, ngl_main);
        } else{
            genotypes[0] = bcf_gt_missing;
            genotypes[1] = bcf_gt_missing;
            bcf_update_genotypes(hdr, rec, genotypes, 2);
//            bcf_update_format_int32(hdr, rec, "PL", 0, 0);
        }
        bcf_write1(fp, hdr, rec);
    }
    bcf_sr_destroy(sr);
    free(genotypes);
    bcf_destroy1(rec);
//    bcf_hdr_destroy(hdr);
    if (hts_close(fp)) fprintf(stderr, "Non zero status when closing VCF/BCF file descriptor");
    printf("n_variants       : %d\n", n_variants);
    printf("n_variants_prior : %d\n", n_variants_prior);
    printf("n_variants_behind: %d\n", n_variants_behind);
}



