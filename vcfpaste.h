//
// Created by 黄志博 on 2020/12/31.
//

#ifndef VCFPASTE_VCFPASTE_H
#define VCFPASTE_VCFPASTE_H


#include <htslib/hts.h>
#include <htslib/sam.h>
#include <htslib/synced_bcf_reader.h>
#include <htslib/vcf.h>


void scan_vcf(char* fmain, char* f1, char* f2, char* out);

#endif //VCFPASTE_VCFPASTE_H
