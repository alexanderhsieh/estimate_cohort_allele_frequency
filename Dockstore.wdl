## Copyright Broad Institute, 2020
## This script takes as input an array of single-sample gvcfs and merges them into a single gvcf
## using bcftools merge, then calculates cohort allele frequency for each variant position based on 
## genotypes
##  
## TESTED: 
## Versions of other tools on this image at the time of testing:
##
## LICENSING : This script is released under the WDL source code license (BSD-3) (see LICENSE in https://github.com/broadinstitute/wdl). 
## Note however that the programs it calls may be subject to different licenses. Users are responsible for checking that they are authorized to run all programs before running this script. 
## Please see the docker for detailed licensing information pertaining to the included programs.
##

###########################################################################
#WORKFLOW DEFINITION
###########################################################################
workflow estimate_cohort_allele_frequency {

  Array[File] gvcfs
  Array[File] index_files
  String outprefix

  File estimation_script

  parameter_meta {
    gvcfs: "array of files representing single-sample gvcfs"
    index_files: "array of corresponding index files"
    outprefix: "string; output file prefix (to be appended to .g.vcf.gz)"
    estimation_script: "estimate_cohort_AF.py script"
  }
  meta{
    author: "Alex Hsieh"
    email: "ahsieh@broadinstitute.org"
  }

  call bcftools_merge {
    input:
    gvcfs = gvcfs,
    index_files = index_files,
    outprefix = outprefix

  }

  call split_gvcf {
    input:
    gvcf = bcftools_merge.out
  }

  scatter (idx in range(length(split_gvcf.out))) {

    call estimate_cohort_AF {
      input:
      script = estimation_script,
      gvcf = split_gvcf.out[idx],
      shard = "${idx}"
    }

  }

  call gather_shards {
    input:
    shards = estimate_cohort_AF.out,
    outprefix = outprefix
  }


  #Outputs 
  output {
      File cohort_AC_file = gather_shards.out
  }

}





###########################################################################
#Task Definitions
###########################################################################

# merges array of single-sample gvcfs into a single cohort gvcf
task bcftools_merge {
  Array[File] gvcfs
  Array[File] index_files
  String outprefix
  String outfname = "${outprefix}.g.vcf"

  command {
    bcftools merge -l ${write_lines(gvcfs)} -o ${outfname} -O v
  }

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File out = "${outfname}"
  }
}


## splits vcf by chromosome
task split_gvcf {

  File gvcf # input gvcf
  String outprefix = basename(gvcf, 'g.vcf.gz')

  command <<<
    # pull header lines
    zgrep "^#" ${gvcf} > header.txt

    # sort input vcf and bgzip
    sort -k1,1V -k2,2n ${gvcf} | bgzip -c > "${gvcf}.gz"

    # tabix index input vcf
    tabix -p vcf "${gvcf}.gz"

    # split vcf by chromosome - use tabix -l to get all contig names from tabix index
    for i in $(tabix -l ${gvcf}.gz)
    do 
      (cat header.txt; tabix ${gvcf}.gz $i) > "${outprefix}.$i.g.vcf"
    done


  >>>

  output {
    Array[File] out = glob("*.vcf") 
  }

  runtime {
    docker: "gatksv/sv-base-mini:cbb1fc"
  }

}

# estimates cohort AF
task estimate_cohort_AF{
  
  File script

  File gvcf
  String shard

  String outprefix = basename(gvcf, '.g.vcf')
  String outfname = "AC.${outprefix}.${shard}.txt"

  command <<<
    python ${script} -i ${gvcf} -o ${outfname}
  >>>

  runtime {
    docker: "mwalker174/sv-pipeline:mw-00c-stitch-65060a1"
  }

  output {
    File out = "${outfname}"
  }
}

#Gathers shards of cohort AC files into single file
task gather_shards {

  Array[File] shards 
  String outprefix

  command <<<

    while read file; do
      cat $file >> "tmp.cat.txt"
    done < ${write_lines(shards)};

    grep "^var_id" "tmp.cat.txt" > "header.txt"

    (cat header.txt; grep -v "^var_id" "tmp.cat.txt") > AC.${outprefix}.txt"

  >>>

  output {
    File out = "AC.${outprefix}.txt"
  }
}
