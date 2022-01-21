#!/usr/bin/env nextflow

params.input_dir = './'

ab1 = Channel.watchPath( "${params.input_dir}*.ab1", 'create,modify' )

//script load
run_nextflow = file('./run_nextflow.R')

params.output_dir = './'

process histogram {

  publishDir "${params.output_dir}/SangerResults", mode:'copy'
  container 'kaischmid/sange_r'
  tag "${ab1}"
  echo true

  input:
  file(ab1) from ab1
  file(run_nextflow)

  output:
  file('*.png') optional true into png
  file('*.csv') into csv

  """
  Rscript run_nextflow.R $ab1

  """

}

/*
* process move{
*   input:
*   set basename, file (mosaic), file (mutated), file (csv) from png
*
*   output:
*   set basename, file (mosaic), file (mutated), file (csv) into moved
*
*   """
*   mv /mnt/miracum_archiv/ImportEPIC/sanger/${basename}.ab1 /mnt/miracum_archiv/SangerDone/
*   """
*
*
* process submit{
*
*   queueSize = 1
*
*   input:
*   file(submit)
*   val basename, file (mosaic), file (mutated), file (csv) from moved
*
*   shell:
*   '''
*   ./submit.pl !{basename}_mosaic.png !{basename}_mutated.png !{basename}.csv !{basename}.ab1
*   '''
}*/
