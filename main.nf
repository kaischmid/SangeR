#!/usr/bin/env nextflow

params.input_dir = './'
params.output_dir = './'
params.POI = './'

//submit
aubmit = file('home/schmid/sanger/submit.pl')

ab1 = Channel.watchPath( "${params.input_dir}*.ab1", 'create,modify' )

//script load
run_nextflow = file('./run_nextflow.R')

process histogram {

  publishDir "${params.output_dir}/SangerResults", mode:'copy'
  tag "${ab1}"
  echo true


  input:
  set basename, file(ab1) from ab1
  file(run_nextflow)
  file(${params.POI})

  output:
  set basename, file('*.png') optional true into png
  file('*.csv') into png

  """
  Rscript run_nextflow.R $ab1 ${params.POI}
  """

}


 process move{
   input:
   set basename, file(mutated) from png
   file(csv) from csv

   output:
   val basename into moved

   """
   mv /mnt/miracum_archiv/ImportEPIC/sanger/${basename}.ab1 /mnt/miracum_archiv/SangerResults/
   """


 process submit{

   queueSize = 1
   file(submit)

   input:
   file(submit)
   val basename from moved

   shell:
   '''
   ./submit.pl !basename !{basename}.png !{basename}.csv
   '''
}
