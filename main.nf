#!/usr/bin/env nextflow

params.input_dir = './'
params.output_dir = './'
params.POI = './'

//submit
submit = file('home/schmid/sanger/submit.pl')

ab1 = Channel.watchPath( "${params.input_dir}*.ab1", 'create,modify' )

//script load
run_nextflow = file('./run_nextflow.R')
POI = file("${params.POI}")

process histogram {

  publishDir "${params.output_dir}/SangerResults", mode:'copy'
  container 'kaischmid/sange_r'
  tag "${ab1}"
  echo true


  input:
  file(ab1) from ab1
  file(run_nextflow)
  file(POI)

  output:
  set basename, file('*.png') optional true into png
  file('*.csv') into csv

  """
  Rscript run_nextflow.R $ab1 ./POI
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
}

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
