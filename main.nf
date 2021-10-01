#!/usr/bin/env nextflow

ab1 = Channel.watchPath( './data/*.ab1', 'create,modify' )

//script load
sanger_script = file('/home/dir/scripts/R/sanger.R')


process histogram {

  publishDir "/mnt/miracum_archiv/SangerDone", mode:'copy'
  container 'kaischmid/sange_r'
  tag "${ab1}"
  echo true

  input:
  set basename, file(ab1) from ab1
  file(sanger_script)

  output:

  set basename, file('*mosaic.png'), file('*mutated.png'), file('*.csv') into png

  """
  ./sanger_script $ab1

  """

}


process move{

   input:
   set basename, file (mosaic), file (mutated), file (csv) from png

   output:
   set basename, file (mosaic), file (mutated), file (csv) into moved

   """
   mv /mnt/miracum_archiv/ImportEPIC/sanger/${basename}.ab1 /mnt/miracum_archiv/SangerDone/
   """

}

process submit{

   queueSize = 1

   input:
   file(submit)
   val basename, file (mosaic), file (mutated), file (csv) from moved

   shell:
   '''
   ./submit.pl !{basename}_mosaic.png !{basename}_mutated.png !{basename}.csv !{basename}.ab1
   '''

}
