#!/usr/bin/env nextflow

ab1 = Channel.watchPath( './data/*.ab1', 'create,modify' )

//script load
sanger_script = file('/home/dir/scripts/R/sanger.R')


process histgram {

  container 'kaischmid/sange_r'
  tag "${ab1}"
  echo true

  input:
  file ab1 from ab1
  file(sanger_script)

  output:

  """
  ./sanger_script $ab1

  """

}
