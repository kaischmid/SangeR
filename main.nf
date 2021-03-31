#!/usr/bin/env nextflow

ab1 = Channel.watchPath( './data/*.ab1', 'create,modify' )

process histgram {

  container 'kaischmid/sange_r'
  tag "${ab1}"
  echo true

  input:
  file ab1 from ab1

  output:


  """
  Rscript $workflow.projectDir/scripts/Sanger.R $ab1

  """

}
