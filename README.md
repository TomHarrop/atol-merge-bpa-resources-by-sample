# atol-merge-bpa-resources-by-sample

1. Query the DP to get a list of all packages
2. Function:
   1. For each package, get the `sample_id` and/or `bioplatforms_sample_id`
   2. Write output to a queue (see
      https://github.com/TomHarrop/tc_alignment_params/blob/main/src/process_trimal_files.py)
3. Function
   1. Pull `sample_id`s from the queue
   2. Query the API by sample ID to get a list of resources
   3. Write this somewhere
4. Have some logic to go through the output of 3 and match samples in a way
   that they can be used as input for assembly pipelines
