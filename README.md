# atol-merge-bpa-resources-by-sample

## Aim

Produce `yaml` files from the BPA data portal that can be passed directly to
`sanger-tol/genomeassembly`.

## Overview

1. Download the BPA data portal (`ckanapi` has a "dump" function, [used here](https://github.com/TomHarrop/atol-merge-bpa-resources-by-sample/blob/a4ee60bab236296bed3b08369cc2b525f8a9687c/workflow/Snakefile#L20-L35))
2. Process the data into dictionaries of `{sample_id: resource}`. There is more than one `sample_id`, see [the script](https://github.com/TomHarrop/atol-merge-bpa-resources-by-sample/blob/main/workflow/scripts/parse_bpa_json.py) 
3. Set up logic to go through the output of 2 and match samples in a way
   that they can be used as input for assembly pipelines
4. Generate a yaml file for each `query_id`


## TODO

- [ ] define `query_id` (it's the identifier that groups datasets into an
  assembly, but we don't know what it is yet)
- [ ] package
- [ ] container / conda recipe
- [ ] automate (e.g. systemd unit)

## Notes

This was an initial draft

1. Query the DP to get a list of all packages
2. Function:
   1. For each package, get the `sample_id` and/or `bioplatforms_sample_id`
   2. Write output to a queue (see
      https://github.com/TomHarrop/tc_alignment_params/blob/main/src/process_trimal_files.py)
3. Function
   1. Pull `sample_id`s from the queue
   2. Query the API by sample ID to get a list of resources
   3. Write this somewhere

