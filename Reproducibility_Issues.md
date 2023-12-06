# Reproducibility Notes For Monocle3 Analysis

## Summary

There are a couple of functions within the Monocle3 package that employ the use of a random number generator which can cause variation in the analysis each time a script is run. In this document those commands are identified and work-arounds are given to ensure reproducibility and that outputs from the commands is consistent between runs.

This document will be updated as necessary if patches are introduced to the Monocle3 package. For more information, and to check if there are new patches, please visit the [monocle3 github page](https://cole-trapnell-lab.github.io/monocle3/).


### reduce\_dimensions()

The following options need to be added to the reduce\_dimensions() command:
- umap.fast\_sgd = FALSE
- cores=1
- n\_sgd\_threads=1

Here is an example:

``reduce_dimension(cds,umap.fast_sgd = FALSE,cores=1,n_sgd_threads=1)``


### cluster\_cells()

The following commands/options need to be added when running cluster\_cells():
- the command set.seed('integer') needs to be immediately prior to cluster\_cells()
- within the cluster\_cells() the option random\_seed = 'integer' needs to be included

Here is an example:

``set.seed(17)``

``cds = cluster_cells(cds, random_seed = 17)``


### find\_gene\_modules()

The following commands/options need to be added when running cluster\_cells():
- the command set.seed('integer') needs to be immediately prior to cluster\_cells()
- within the cluster\_cells() the option random\_seed = 'integer' needs to be included

Here is an example:

``set.seed(17)``

``gene_module_df <- find_gene_modules(cds[pseu_pr_deg_ids,], random_seed = 17)``


### Additional Notes

The initial basis for these work-arounds is described in [this issue page](https://github.com/cole-trapnell-lab/monocle3/issues/277) on the Monocle3 Github.

For the reduce\_dimensions() command, setting the cores and sgd\_threads to 1 can increase the computational time it takes to run this command on very large datasets.

The set.seed command needs to be immediately prior to the given commands above. It cannot be listed at the start of a script once. Although it is not required that all integers be the same value throughout the script, it is recommended for simplicity.
