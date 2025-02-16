# HPC-Project

Various versions of the program are available:
* `versions/main_local.c`: fishes advance locally to their process and communicate sysnthetic
  information with other processes only periodically. The information shared are
  the baricenter of local fishes, their average weight and their total weight
  improvement. All these information and then aggregated to that of other
  processes to compute a new global baricenter which is used as reference for
  an additional `collective_volitive_move`
* `versions/main_b.c`: 1 AllGather before collective_istinctive_move
* `versions/main_a.c`: 3 AllGather before collective_istinctive_move
* `versions/main_old.c`: 1 AllGather before collective_istinctive_move, but fishes
  have a less efficient representation
* `versions/main_allreduce.c`: Multiple Allreduce to compute global sums and maximum values 