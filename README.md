Clustering of Long Noisy Reads

Long reads can be quite noisy, and quality values are used to detect wrong bases.
Q-cluster is a tool to cluster short reads based on k-mers and quality values, here the code:
https://github.com/CominLab/QCluster
The idea is to extend this tool to be able to handle larger values of k,
and test its effectiveness as in this paper.


TEST
The command for the tests are the following.

TEST A:    QCluster -d d -c 5 -k 5 -t 5 -S 0 -w Input/sequences.fastq

TEST B:    QCluster -d d -c 5 -k 5 -t 5 -S 0 -w Input/simulated.fastq

TEST C:    QCluster -d d -c 5 -k 20 -t 5 -S 0 -w Input/sequences.fastq


FUTURE APPROACH
At the current state of the project, the software does not work properly. We have observed the following wrong behaviours:
1) the final assignment of clusters is incomplete. Only a small part of input sequences is assigned to a cluster (around 57 sequences) while the remaining sequences are assigned to 0;
2) "out_of_range" exception when d2ast method is used ("a" option related to distance computing methods)
3) In different cases, we have encountered a "segmentation fault (core dumped)" caused by a wrong deallocation of memory.

We have tried to formilize some possibile solution to the previous errors:
1) it may be caused from a bad computing of the centroids and distances through appropriate functions (centroid.cc and dists.cc). In particular, we have some doubts about our use of expected_freq and expected_qual since we have splitted them into expected_freq_1 nad expected_qual_1 in order to distinguish the results between freq and freq_1. Thus, we think the right approach is to study the behaviours of managed data for each data structure in order to check the correct computing of these elements.
2) it should be an error caused by the main "for" loop (dists.cc file, line 109), where the copy of the quality structure (x variable) is recomputed. The solution might be the checking of map access.
3) it is caused by some errors during the deallocation of data structures. In this case, it should be necessary only to check if all data structures are deallocated correctly.

Additional informations:
- mm_centroid function (centroid.cc file, 82 line) is not completed, since we were not able to interpret the algorithm in the right way. It shouldn't cause the previous second error since it does not throw the exception using the distance computing method indicated with "c" option.
- it may be a good idea try some test using ordered map structures
