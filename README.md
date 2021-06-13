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
