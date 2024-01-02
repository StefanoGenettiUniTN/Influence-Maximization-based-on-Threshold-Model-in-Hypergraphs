# Influence Maximization based on Threshold Model in Hypergraphs


## HCI-TM-algorithm
The HCI-TM algorithm (Hypergraph-based Threshold Model) is a threshold propagation algorithm based on the hypergraph model, designed to identify the node set with maximum influence. And the open source code includes hypergraph building module, threshold propagation module, HCI1-TM algorithm module and HCI2-TM algorithm module.

For more details, please refer the following paper:
> Renquan Zhang, Xilong Qu, Qiang Zhang, Xirong Xu, Sen Pei: Influence Maximization based on Threshold Model in Hypergraphs.
> arXiv:2306.13458


## About File

Files are presented in there folder, where the folder 'code' contains file 'OpenCode.cpp', including hypergraph building module, threshold propagation module, HCI1-TM algorithm module and HCI2-TM algorithm module, the folder 'data' contains 8 real-world networks using in this paper and the folder 'test' is the test file for the project. In folder 'data', the hypergraph is stored in the txt file. Each line in the file represents a hyperedge, the fisrst number represent the hyperedge id, other positive numbers represent the node id, and -1 represents the end of the hyperedge. In folder 'test', the file 'hypergraph-2-ER-10000-3000.txt' stored a ER hypergraph with 8751 nodes and 2993 hyperedges, the file 'Record-Detail-Tr0.5.txt' stored the proportion of seed nodes and the activation proportion of nodes in the activation process of HCI1-TM algorithm and HCI2-TM algorithm when threshold is set to 0.5, the file 'Record-HCI-Tr0.5.txt' stored the HCI1 and HCI2 of the seed node, the file 'Record-Tr0.5.txt' stored the number of seed, the runing time and the proportion of seed of HCI1-TM and HCI2-TM algorithm, the file 'Record-HCI1TM-SeedSet-Tr0.5.txt' and 'Record-HCI2TM-SeedSet-Tr0.5.txt' stored seed set selected by HCI1-TM algorithm and HCI2-TM algorithm.


## Contact


If you have any question about the paper or the code, please contact us. **Xilong Qu**, **[quxilong@mail.dlut.edu.cn]**

Please cite that paper if you use this code. Thanks!
