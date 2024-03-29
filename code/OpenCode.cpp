#include <iostream>
#include <fstream>
#include <math.h>
#include <deque>
#include <string>
#include <cassert>
#include <vector>
#include <time.h> 
using namespace std;

#define Tr 0.8 //HYPEREDGE THRESHOLD

#define MAX_SEED_SET_SIZE 100

// DEFINES THE MAX NUMBER OF NODES IN THE HYPERGRAPH   
const int Max_Node_Size = 100000;
// DEFINES THE MAX NUMBER OF HYPEREDGES IN HYPERGRAPH  
const int Max_Edge_Size = 100000;

// Keeps track of the current number of activated nodes  
int AC_NODE_NUM = 0;
// Keeps track of the current number of activated hyperedges  
int AC_HYEDGE_NUM = 0;

// Records detailed information about the activation process. 
double RECORD[2][2500] = { 0 };

// Stores the HCI1 value for each seed node  
int HCI1_Rec[20000] = { 0 };
// Stores the HCI2 value for each seed node  
int HCI2_Rec[20000] = { 0 };

// Stores the seed set for HCI1-TM and HCI2-TM
int SeedSet[2][20000];
// Defines the precision or accuracy of the recorder  
const double Scale = 0.0004;
// Termination criteria parameter  
const double ActiveRadio = 0.9;

// Queue to store the activated nodes  
deque<int> ActiveNodeQueue;
// Queue to store the activated hyperedges  
deque<int> ActiveEdgeQueue;

// STRUCTURE TO HOLD INFORMATION ABOUT A NODE IN THE HYPERGRAPH  
struct NODE
{
	// Tells if the node is currently active or not  
	bool Active = false;
	// Degree of the node (number of hyperedges it is part of)  
	int Degree = 0;
	// IDs of the hyperedges this node belongs to  
	deque<int> Edge_id;
	// Marks if the node has been added to a queue before  
	bool Mark = false;
};

// STRUCTURE TO HOLD INFORMATION ABOUT A HYPEREDGE IN THE HYPERGRAPH  
struct HYPEREDGE
{
	// Tells if the hyperedge is currently active or not  
	bool Active = false;
	// Tells if the hyperedge is in a subcritical state or not  
	bool Subcritical = false;
	// Number of activated nodes in this hyperedge  
	int Activenum = 0;
	// Threshold value for the hyperedge activation  
	int Threshold = 0;
	// Hyperedge's cardinality  
	int Cardinality = 0;
	// IDs of the nodes associated with the hyperedge
	deque <int> Nodeinfo_id;
	// States of the nodes associated with this hyperedge (0 for activated, 1 for not activated)  
	deque <int> Nodeinfo_state;
};
// Definition of function 'constructhypergraph'  
int ConstructHypergraph(vector<NODE>& node, vector<HYPEREDGE>& hyperedge, string filename)
{
	// Open the specified file for reading  
	ifstream infile(filename);

	// Read the first integer from the file, representing the number of edges to be processed  
	int iter_num = 0;
	infile >> iter_num;

	// Initialization of variables to keep track of edge status  
	int edgeid = iter_num;  // initialization of edge id as the first integer  
	int iter_over = 0;  // initialization of over status as 0 (incomplete edge)  

	// Loop to process each integer in the file  
	while (infile >> iter_num)
	{
		// Check if the current integer is -1, indicating that the edge is completed  
		if (iter_num == -1)
		{
			iter_over = 1;  // setting over status to 1 (complete edge)  
		}
		// Check if the over status is 1, meaning a new edge should be started  
		else if (iter_over == 1)
		{
			edgeid = iter_num;  // updating edge id with new integer  
			iter_over = 0;  // resetting over status to 0 (new node-edge relationship being established)  
		}
		// Check if the over status is 0, indicating processing of node-edge relationship  
		else if (iter_over == 0)
		{
			// Updating edge id list for current node and increasing its degree  
			node[iter_num].Edge_id.push_back(edgeid);
			node[iter_num].Degree++;

			// Updating node id list for current edge and increasing its Cardinality (hyperdegree)  
			hyperedge[edgeid].Nodeinfo_id.push_back(iter_num);
			hyperedge[edgeid].Cardinality++;  // increasing hyperdegree of the edge  

			// Updating node state list for current edge and setting it to 0 (node not activated yet)  
			hyperedge[edgeid].Nodeinfo_state.push_back(0);  // node not activated yet  
		}
	}

	// Calculating threshold for each edge using hyperdegree and global variable tr  
	for (int i = 0; i < Max_Edge_Size; i++)
	{
		hyperedge[i].Threshold = int(ceil(hyperedge[i].Cardinality * Tr));  // calculating threshold and rounding up  

		// Checking if threshold of edge matches its activated node count plus 1, marking it subcritical  
		if (hyperedge[i].Threshold == hyperedge[i].Activenum + 1)
		{
			hyperedge[i].Subcritical = true;  // setting subcritical status to true  
		}
	}

	// Counting number of valid nodes (nodes with degree not zero)  
	int node_num = 0;  // initialization of valid node count as zero  
	for (int i = 0; i < Max_Node_Size; i++)
	{
		if (node[i].Degree != 0)  // checking if node degree is not zero  
		{
			node_num++;  // incrementing valid node count by 1  
		}
	}
	return node_num;  // returning number of valid nodes as function output  
}

void Activation_rule(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge, int seed[], int seedsize)//A function about activation rules
{
	// Initialize a queue to store the IDs of nodes that need to be activated
	deque<int> queue;

	// Add seed nodes to the queue
	for (int i = 0; i < seedsize; i++)
	{
		Node[seed[i]].Mark = true; // Mark the seed node as activated
		Node[seed[i]].Active = true; // Set the seed node to an active state
		queue.push_back(seed[i]); // Add the seed node ID to the queue
		ActiveNodeQueue.push_back(seed[i]); // Add the seed node ID to the active node queue
		AC_NODE_NUM++; // Increment the number of activated nodes by 1
	}
	// While there are still nodes in the queue  
	while (queue.size() > 0) 
	{
		// Get the ID of the node at the front of the queue  
		int nodeid = queue.front();
		// Remove the node from the front of the queue  
		queue.pop_front();

		// Iterate through all hyperedges related to the node with ID nodeid  
		for (int i = 0; i < Node[nodeid].Degree; i++) {
			// If the hyperedge related to nodeid is not active  
			if (Hyperedge[Node[nodeid].Edge_id[i]].Active == false) {
				// Iterate through all nodes in the hyperedge  
				for (int j = 0; j < Hyperedge[Node[nodeid].Edge_id[i]].Cardinality; j++) {
					// If the current node is nodeid  
					if (nodeid == Hyperedge[Node[nodeid].Edge_id[i]].Nodeinfo_id[j]) {
						// If the state of the node in the hyperedge is 0 (not active)  
						if (Hyperedge[Node[nodeid].Edge_id[i]].Nodeinfo_state[j] == 0) {
							// Increment the active count of the hyperedge  
							Hyperedge[Node[nodeid].Edge_id[i]].Activenum++;
							// Change the state of the node in the hyperedge to 1 (active)  
							Hyperedge[Node[nodeid].Edge_id[i]].Nodeinfo_state[j] = 1;

							// If the threshold of the hyperedge is reached  
							if (Hyperedge[Node[nodeid].Edge_id[i]].Threshold == Hyperedge[Node[nodeid].Edge_id[i]].Activenum + 1) {
								// Set the subcritical flag of the hyperedge to true  
								Hyperedge[Node[nodeid].Edge_id[i]].Subcritical = true;
							}
							else if (Hyperedge[Node[nodeid].Edge_id[i]].Activenum == Hyperedge[Node[nodeid].Edge_id[i]].Threshold) {
								// If the number of active nodes in the hyperedge reaches its threshold, activate the hyperedge  
								Hyperedge[Node[nodeid].Edge_id[i]].Active = true;
								Hyperedge[Node[nodeid].Edge_id[i]].Subcritical = false;
								AC_HYEDGE_NUM++; // Increment the number of activated hyperedges  
								ActiveEdgeQueue.push_back(Node[nodeid].Edge_id[i]); // Add the ID of the activated hyperedge to the active edge queue  

								// Iterate through all nodes in the hyperedge  
								for (int k = 0; k < Hyperedge[Node[nodeid].Edge_id[i]].Cardinality; k++) {
									// If a node in the hyperedge is not active  
									if (Hyperedge[Node[nodeid].Edge_id[i]].Nodeinfo_state[k] == 0) {
										// Activate the node and set its state to 1 (active)  
										Hyperedge[Node[nodeid].Edge_id[i]].Nodeinfo_state[k] = 1;
										Node[Hyperedge[Node[nodeid].Edge_id[i]].Nodeinfo_id[k]].Active = true; // Set the node to active  

										// If the node has never been added to the queue before, add it to the queue and mark it as active for further propagation  
										if (Node[Hyperedge[Node[nodeid].Edge_id[i]].Nodeinfo_id[k]].Mark == false) {
											AC_NODE_NUM++; // Increment the number of activated nodes  
											queue.push_back(Hyperedge[Node[nodeid].Edge_id[i]].Nodeinfo_id[k]); // Add the activated node to the queue  
											ActiveNodeQueue.push_back(Hyperedge[Node[nodeid].Edge_id[i]].Nodeinfo_id[k]); // Add the ID of the activated node to the active node queue  
											Node[Hyperedge[Node[nodeid].Edge_id[i]].Nodeinfo_id[k]].Mark = true;//Mark the activated node has been added to the queue
										}
									}
								}
							}
						}
					}
				}
			}
		}
	}
}


// The HGRecover function is used to recover the hypergraph for reuse.  
void HGRecover(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge)
{
	// Restore node states  
	for (int i = 0; i < Max_Node_Size; i++)
	{
		Node[i].Active = false;   // Set the node to inactive status  
		Node[i].Mark = false;     // Mark the node as unmarked  
	}

	// Restore hyperedge states  
	for (int i = 0; i < Max_Edge_Size; i++)
	{
		Hyperedge[i].Active = false;   // Set the hyperedge to inactive status  
		Hyperedge[i].Activenum = 0;     // Reset the number of active nodes in the hyperedge to 0  

		// Restore node states within the hyperedge  
		for (int j = 0; j < Hyperedge[i].Cardinality; j++)
		{
			Hyperedge[i].Nodeinfo_state[j] = 0;   // Reset the state of nodes in the hyperedge to 0  
		}

		Hyperedge[i].Subcritical = false;  // Mark the hyperedge as non-subcritical  
	}

	// Clear the active edge and node queues  
	ActiveEdgeQueue.clear();
	ActiveNodeQueue.clear();

	// Reset the counts of active edges and nodes to 0  
	AC_HYEDGE_NUM = 0;
	AC_NODE_NUM = 0;
}


// This function checks if the hyperedge with ID hyperedgeid becomes subcritical after removing nodes with IDs nodeid1 and nodeid2.  
int IsSubcritical(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge, int nodeid1, int nodeid, int hyperedgeid)
{
	// If either nodeid1 or nodeid is active or if they are the same node, return 0 (indicating the hyperedge is not subcritical).  
	if ((Node[nodeid1].Active == true) || (Node[nodeid].Active == true) || (nodeid1 == nodeid))
	{
		return 0;
	}
	// If none of the above conditions are met, return 1 (indicating the hyperedge becomes subcritical).  
	return 1;
}

// This function calculates the HCI1 for a single node.  
int Calculate_Node_HCI1(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge, int Nodeid)
{
	int iter_sum = 0; // Variable to record the count of subcritical paths of length 2.  
	// Iterate through the hyperedges connected to node Nodeid.  
	for (int j = 0; j < Node[Nodeid].Degree; j++)
	{
		// Check if the current hyperedge is subcritical.  
		if (Hyperedge[Node[Nodeid].Edge_id[j]].Subcritical == true)
		{
			// Iterate through the nodes in the current hyperedge.  
			for (int k = 0; k < Hyperedge[Node[Nodeid].Edge_id[j]].Cardinality; k++)
			{
				// Check if the path formed by the current node and the node in the hyperedge is subcritical.  
				int m = IsSubcritical(Node, Hyperedge, Hyperedge[Node[Nodeid].Edge_id[j]].Nodeinfo_id[k], Nodeid, Node[Nodeid].Edge_id[j]);
				iter_sum = iter_sum + m; // Add the result to the count of subcritical paths.  
			}
		}
	}
	return  Node[Nodeid].Degree + iter_sum; // Return the HCI1 value for the node, which is the sum of its hyperdegree and the count of subcritical paths.  
}

// This template function initializes the calculation of HCI1 values for all valid nodes in the hypergraph.  
// It takes two integers n1 and n2 as template parameters.
template<int n1, int n2>
void In_Calculate_HCI1(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge, int(&HCI1)[n1][n2])
{
	for (int i = 0; i < Max_Node_Size; i++)
	{
		// Set the node's index as its HCI1 for the first column.  
		HCI1[i][0] = i;
		// If the node has a non-zero degree, calculate its HCI1 value using the Calculate_Node_HCI1 function.  
		if (Node[i].Degree != 0)
		{
			HCI1[i][1] = Calculate_Node_HCI1(Node, Hyperedge, i);
		}
	}
}


// This template function re-calculates the HCI1 values for the 1-hop neighbors of newly activated nodes.  
// It takes two integers n3 and n4 as template parameters.  
template<int n3, int n4>
void ReCalculate_HCI1(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge, int(&HCI1)[n3][n4])
{
	int IterMark[Max_Node_Size] = { 0 }; // An array to mark if a node's HCI1 value has been updated.  
	while (ActiveNodeQueue.size() > 0) // Processes only the 1-hop neighbors of newly activated nodes once.  
	{
		int ActiveNode = ActiveNodeQueue[0]; // Retrieve the first node from the active node queue.  
		ActiveNodeQueue.pop_front(); // Remove the node from the queue.  
		for (int i = 0; i < Node[ActiveNode].Degree; i++) // Iterate through each neighbor of the active node.  
		{
			for (int j = 0; j < Hyperedge[Node[ActiveNode].Edge_id[i]].Cardinality; j++)
			{
				if (Hyperedge[Node[ActiveNode].Edge_id[i]].Nodeinfo_state[j] == 0 && IterMark[Hyperedge[Node[ActiveNode].Edge_id[i]].Nodeinfo_id[j]] == 0) // If the neighbor is inactive and its HCI1 value hasn't been updated, update it.  
				{
					HCI1[Hyperedge[Node[ActiveNode].Edge_id[i]].Nodeinfo_id[j]][1] = Calculate_Node_HCI1(Node, Hyperedge, Hyperedge[Node[ActiveNode].Edge_id[i]].Nodeinfo_id[j]); // Calculate and update the HCI1 value for the neighbor.  
					IterMark[Hyperedge[Node[ActiveNode].Edge_id[i]].Nodeinfo_id[j]] = 1; // Mark the neighbor as updated.  
				}
			}
		}
	}
}

// The HCI-TM algorithm to select Influence Maximization Set in a hypergraph.  
// It takes the arrays of nodes and hyperedges, the total number of nodes N, and the number of rounds as input.  
int HCI1_TM_Algorithm(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge, int N, int round)
{
	int HCI1[Max_Node_Size][2] = { 0 }; // A 2D array to store the IDs of nodes and their corresponding HCI1 values.  
	bool Seed[Max_Node_Size] = { false }; // An array to track which nodes are selected as seeds.  
	In_Calculate_HCI1(Node, Hyperedge, HCI1); // Call the function to initially calculate the HCI values for all nodes.  
	int iter_percent = 1; // The percentage interval to record the activation ratio.  
	int SeedNum = 0; // The count of selected seed nodes.  
	for (int i = 0; i < N; i++) // Iterate N times, selecting one node as a seed each time.  
	{
		int MaxCI = 0;
		int Seedid[1] = { 0 };
		for (int j = 0; j < Max_Node_Size; j++)
		{
			if (MaxCI < HCI1[j][1] && Node[HCI1[j][0]].Active == false) // Find the node with the highest CITM value that hasn't been activated yet.  
			{
				MaxCI = HCI1[j][1];
				Seedid[0] = HCI1[j][0];
			}
		}
		SeedNum++; // Increment the count of selected seed nodes.  
		Seed[Seedid[0]] = true; // Mark the selected node as a seed.  
		HCI1_Rec[SeedNum - 1] = MaxCI; // Record the maximum CITM value for the current seed.  
		SeedSet[0][SeedNum - 1] = Seedid[0]; // Add the selected node to the seed set.  
		Activation_rule(Node, Hyperedge, Seedid, 1); // Propagate activation based on the hypergraph threshold rule.  
		ActiveEdgeQueue.clear(); // Clear the active edge queue.  
		int iter_m = int(double(iter_percent * N) * Scale); // Calculate the number of seeds for the next recording interval.  
		if (SeedNum == iter_m) // If it's time to record the activation ratio, do so.  
		{
			RECORD[0][iter_percent - 1] = double(AC_NODE_NUM) / double(N);
			iter_percent = iter_percent + 1;
		}
		//if (AC_NODE_NUM >= int(double(N) * ActiveRadio)) // Check if the termination condition is met.
		if(SeedNum>=MAX_SEED_SET_SIZE)
		{
			RECORD[0][iter_percent - 1] = double(AC_NODE_NUM) / double(N);
			break; // Exit the loop if the condition is met.  
		}
		ReCalculate_HCI1(Node, Hyperedge, HCI1); // Recalculate the HCI values for the 1-hop neighbors of newly activated nodes.   
	}
	return SeedNum; // Return the total number of seed nodes selected.  
}

// Calculate the HCI2 value for a single node in the hypergraph.  
// The function takes the arrays of nodes and hyperedges, the node ID, and an array of seed statuses as input.  
template<int n9>
int Calculate_Node_HCI2(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge, int Nodeid, int(&Is_Seed)[n9])
{
	int iter_sum = 0; // The sum of length-2 subcritical paths.
	int iter_sum1 = 0; // The sum of length-3 subcritical paths.
	for (int j = 0; j < Node[Nodeid].Degree; j++) // Iterate through the hyperedges of node 'Nodeid'.
	{
		if (Hyperedge[Node[Nodeid].Edge_id[j]].Subcritical == true)
		{
			for (int k = 0; k < Hyperedge[Node[Nodeid].Edge_id[j]].Cardinality; k++)
			{
				int iter1 = Hyperedge[Node[Nodeid].Edge_id[j]].Nodeinfo_id[k]; // Iterate through the nodes in the hyperedge.  
				int m = IsSubcritical(Node, Hyperedge, iter1, Nodeid, Node[Nodeid].Edge_id[j]); // Determine if the path is subcritical.  
				iter_sum = iter_sum + m; // Update the sum of length-2 subcritical paths.  
				iter_sum1 = iter_sum1 + abs(Is_Seed[iter1] - 1) * m * (Node[iter1].Degree - 1); // Calculate the contribution of length-3 subcritical paths.  
			}
		}
	}
	return Node[Nodeid].Degree + iter_sum + iter_sum1; // Return the final HCI2 value.  
}


// Initialize the calculation of HCI2 values for all nodes in the hypergraph.  
// The function takes the arrays of nodes and hyperedges, and a 2D array to store the HCI2 values.  
// It also returns the number of valid nodes.  
template<int n1, int n2>
void In_Calculate_HCI2(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge, int(&HCI2)[n1][n2])
{
	int Is_Seed[Max_Node_Size] = { 0 }; // Initialize the seed status array.  
	for (int i = 0; i < Max_Node_Size; i++) // Iterate through each node in the hypergraph.  
	{
		HCI2[i][0] = i; // Set the node ID as the HCI2 value for the initial iteration.  
		if (Node[i].Degree != 0) // Check if the node has hyperedges.  
		{
			HCI2[i][1] = Calculate_Node_HCI2(Node, Hyperedge, i, Is_Seed); // Calculate the HCI2 value for the node.  
		}
	}
}



// Re-calculate the HCI2 values for the 1-layer neighbors of newly activated nodes.  
template<int n3, int n4>
void ReCalculate_HCI2(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge, int(&Node_CI1)[n3][n4], int(&Is_Seed)[n3])
{
	int IterMark[Max_Node_Size] = { 0 }; // Marker to ensure each node's HCI1 value is calculated only once.  

	while (ActiveNodeQueue.size() > 0) // Only calculate the HCI2 values for the 1-layer neighbors of newly activated nodes (once).  
	{
		int ActiveNode = ActiveNodeQueue[0]; // Get the first node from the queue.  
		ActiveNodeQueue.pop_front(); // Remove the first node from the queue.  

		for (int i = 0; i < Node[ActiveNode].Degree; i++) // Calculate the HCI2 value for this node's neighbors.  
		{
			for (int j = 0; j < Hyperedge[Node[ActiveNode].Edge_id[i]].Cardinality; j++)
			{
				// If the 1-layer neighbor is not activated and its CI value hasn't been updated.  
				if ((Hyperedge[Node[ActiveNode].Edge_id[i]].Nodeinfo_state[j] == 0) && (IterMark[Hyperedge[Node[ActiveNode].Edge_id[i]].Nodeinfo_id[j]] == 0))
				{
					// Calculate the HCI2 value for this neighbor and update its CI value.  
					Node_CI1[Hyperedge[Node[ActiveNode].Edge_id[i]].Nodeinfo_id[j]][1] = Calculate_Node_HCI2(Node, Hyperedge, Hyperedge[Node[ActiveNode].Edge_id[i]].Nodeinfo_id[j], Is_Seed);
					IterMark[Hyperedge[Node[ActiveNode].Edge_id[i]].Nodeinfo_id[j]] = 1; // Mark this neighbor as processed.  
				}
			}
		}
	}
}

// HCI-TM2 Algorithm  
int HCI2_TM_Algorithm(vector<NODE>& Node, vector<HYPEREDGE>& Hyperedge, int N, int round)
{
	int HCI2[Max_Node_Size][2] = { 0 }; // Array to store the IDs of nodes and their CITM values.  
	int iter_percent = 1;
	int SeedNum = 0;
	int Is_Seed[Max_Node_Size] = { 0 }; // Array to mark if a node is a seed.  
	In_Calculate_HCI2(Node, Hyperedge, HCI2); // Calculate the initial HCI values for all nodes.  

	for (int i = 0; i < N; i++) // Select a seed node each round.  
	{
		int MaxCI = 0;
		int Seedid[1] = { 0 };

		// Find the node with the highest CITM value that hasn't been activated yet.  
		for (int j = 0; j < Max_Node_Size; j++)
		{
			if (MaxCI < HCI2[j][1] && Node[HCI2[j][0]].Active == false)
			{
				MaxCI = HCI2[j][1];
				Seedid[0] = HCI2[j][0];
			}
		}

		SeedNum++; // Increment the number of seed nodes.  
		Is_Seed[Seedid[0]] = 1; // Mark the selected node as a seed.  
		HCI2_Rec[SeedNum - 1] = MaxCI; // Record the maximum CITM value for this round.  
		SeedSet[1][SeedNum - 1] = Seedid[0]; // Add the selected node to the seed set.  
		Activation_rule(Node, Hyperedge, Seedid, 1); // Propagate activation based on hypergraph threshold rules.  
		ActiveEdgeQueue.clear(); // Clear the active edge queue.  

		int iter_m = int(double(iter_percent * N) * Scale); // Calculate the number of seed nodes for the next round.  
		if (SeedNum == iter_m) // Check if it's time to record the activation ratio.  
		{
			RECORD[1][iter_percent - 1] = double(AC_NODE_NUM) / double(N);
			iter_percent = iter_percent + 1;
		}

		//if (AC_NODE_NUM >= int(double(N) * ActiveRadio)) // Check if enough nodes have been activated.
		if(SeedNum>=MAX_SEED_SET_SIZE)  
		{
			RECORD[1][iter_percent - 1] = double(AC_NODE_NUM) / double(N);
			break; // Exit the loop if enough nodes have been activated.  
		}
		ReCalculate_HCI2(Node, Hyperedge, HCI2, Is_Seed); // Recalculate HCI2 values for the 1-hop neighbors of newly activated nodes.   
	}
	return SeedNum; // Return the total number of seed nodes.  
}
int main()
{
	// File path for the hypergraph data  
	string File_path_name = "data/Algebra-question.txt";

	// Arrays to store the number of seed nodes and the running time of algorithms  
	double Record_Seed_Number[2] = { 0 };
	double Record_T[2] = { 0 };
	double Record_Seed_Radio[2] = { 0 };

	// Nodes and edges of the hypergraph  
	vector<NODE> Node(Max_Node_Size);
	vector<HYPEREDGE> Hyperedge(Max_Edge_Size);

	// Build the hypergraph and get the number of nodes  
	int N = ConstructHypergraph(Node, Hyperedge, File_path_name);

	// Output the file path, construction success message, and the number of nodes  
	cout << File_path_name << endl;
	cout << "Construct successfully!" << endl;
	cout << "The number of nodes are " << N << endl;

	// Start and end time for measuring algorithm runtime  
	double start = clock();

	// Run the HCI1_TM_Algorithm and store the number of selected seed nodes  
	int s1 = HCI1_TM_Algorithm(Node, Hyperedge, N, 1);
	double end = clock();

	// Output the running time of the HCI1-TM algorithm and the number of selected seed nodes  
	cout << "The running time of HCI1-TM algorithm: " << (end - start) / 1000 << " s" << endl;
	cout << "When 90% nodes in the hypergraph are activated, the number of seeds selected by HCI1-TM algorithm are:" << s1 << endl;
	Record_Seed_Number[0] = s1; // Store the number of selected seed nodes  
	Record_Seed_Radio[0] = double(s1) / N; // Calculate and store the ratio of selected seed nodes to total nodes  
	Record_T[0] = (end - start) / 1000; // Store the running time of the algorithm  

	// Call the HGRecover function to reset or recover the hypergraph state (possible implementation detail)  
	HGRecover(Node, Hyperedge);

	// Restart and end time for measuring algorithm runtime  
	start = clock();
	int s2 = HCI2_TM_Algorithm(Node, Hyperedge, N, 1); // Run the HCI2-TM algorithm and store the number of selected seed nodes  
	end = clock();

	// Output the running time of the HCI2-TM algorithm and the number of selected seed nodes  
	cout << "The running time of HCI2-TM algorithm: " << (end - start) / 1000 << " s" << endl;
	cout << "When 90% nodes in the hypergraph are activated, the number of seeds selected by HCI2-TM algorithm are:" << s2 << endl;
	Record_Seed_Number[1] = s2; // Store the number of selected seed nodes  
	Record_Seed_Radio[1] = double(s2) / N; // Calculate and store the ratio of selected seed nodes to total nodes  
	Record_T[1] = (end - start) / 1000; // Store the running time of the algorithm  

	// Call the HGRecover function again to reset or recover the hypergraph state (possible implementation detail)  
	HGRecover(Node, Hyperedge);
	cout << endl; // Output a new line for better formatting

	// Open a file stream to write to "D://QXL//Hypergraph-code//revise//test//Record-0.5.txt"  
	// This file will record the number of seed nodes, the running time of algorithms, and the ratio of seed nodes  
	ofstream out("output0.txt");

	// Loop through the array to record the results of two algorithms (HCI1-TM and HCI2-TM)  
	for (int i = 0; i < 2; i++)
	{
		// Write the results to the file  
		out << "HCI" << i + 1 << "-TM " << endl;
		out <<"The number of seed: " << Record_Seed_Number[i] << endl;
		out <<"The running time of the algorithm: " << Record_T[i] << endl;
		out<< "The proportion of seed nodes: "<<Record_Seed_Radio[i] << endl;
	}
	out.close(); // Close the file stream  



	// Loop through the RECORD array to smooth out any decrease in values  
	for (int i = 0; i < 2; i++)
	{
		for (int k = 1; k < 2500; k++)
		{
			if (RECORD[i][k] < RECORD[i][k - 1])
			{
				RECORD[i][k] = RECORD[i][k - 1];
			}
		}
	}

	// Open a file stream to write to "D://QXL//Hypergraph-code//revise//test//Record-Detail-Tr0.5.txt"  
	// This file will record detailed information about the activation process (ratio of seed nodes and activated nodes)  
	ofstream out1("output1.txt");

	// Loop through the RECORD array and write the detailed information to the file  
	for (int k = 0; k < 2500; k++)
	{
		out1 << double(double(k + 1) * Scale) << " " << RECORD[0][k] << " " << RECORD[1][k] << endl;
	}
	out1.close(); // Close the file stream  

	// Open a file stream to write to "D://QXL//Hypergraph-code//revise//test//Record-HCI-Tr0.5.txt"  
	// This file will record the evolution of the maximum HCI value in the hypergraph  
	ofstream out2("output2.txt");

	// Loop through the HCI1_Rec and HCI2_Rec arrays and write their values to the file  
	for (int j = 0; j < 20000; j++)
	{
		out2 << j << " " << HCI1_Rec[j] << " " << HCI2_Rec[j] << endl;
	}
	out2.close(); // Close the file stream  

	// Open a file stream to write to "D://QXL//Hypergraph-code//revise//test//Record-HCI1TM-SeedSet-Tr0.5.txt"  
	// This file will record the seed sets selected by the algorithms  
	ofstream out3("output3.txt");

	// Loop through the SeedSet arrays and write their values to the file, separated by spaces  
	for (int j = 0; j < s1; j++)
	{
		out3 << SeedSet[0][j] << " ";
	}
	out3.close(); // Close the file stream
	// Open a file stream to write to "D://QXL//Hypergraph-code//revise//test//Record-HCI2TM-SeedSet-Tr0.5.txt"  
	// This file will record the seed sets selected by the algorithms  
	ofstream out4("output4.txt");
	for (int j = 0; j < s2; j++)
	{
		out4 << SeedSet[1][j] << " ";
	}
	out4 << endl; // New line after the second set of seed nodes  
	return 0;
}

