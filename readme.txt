The folder includes 2 source code files (cluster.cpp, hashtable.cpp), 2 header files (cluster.h, hashtable.h) along with 1 .txt file with the requested results of the algorithms and the makefile of the program plus one of the indicative test input files provided.
The first part of the work, which concerned the grouping of vectors (Vectors), was completed.
The Cluster folder includes:

cluster.cpp ----> includes the constructor and destructor of the PointSet point table (Here we keep the coordinates of each point read). First we call the lsh_init() function which takes as arguments the number of arguments from the command line and the value of each string stored in the argv table. It initializes L hashtables, the number of hashfunctions and anything else required. Then, after the command prompt input, it inserts the required filenames into the corresponding cluster.h structure. Lsh_init() calls readFile(), which takes as its argument the file name and an object of the Pointset class (to fill an array of input points accordingly). After we open the file and check if it has a wrong name, through input.peek() we read line by line each point (up to the blank line) and pass it to the function read_file_line_input(), read_file_line_config() (for the input file and config respectively ) which takes as an argument each line we read (that is, each point). Read_file_line() (for input and config respectively)  "breaks" each line into pieces via strtok with space, tab and endline as separators. If the first token is not equal to "vectors", then we store each token in the Vector of Point and then we insert the point in the table from PointSet tables. Then we close each file and then the LSH_calculation() function is called to initialize the hashtables.
LSH_calculation briefly allocates memory for the number of hashtables that will be created and inserts in them each point through the insert_training() function that we will analyze below, finally releasing the hashtables from memory. The initial idea of ​​the algorithmic problem starts with the function calculate_all_distances() (line 784) which calculates the distances of all points of the dataset between them. Then we call 8 functions to cover all cases of interaction of the algorithms for initialization, assignment and update by keeping in 2 variables and the clustering time and finally freeing the memory from the clusters variable which is essentially the table of clusters grouped as each pattern function returns.
The following is a brief description of the cluster.cpp functions:

a) random_Kpoints () (line 675) -> initializes clusters by randomly selecting centroids from the database.
b) is_exist () (line 662) -> finds if one position in the cluster is the centroid position by returning an appropriate value.
c) D () (line 577) -> Returns the minimum distance of a point from the Clusters.
d) K_means_pp () (line 645) -> initializes clusters to centroids based on the probability algorithm of the theory.
e) get_new_centroid () (line 616) -> Returns the position of the table where we will insert the centroid.
f) binarySearch () (line 593) -> Binary search to find the position of the centroid.
g) get_nearest_cluster () (line 310) -> Returns where a point in a cluster is located.
h) Lloyds_assignment () (line 295) -> Creates the Clusters group based on the LLoyds_assignment algorithm.
i) LSH_assignments () (line 492) -> Creates the Clusters group based on the LSH_assignments algorithm.
j) get_ht_results () (line 329) -> Fills the ht_r table with bucket number and clusters whose centroids are in buckets.
k) assign_points_to_clusters () (line 358) -> Distributes ht_r points in clusters, making sure they are unique.
l) get_cluster_for_point () (line 397) -> Returns the cluster number to the clusters to which the point belongs.
m) merge_clusters () (line 453) -> We put the points in the teams of each cluster without duplicates.
n) clear_clusters_temp_teams () (line 474) -> Clear the temps from the clusters_temp, leaving the centroids.
o) init_clusters_temp () (line 483) -> Initiates clusters_temp with only the centroids of the clusters.
p) put_point_in_proper_cluster () (line 429) -> Puts a dot in the correct cluster.
q) let_or_remove_point () (line 415) -> To determine if a Point should be moved from one cluster to another according to the algorithms.
r) Lloyds_update () (line 535) -> Upgrading centroids according to LLoyd.
s) mean_vector_update () (line 152) -> Upgrading clusters by creating fake points according to theory.
t) avg_new_centroids () (line 198) -> Calculates the new centroids and adds them to the end of the ps.points table.
u) get_similarity () (line 180) -> Returns an index of similarity of the new clusters with the old ones according to some features. If the similarity index is conventionally below 10, then appropriate actions are taken (see implementation).
v) remove_fake_points () (line 283) -> Removes the fake points that are not needed.
w) assign_new_clusters () (line 231) -> Rearranges the clusters in such a way that their new copies are created correctly.
x) pattern () 1-8 (line 704-780) -> They call the basic functions init, assign, update with different sequences at a time for each algorithm.
y) silhouette () (line 858) -> "Evaluates" how efficient the algorithms are each time for each sequence using the functions a_i, b_i, s_i (Respectively lines 801, 816, 847).
z) free_Clusters () (line 892) -> Releases cluster tables from memory.

cluster.h ----> Includes all the necessary libraries needed for the whole program along with the cluster.cpp function statements and the corresponding classes. The Point class is the data of each point, it's name and it's coordinates together with the minimum distance described below. The PointSet class contains a table of Point type pointers and lsh_Arguments holds the input files and the corresponding values. In addition we have the Cluster class that contains the positions in the corresponding Point table of the points that are clustered and 2 corresponding structures, the content that contains the attributes of the table with the partial sums and the ht_results that contains the positions of the buckets and clusters for LSH_assignments algorithm's needs.

hashtable.cpp ----> Includes the Hashtable constructor, with the initialization of it's features, the allocation of memory locations and the initialization of the starting points through the corresponding function. Create_sPoints() briefly creates a table of points, where by various transformations of creating pseudo-random numbers (based on the theory) it randomly inserts various values ​​into it (for a number 5 times greater than the number of "small" hashfunctions ).It then frees up table locations that are not needed so that there is no unnecessary free memory. The function h takes 2 Point type pointers and based on the theory makes transformations to calculate the equation -h (x) = ad-1 + m * ad-2 ... + m ^ d-1a0modM- calling the auxiliary function calc_mod_power() and returns (sum % M). Then we have g, which is the basic hashfunction and returns the (concatenation_of_h) Mod (array size), ie the final position where we will enter the point. Insert_point() puts a dot in the appropriate position and insert_training() puts all points in the hashtable using insert_point(). Then we have Manhattan_dist which takes 2 points and calculates their Manhattan distance assigning them to a final distance sum. Finally Nearest_Neighb() finds the nearest neighbor of each query point, checking it's buckets in each hashtable position with additional functions for it's needs( the get_bucket() and get_bucket_points_positions() that assign clusters and points according to the theory and LSH algorithm).

hashtable.h ----> Includes hashtable attribute statements, along with all functions used in hashtable.cpp.


Command for compiling = 'make' and = 'make clean' for cleaning object files of the program
Command line command with arguments: './cluster -i test.csv -c cluster.conf -o output'
Input files: 'test.csv' 'cluster.conf'
Output file: 'output'

Implementation Language: C++
Implementation Environment: Linux