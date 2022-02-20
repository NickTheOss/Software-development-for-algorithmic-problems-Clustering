#ifndef CLUSTER_H_
#define CLUSTER_H_

#include <iostream>         
#include <string>          
#include <vector>			
#include <algorithm>
#include <string.h>         
#include <stdlib.h>        
#include <stdio.h>         
#include <stdbool.h>       
#include <unistd.h>         
#include <time.h>          
#include <fstream>          
#include <map>
#include "hashtable.h"


#define MAXSIZE 5000        
#define BUFFSIZE 50         

using namespace std;  


class Point{        

public:
	string name;       
	vector<float> coords;   
	int pos_in_pointset;	
	Point():name("none"){}             
};


class PointSet{           

public:
	
	vector<Point *> points;       
	~PointSet();                 
	void deletePoints();         
	void insertPoint(Point * );	  
};




class lsh_Arguments{  

public:
	char * input;    
	char * config;           
	string output;   
	int c, L, k, g, counter;  
	
};




class Cluster{

public:
	int centroid;
	vector<int> team;  
	Cluster(int pos): centroid(pos){}
	Cluster(){}  
};


typedef struct{  
		
		float partsum;
		int pos;

}content;




typedef struct{ 

	int bucket_position;
	vector<int> clusters_pos;

}ht_results;




bool readFile_input(char *);
bool readFile_config(char *);
void read_file_line_input(char *);
void read_file_line_config(char *);
void lsh_init(int ,char **);
void LSH_calculation();
vector<Cluster *> random_Kpoints();
bool is_exist(int pos,vector<Cluster*> vec);
float D(vector<Cluster*> clusters,Point * x);
vector<Cluster *> K_means_pp();
int get_new_centroid(vector<Cluster *>);
int binarySearch(vector<content>  contents, int p, int r, float x);
int get_nearest_cluster(vector<Cluster*> clusters,Point * x);
void Lloyds_assignment(vector<Cluster *> & clusters);
void LSH_assignments(vector<Cluster *> & clusters);
void get_ht_results(vector<ht_results> & ht_r, int ht, vector<Cluster *> & clusters); 
void assign_points_to_clusters(vector<Cluster *> & clusters, vector<ht_results> & ht_r, int ht); 
int get_cluster_for_point(vector<Cluster *> & clusters, vector<int> & clusters_pos, int point);
void merge_clusters(vector<Cluster *> & clusters, vector<Cluster *> & clusters_temp); 
void clear_clusters_temp_teams(vector<Cluster *> & clusters_temp); 
void init_clusters_temp(vector<Cluster *> & clusters_temp, vector<Cluster *> & clusters); 
void put_point_in_proper_cluster(vector<Cluster *> & clusters, int suggestedCluster, int point);
void let_or_remove_point(vector<Cluster *> & clusters, int suggested, int othercluster, int point);
void calculate_all_distances(); 
void Lloyds_update(vector<Cluster *> & clusters);
void mean_vector_update(vector<Cluster *> & clusters);
void avg_new_centroids(vector<Cluster *> & clusters, vector<Cluster *> & new_clusters);
int get_similarity(vector<Cluster *> & clusters, vector<Cluster *> & new_clusters);
void remove_fake_points();
void assign_new_clusters(vector<Cluster *> & clusters, vector<Cluster *> & new_clusters);
vector<Cluster *> pattern1();
vector<Cluster *> pattern2();
vector<Cluster *> pattern3();
vector<Cluster *> pattern4();
vector<Cluster *> pattern5();
vector<Cluster *> pattern6();
vector<Cluster *> pattern7();
vector<Cluster *> pattern8();
void clusters_output(FILE * fp ,vector<Cluster *> & clusters, double algo_time, string algo_name );
double a_i(Cluster * cluster, int point);
double b_i(vector<Cluster *> & clusters, int cluster_pos, int point);
double s_i(vector<Cluster *> & clusters, int cluster_pos, int point);
void silhouette(FILE * fp, vector<Cluster *> & clusters);
void free_Clusters(vector<Cluster *> & clusters);
#endif