#ifndef HASHTABLE_H_
#define HASHTABLE_H_

#include <iostream>
#include <string>
#include <cmath>
#include "cluster.h"

using namespace std; 

class LSHhashtable{

	public:
		float w;  
		int size,k,d;  
		vector< vector<class Point *> *> hashtable;  
		vector<class Point *> sPoints; 



		LSHhashtable(int size,int d,int k);  
		LSHhashtable(){};
		void createSpoints();
		~LSHhashtable(){} 
		
		int h(Point * x, Point * s);
		int calc_mod_power(int m, int power, int M);
		int g(Point * x);
		void insert_point(Point * p); 
		void insert_training(class PointSet & ps);
		float Manhattan_dist(Point  &x,Point  &y);
		Point * Nearest_Neighb(Point * q);
		int get_bucket(Point * p);
		vector<int> get_bucket_points_positions(int bucket_num, vector<int> & centroids);
};


#endif