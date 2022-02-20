#include "hashtable.h"
#include "cluster.h"



LSHhashtable::LSHhashtable(int size,int d,int k){ 

	this->size = size; 
	this->k = k;  
	this->w = 4500; 
	this->d = d;   
	int i;
	for(i = 0;i<size;i++){   
		hashtable.push_back(new vector<Point *>);  
	}

	srand(time(NULL));   
	createSpoints();   
}


void LSHhashtable::createSpoints(){

	vector<Point *> temp_sPoints;  
	int i,j;
	Point * p;            
	int ww = w; 
	ww--; 
	for(i=0;i<5*k;i++){   
		p = new Point; 
		for(j=0;j<d;j++){    
			p->coords.push_back(((float) rand()/RAND_MAX)+ rand()%ww);   // (pv+t/w) 
		}

		temp_sPoints.push_back(p); 

	}
	for(i = 0;i<k;i++){    
		int pos = rand()%(5*k);  
		sPoints.push_back(temp_sPoints[pos]); 
		temp_sPoints[pos]->name = "valid";
	}
	string vld = "valid";
	for(i = 0;i<temp_sPoints.size();i++){
		if(temp_sPoints[i]->name.compare(vld) != 0){  
			delete temp_sPoints[i]; 
		}
	}
}




int LSHhashtable::h(Point * x, Point * s){ // h h() g(x)=h1|h2|h3...
	Point a;  
	int sum = 0 , M = pow(2, 32/k);   //M = 2^32-5
	int i, j, num, m = (x->coords[0] - s->coords[0])/w;  //a = (xi-si/w) 
	
	for(i = 0 ; i<d; i++){  

		num = ((x->coords[i] - s->coords[i])/w) * 100;  
		a.coords.push_back(num); 
		if(m < a.coords[i]) m = a.coords[i];   
	}
	m+=10; 

	for(i = a.coords.size()-1,j = 0 ; i>=0; i--,j++){   //j =power nnlsh
		int temp = a.coords[i]; 
		sum += (temp%M * calc_mod_power(m, j, M))%M;   
	}
	return sum%M;  //h(x) = ad􀀀1 + m  ad􀀀2 +    + md􀀀1  a0 mod M 2 N;
}




int LSHhashtable::calc_mod_power(int m, int power, int M){  

	if(power == 2){
		return ((m%M)*(m%M))%M;
	}
	else if(power == 1) return m;
	else if(power == 0) return 1;

	return (calc_mod_power(m,power-1,M) * (m % M)) % M; 

}




int LSHhashtable::g(Point * x){ 

	
	int i ,sum = h( x, sPoints[0]);

	for(i = 1; i<k;i++){  
		sum |= h( x, sPoints[i]); //OR level

	}
	return sum%hashtable.size();

}






void LSHhashtable::insert_point(Point * p){


	int position = g(p);    
	hashtable[position]->push_back(p); 

}


int LSHhashtable::get_bucket(Point * p){

	return g(p);

}


vector<int> LSHhashtable::get_bucket_points_positions(int bucket_num, vector<int> & centroids){


	vector<int> vec;
	for (int i = 0; i < hashtable[bucket_num]->size(); ++i)
	{
		
		if(find(centroids.begin(), centroids.end(), hashtable[bucket_num]->at(i)->pos_in_pointset) != centroids.end()){
			continue;
		}

		vec.push_back(hashtable[bucket_num]->at(i)->pos_in_pointset); 

	}
	return vec;
}




void LSHhashtable::insert_training(PointSet & ps){  

	int i;
	for(i = 0; i<ps.points.size() ; i++){
		insert_point(ps.points[i]); 

	}



}



float LSHhashtable::Manhattan_dist(Point  &x,Point  &y){ 

	float distance = 0.;
	int d = x.coords.size();
	int i;
	for(i = 0; i< d ;i++){ 
		distance += abs(x.coords[i]-y.coords[i]); 
		
	}
	return distance;

}





Point * LSHhashtable::Nearest_Neighb(Point * q){

	int bucket = g(q);  
	int i;
	Point * nn;
	float distance;
	float min_distance;
	if(hashtable[bucket]->size() == 0) return NULL;  
	nn = hashtable[bucket]->at(0); 
	min_distance = Manhattan_dist( *q, *nn); 

	for(i = 1 ;i<hashtable[bucket]->size();i++){ 
		distance = Manhattan_dist(*q,* (hashtable[bucket]->at(i)));  
		if(min_distance > distance ) {  
			min_distance = distance;   
			nn = hashtable[bucket]->at(i);  
		}

	}
	return nn; 
}