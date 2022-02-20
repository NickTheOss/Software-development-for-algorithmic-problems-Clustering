#include "cluster.h"


PointSet ps;   
Point p;
lsh_Arguments arguments;  
vector<LSHhashtable *> hashtables; 
map<pair<int,int> , double> all_distances;  

void PointSet::insertPoint(Point * p){      
		p->pos_in_pointset = points.size();   
		points.push_back(p);
	}


void PointSet::deletePoints(){         

	int i;
	for(i = 0; i<points.size(); i++){
		delete points[i];
	}
	points.clear();

}

  
PointSet::~PointSet(){             

	deletePoints();

}


bool readFile_input(char * fileName){     
	ifstream input(fileName);             

	if(!input) {                          
		cout << "File Not Found"<<endl;
		return false;
	}


	char line[MAXSIZE];                 
	input.getline(line, MAXSIZE);
	line[7] = 0;
	if(strcmp(line,"vectors"))   
		exit(0);
	while((input.peek()!=EOF)) {        
		input.getline(line, MAXSIZE);   
		read_file_line_input(line);    
    }

    input.close();                 
    arguments.counter++;
   
	return true;            
}


bool readFile_config(char *fileName){      
	ifstream input(fileName);             

	if(!input) {                          
		cout << "File Not Found"<<endl;
		return false;
	}

	char line[MAXSIZE];                
	while((input.peek()!=EOF)) {        
		input.getline(line, MAXSIZE);   
		read_file_line_config(line);    
    }

    input.close();                 
    arguments.counter++;
   
	return true;                
}


void read_file_line_config(char * line){       
	char * token = strtok(line," \t\n");               
	if (strcmp(token,"number_of_clusters:") == 0){ 
		token = strtok(NULL," \t\n");
		arguments.c = atoi(token);         
	}

	if (strcmp(token,"number_of_grids:") == 0){ 
		token = strtok(NULL," \t\n");
		arguments.g = atoi(token);
	}

	if (strcmp(token,"number_of_vector_hash_tables:") == 0){ 
		token = strtok(NULL," \t\n");
		arguments.L = atoi(token);
	}

	if (strcmp(token,"number_of_vector_hash_functions:") == 0){ 
		token = strtok(NULL," \t\n");
		arguments.k = atoi(token);
	}
}


void read_file_line_input(char * line){      

	Point * temp = new Point; 
	char * token = strtok(line," \t\n");   
	if(arguments.counter != 1){
	    if((token!=NULL) && (strcmp(token,"curves") != 0)){                    
	    temp->name = token;        
	   	 
	   	} 
	   	while((token!=NULL)){     
	    	token = strtok(NULL," \t\n"); 
	  		if(token!=NULL)
	  		{             
	    		temp->coords.push_back(atof(token)); 
	    	}
		}

		ps.insertPoint(temp);    
	}	
}


void lsh_init(int c,char ** v){       
   
    int opt_value = 0;
    arguments.g = 2;            
    arguments.L = 3;  

   		while ((opt_value = getopt (c, v, "i:c:o:")) != -1){       
	     	switch(opt_value)
	     	{
	      	case 'i':                     
	          	arguments.input = optarg;  
	          	if(readFile_input(arguments.input) == true) continue; 
	          	else exit(0);
	      	case 'c':                      
	          	arguments.config = optarg; 
	          	if(readFile_config(arguments.config) == true) continue;  
	          	else exit(0);
	      	case 'o':                           
	          	arguments.output = optarg;    
	          	break;
     		}
		} 	
}


void mean_vector_update(vector<Cluster *> & clusters){  

	vector<Cluster *> new_clusters;

	for (int i = 0; i < clusters.size(); ++i)        
	{
		new_clusters.push_back(new Cluster());       
	} 

	avg_new_centroids(clusters, new_clusters); 

	Lloyds_assignment(new_clusters);    
	
	while(get_similarity(clusters,new_clusters) > 10){  

		assign_new_clusters( clusters, new_clusters); 
		avg_new_centroids(clusters, new_clusters);  
		Lloyds_assignment(new_clusters);      

	}
	assign_new_clusters( clusters, new_clusters);
	for (int i = 0; i < new_clusters.size(); ++i)
	{
		delete new_clusters[i];          
	}
}


int get_similarity(vector<Cluster *> & clusters, vector<Cluster *> & new_clusters){ 

	int counter = 0;
	for (int i = 0; i < clusters.size(); ++i) 
	{
		for (int j = 0; j < clusters[i]->team.size(); ++j) 
		{
			if ( find(new_clusters[i]->team.begin(), new_clusters[i]->team.end(), clusters[i]->team[j]) == new_clusters[i]->team.end() ) 
				counter++;
		}
	}
	return counter; 
}




 
void avg_new_centroids(vector<Cluster *> & clusters, vector<Cluster *> & new_clusters){

	for (int i = 0; i < clusters.size(); ++i)  
	{
		Point * new_avg_point = new Point();   

		for (int j = 0; j < ps.points[0]->coords.size(); ++j)  
		{
			new_avg_point->coords.push_back(ps.points[clusters[i]->centroid]->coords[j]); 

		}
		
		for (int j = 0; j < clusters[i]->team.size(); ++j)  
		{
			for (int z = 0; z < ps.points[0]->coords.size(); ++z)  
			{
				new_avg_point->coords[z] += ps.points[clusters[i]->team[j]]->coords[z]; 
			}			
		}
		for (int j = 0; j < ps.points[0]->coords.size(); ++j)
		{
			new_avg_point->coords[j] = new_avg_point->coords[j]/(clusters[i]->team.size()+1);  
		}
		new_avg_point->pos_in_pointset = ps.points.size();
		ps.points.push_back(new_avg_point);
		new_clusters[i]->centroid = ps.points.size()-1; 
	}
}





void assign_new_clusters(vector<Cluster *> & clusters, vector<Cluster *> & new_clusters){ 

	for (int i = 0; i < clusters.size(); ++i)  
	{
		delete clusters[i];    
		clusters[i] = new_clusters[i]; 
		new_clusters[i] = new Cluster;  
	}

	
	LSHhashtable ht;
	for (int i = 0; i < clusters.size(); ++i)
	{	
		if(ps.points[clusters[i]->centroid]->name == "none"){  
			if(clusters[i]->team.size() == 0){ 
				clusters[i]->centroid = 1; 
				continue;
			}
			
			float min = ht.Manhattan_dist(*ps.points[clusters[i]->centroid], *ps.points[clusters[i]->team[0]]);
			int new_centroid = clusters[i]->team[0];   
			float distance;

			for (int j = 1; j < clusters[i]->team.size(); ++j)  
			{	
				distance = ht.Manhattan_dist(*ps.points[clusters[i]->centroid],*ps.points[clusters[i]->team[j]]);
				if (min > distance && ps.points[clusters[i]->team[j]]->name != "none"){ 
					min = distance;  
					new_centroid = clusters[i]->team[j];  
				}
			}
			
			clusters[i]->centroid = new_centroid;  
			clusters[i]->team.erase(remove(clusters[i]->team.begin(), clusters[i]->team.end(), new_centroid), clusters[i]->team.end());
		}

	}
	
	for (int i = 0; i < clusters.size(); ++i)
	{
		for(int j = 0; j< clusters[i]->team.size(); j++){
			if(ps.points[clusters[i]->team[j]]->name == "none"){
				int avg_point = clusters[i]->team[j];
				clusters[i]->team.erase(remove(clusters[i]->team.begin(), clusters[i]->team.end(), avg_point), clusters[i]->team.end());
			}
		}
	}

	remove_fake_points(); 
}


void remove_fake_points(){

	for (int i = 0; i < arguments.c; ++i)   
	{
		delete ps.points[ps.points.size()-1];      
		ps.points.pop_back();	
	}
}




void Lloyds_assignment(vector<Cluster *> & clusters){ 

	int pos_cluster;
	for (int i = 0; i < ps.points.size(); i++)       
	{
		
		if(is_exist(i,clusters))     
			continue;
		pos_cluster = get_nearest_cluster(clusters,ps.points[i]);  
		clusters[pos_cluster]->team.push_back(i);   

	}
}


int get_nearest_cluster(vector<Cluster*> clusters,Point * x){  


	LSHhashtable ht;
	float min = ht.Manhattan_dist(*ps.points[clusters[0]->centroid],*x);  
	int cluster_pos = 0;
	float distance;
	for (int i = 1; i < clusters.size(); ++i) 
	{	
		distance = ht.Manhattan_dist(*ps.points[clusters[i]->centroid],*x); 
		if (min > distance){
			min = distance;
			cluster_pos = i;  
		}
	}
	return cluster_pos; 
}


void get_ht_results(vector<ht_results> & ht_r, int ht, vector<Cluster *> & clusters){  
	int flag;
	LSHhashtable * hashtable = hashtables[ht];  
	ht_results temp;  
	temp.clusters_pos.push_back(0);   
	temp.bucket_position = hashtable->get_bucket(ps.points[clusters[0]->centroid]); 
	ht_r.push_back(temp);
	
	for (int i = 1; i < clusters.size(); ++i)
	{
		flag = 1;
		temp.clusters_pos.clear();  
		temp.clusters_pos.push_back(i);
		temp.bucket_position = hashtable->get_bucket(ps.points[clusters[i]->centroid]);
		for (int j = 0; j < ht_r.size(); ++j)  
		{
			if(ht_r[j].bucket_position == temp.bucket_position){ 
				ht_r[j].clusters_pos.push_back(temp.clusters_pos[0]);
				flag = 0;
				break;
			}
		}
		if(flag){ 
			ht_r.push_back(temp);
		}
	}
}


void assign_points_to_clusters(vector<Cluster *> & clusters, vector<ht_results> & ht_r, int ht){

	vector<int> points;
	vector<int> centroids;
	for (int i = 0; i < clusters.size(); ++i)
	{
		centroids.push_back(clusters[i]->centroid);
	}
	for (int i = 0; i < ht_r.size(); ++i)
	{
		points.clear();
		points = hashtables[ht]->get_bucket_points_positions(ht_r[i].bucket_position, centroids); 
		
		if(ht_r[i].clusters_pos.size() == 1){
			Cluster * cluster = clusters[ht_r[i].clusters_pos[0]];
			for (int j = 0; j < points.size(); ++j)
			{	
				if(find(cluster->team.begin(), cluster->team.end(), points[j]) == cluster->team.end()){
					cluster->team.push_back(points[j]);
				}
			}
		}
		else { 
			int cluster_pos;
			vector<int> & clusters_pos = ht_r[i].clusters_pos;
			for (int j = 0; j < points.size(); ++j) 
			{
				cluster_pos = get_cluster_for_point(clusters, clusters_pos, points[j]);
				Cluster * cluster = clusters[cluster_pos];

				if(find(cluster->team.begin(), cluster->team.end(), points[j]) == cluster->team.end()){
					cluster->team.push_back(points[j]);
				}
			}
		}
	}
}


int get_cluster_for_point(vector<Cluster *> & clusters, vector<int> & clusters_pos, int point){ 

	LSHhashtable ht;
	float min = ht.Manhattan_dist(*(ps.points[clusters[clusters_pos[0]]->centroid]),*(ps.points[point]));
	int cluster_pos = clusters_pos[0];
	float distance;
	for (int i = 1; i < clusters_pos.size(); ++i)
	{	
		distance = ht.Manhattan_dist(*(ps.points[clusters[clusters_pos[i]]->centroid]),*(ps.points[point]));
		if (min > distance){
			min = distance;
			cluster_pos = clusters_pos[i];
		}
	}
	return cluster_pos;
}


void let_or_remove_point(vector<Cluster *> & clusters, int suggested, int othercluster, int point){

	LSHhashtable ht;
	float sug_dist = ht.Manhattan_dist(*ps.points[clusters[suggested]->centroid],*ps.points[point]);
	float oth_dist = ht.Manhattan_dist(*ps.points[clusters[othercluster]->centroid],*ps.points[point]);

	if(oth_dist < sug_dist)
		return;
	
	clusters[suggested]->team.push_back(point);
	clusters[othercluster]->team.erase(remove(clusters[othercluster]->team.begin(), clusters[othercluster]->team.end(), 5), clusters[othercluster]->team.end());
}


void put_point_in_proper_cluster(vector<Cluster *> & clusters, int suggestedCluster, int point){  
	int flag = 1, othercluster;
	for (int i = 0; i < clusters.size(); ++i)
	{
		if(i == suggestedCluster)
			continue;
		
		if(find(clusters[i]->team.begin(), clusters[i]->team.end(), point) != clusters[i]->team.end()){
			flag = 0;
			othercluster = i;
			break;
		}
	}
	
	if(flag){
		clusters[suggestedCluster]->team.push_back(point);
		return;
	}
	
	let_or_remove_point(clusters, suggestedCluster, othercluster, point);
}



void merge_clusters(vector<Cluster *> & clusters, vector<Cluster *> & clusters_temp){ 
	
	for (int i = 0; i < clusters_temp.size(); ++i)
	{
		
		for (int j = 0; j < clusters_temp[i]->team.size(); ++j)
		{
			
			if(find(clusters[i]->team.begin(), clusters[i]->team.end(), clusters_temp[i]->team[j]) != clusters[i]->team.end()){
				continue;
			}
			
			else{
				put_point_in_proper_cluster(clusters, i, clusters_temp[i]->team[j]);
			}
		}
	}
}



void clear_clusters_temp_teams(vector<Cluster *> & clusters_temp){

	for (int i = 0; i < clusters_temp.size(); ++i)
	{
		clusters_temp[i]->team.clear();        
	}
}


void init_clusters_temp(vector<Cluster *> & clusters_temp, vector<Cluster *> & clusters){

	for (int i = 0; i < clusters.size(); ++i)
	{
		clusters_temp.push_back(new Cluster(clusters[i]->centroid));  
	}
}


void LSH_assignments(vector<Cluster *> & clusters){

	vector<Cluster *> clusters_temp;
	init_clusters_temp(clusters_temp, clusters);   
	vector<ht_results> ht_r;
	
	get_ht_results(ht_r, 0, clusters);    
	assign_points_to_clusters(clusters, ht_r, 0);
	for (int i = 1; i < arguments.L; ++i)
	{
		ht_r.clear();
		get_ht_results(ht_r, i, clusters_temp); 
		assign_points_to_clusters(clusters_temp, ht_r, i);
		merge_clusters(clusters, clusters_temp);
		clear_clusters_temp_teams(clusters_temp);
	}


	for (int i = 0; i < clusters_temp.size(); ++i)
	{
		delete clusters_temp[i];   
	}

	
	vector<int> clusters_points;
	for (int i = 0; i < clusters.size(); ++i) 
	{
		clusters_points.insert(clusters_points.end(), clusters[i]->team.begin(), clusters[i]->team.end());
		clusters_points.push_back(clusters[i]->centroid); 
	}
	int pos;
	
	for (int i = 0; i < ps.points.size(); ++i)
	{
		if(find(clusters_points.begin(), clusters_points.end(), i) == clusters_points.end()){
			pos = get_nearest_cluster(clusters,ps.points[i]);
			clusters[pos]->team.push_back(i);
		}
	}
}



void Lloyds_update(vector<Cluster *> & clusters){  


	double min_sum_dist, sum_distance = 0.;
	int new_centroid;
	for (int i = 0; i < clusters.size(); ++i)
	{
		
		new_centroid = clusters[i]->centroid;  
		for (int j = 0; j < clusters[i]->team.size(); ++j)
		{
			sum_distance += all_distances[pair<int,int> (new_centroid,j)];  
		}
		min_sum_dist = sum_distance;  

		for (int j = 0; j < clusters[i]->team.size(); ++j)  
		{
			sum_distance = 0.;

			for (int z = 0; z < clusters[i]->team.size(); ++z)  
			{	
				if(clusters[i]->team[j] == clusters[i]->team[z])
					continue;
				sum_distance += all_distances[pair<int,int> (clusters[i]->team[j],clusters[i]->team[z])];
			}
			if(min_sum_dist > sum_distance){

				min_sum_dist = sum_distance;
				new_centroid = clusters[i]->team[j]; 
			}
		}
		if(clusters[i]->centroid != new_centroid){
			clusters[i]->team.push_back(clusters[i]->centroid); 
			clusters[i]->centroid = new_centroid;
			clusters[i]->team.erase(remove(clusters[i]->team.begin(), clusters[i]->team.end(), clusters[i]->centroid), clusters[i]->team.end()); 
		}

	}
}



float D(vector<Cluster*> clusters,Point * x){ 

	LSHhashtable ht;
	float min = ht.Manhattan_dist(*(ps.points[clusters[0]->centroid]),*x);
	float distance;
	for (int i = 0; i < clusters.size(); ++i)
	{	
		distance = ht.Manhattan_dist(*(ps.points[clusters[i]->centroid]),*x);
		if (min > distance)
			min = distance;
	}
	return distance;
}



int binarySearch(vector<content>  contents, int p, int r, float x) {  
   	if (p <= r) {    

     	int mid = (p + r)/2;   
    	float dif1 = fabs(x-contents[mid].partsum);

		float difprevious = fabs(contents[mid].partsum-contents[mid-1].partsum);
		float difnext = fabs(contents[mid].partsum-contents[mid+1].partsum);

		if(dif1 <= difprevious && dif1 <= difnext){
			return mid;	 
		}
		if (contents[mid].partsum > x)  
	         return binarySearch(contents, p, mid-1, x);   

	    if (contents[mid].partsum < x)  
	         return binarySearch(contents, mid+1, r, x); 
   } 
   return -1;  
} 



int get_new_centroid(vector<Cluster *> clusters){  

	
	vector<content> partial_sums; 
	content c;  
	c.partsum = 0;
	c.pos = -1;
	partial_sums.push_back(c);
	float min_dist, x;

	for (int i = 0; i < ps.points.size(); ++i)
	{
		
		if(is_exist(i, clusters))
			continue;

		min_dist = D(clusters, ps.points[i]);
		c.pos = i;
		c.partsum = min_dist*min_dist + partial_sums[partial_sums.size()-1].partsum; 
		partial_sums.push_back(c);  

	}
	int pos = random() % partial_sums.size();  
	x = partial_sums[pos].partsum;
	return binarySearch(partial_sums,0,partial_sums.size()-1,x);
}



vector<Cluster *> K_means_pp(){  


	vector<Cluster *> clusters;
	int position = rand()%ps.points.size();   
	clusters.push_back(new Cluster(position));
	while(clusters.size() < arguments.c){  

		position = get_new_centroid(clusters); 
		clusters.push_back(new Cluster(position));
		
	}
	return clusters;
}



bool is_exist(int pos,vector<Cluster *> vec){  

	for (int i = 0; i < vec.size(); ++i)   
	{
		if(pos == vec[i]->centroid) return true;  
	}

	return false;

}



vector<Cluster *> random_Kpoints(){  
	
	vector<Cluster *> clusters;
	while(clusters.size() < arguments.c){

		int position = rand()%ps.points.size();
		if(is_exist(position,clusters)) 
			continue;
		clusters.push_back(new Cluster(position));
		
	}
	return clusters;
}


 
void LSH_calculation(){   

	Point * temp;        
	int i, j;       

	for(i = 0; i<arguments.L; i++){    
		hashtables.push_back(new LSHhashtable(ps.points.size()/8, ps.points[0]->coords.size(),arguments.k));
		hashtables[i]->insert_training(ps);  
		cout<<"Hashtable "<<i+1<<" inserted"<< endl;   
	}
}


vector<Cluster *> pattern1(){         //1 1 1 

	vector<Cluster *> clusters;   
	clusters = random_Kpoints();
	Lloyds_assignment( clusters );
	Lloyds_update( clusters );
	return clusters;
}


vector<Cluster *> pattern2(){   //1 2 1

	vector<Cluster *> clusters;
	clusters = random_Kpoints();
	LSH_assignments( clusters );
	Lloyds_update( clusters );
	return clusters;
}


vector<Cluster *> pattern3(){   //1 2 2

	vector<Cluster *> clusters;
	clusters = random_Kpoints();
	LSH_assignments( clusters );
	mean_vector_update( clusters );
	return clusters;
}


vector<Cluster *> pattern4(){   //1 1 2

	vector<Cluster *> clusters;
	clusters = random_Kpoints();
	Lloyds_assignment( clusters );
	mean_vector_update( clusters );
	return clusters;
}


vector<Cluster *> pattern5(){   //2 1 1 

	vector<Cluster *> clusters;
	clusters = K_means_pp();
	Lloyds_assignment( clusters );
	Lloyds_update( clusters );
	return clusters;
}


vector<Cluster *> pattern6(){   //2 1 2

	vector<Cluster *> clusters;
	clusters = K_means_pp();
	Lloyds_assignment( clusters );
	mean_vector_update( clusters );
	return clusters;
}


vector<Cluster *> pattern7(){   //2 2 1

	vector<Cluster *> clusters;
	clusters = K_means_pp();
	LSH_assignments( clusters );
	Lloyds_update( clusters );
	return clusters;
}


vector<Cluster *> pattern8(){   //2 2 2

	vector<Cluster *> clusters;
	clusters = K_means_pp();
	LSH_assignments( clusters );
	mean_vector_update( clusters );
	return clusters;
}


void calculate_all_distances(){ 

	LSHhashtable ht;
	double dist;
	for (int i = 0; i < ps.points.size(); ++i)
	{
		for (int j = i+1; j < ps.points.size(); ++j)
		{
			dist = ht.Manhattan_dist(*ps.points[i],*ps.points[j]);
            all_distances[pair<int,int> (i,j)] = dist ;
            all_distances[pair<int,int> (j,i)] = dist ;
		}
	}
}



double a_i(Cluster * cluster, int point){  

	double avg_distance = all_distances[pair<int,int> (cluster->centroid,point)];

	for (int i = 0; i < cluster->team.size(); ++i)
	{
		
		avg_distance += all_distances[pair<int,int> (cluster->team[i],point)];

	}

	return cluster->team.size() > 0 ? avg_distance/cluster->team.size() : 0.;
}


double b_i(vector<Cluster *> & clusters, int cluster_pos, int point){  

	double min = cluster_pos > 0 ?  all_distances[pair<int,int> (clusters[0]->centroid, point)] : all_distances[pair<int,int> (clusters[1]->centroid,point)];
	double distance;
	int nearest_cluster = cluster_pos > 0 ? 0 : 1;
	
	for (int i = 0; i < clusters.size(); ++i)
	{
		
		if(i == cluster_pos)
			continue;
		distance = all_distances[pair<int,int> (clusters[i]->centroid, point)];
		if(min > distance){

			min = distance;
			nearest_cluster = i;

		}
	}
	distance = min; 
	for (int i = 0; i < clusters[nearest_cluster]->team.size(); ++i)
	{
		
		distance += all_distances[pair<int,int> (clusters[nearest_cluster]->team[i], point)];

	}
	return clusters[nearest_cluster]->team.size() > 0 ? distance/clusters[nearest_cluster]->team.size() : 0.;
}



double s_i(vector<Cluster *> & clusters, int cluster_pos, int point){  

	double ai = a_i(clusters[cluster_pos],point);
	double bi = b_i(clusters,cluster_pos,point);
	double max = ai > bi ? ai : bi;
	return max == 0 ? 0 : (bi-ai)/max;

}



void silhouette(FILE * fp, vector<Cluster *> & clusters){ 

	double sc = 0., sc_total = 0.;
	fprintf(fp,"Silhouette: [");
	for (int i = 0; i < clusters.size(); ++i)
	{
		sc += s_i(clusters, i, clusters[i]->centroid);

		for (int j = 0; j < clusters[i]->team.size(); ++j)
		{

			sc += s_i(clusters, i, clusters[i]->team[j]);
		}
		sc_total += sc;
		fprintf(fp,"%lf, ",sc/(clusters[i]->team.size()+1));
		sc = 0.;
	}
	fprintf(fp,"%lf]\n",sc_total/ps.points.size());
}


void clusters_output(FILE * fp ,vector<Cluster *> & clusters, double algo_time, string algo_name ){ 

	fprintf(fp,"Algorithm: %s\n", algo_name.c_str());
	for (int i =0; i < clusters.size(); ++i)
	{
		
		fprintf(fp,"CLUSTER-%d size: %ld, centroid: %s\n",(i+1),(clusters[i]->team.size()+ 1), ps.points[clusters[i]->centroid]->name.c_str());

	}
	fprintf(fp,"clustering_time: %lf\n",algo_time/CLOCKS_PER_SEC);
	silhouette(fp, clusters);
}

void free_Clusters(vector<Cluster *> & clusters){  
	for (int i = 0; i < clusters.size(); ++i)
	{
		delete clusters[i];
	}
}

int main( int argc, char ** argv){             

	srand(time(NULL));
	arguments.counter = 0; 
	lsh_init(argc, argv);
	fflush(stdin);  
	int pos = 0;
	FILE * fp =fopen(arguments.output.c_str(),"w+"); 
	clock_t c1,c2;  
	LSH_calculation();  
	calculate_all_distances(); 
	vector<Cluster *> clusters;
	c1 = clock();
	clusters = pattern1();
	c2 = clock();
	clusters_output(fp , clusters, c2-c1, "I1A1U1" );
	free_Clusters(clusters);  
	//algorithm = "1 1 1";
	c1 = clock();
	clusters = pattern2();
	c2 = clock();
	clusters_output(fp , clusters, c2-c1, "I1A2U1" );
	free_Clusters(clusters);
	//algorithm = "1 2 1";
	c1 = clock();
	clusters = pattern3();
	c2 = clock();
	clusters_output(fp , clusters, c2-c1, "I1A2U2" );cout << ps.points.size() << endl;
	free_Clusters(clusters);
	//algorithm = "1 2 2";
	c1 = clock();
	clusters = pattern4();
	c2 = clock();
	clusters_output(fp , clusters, c2-c1, "I1A1U2" );cout << ps.points.size() << endl;
	free_Clusters(clusters);
	//algorithm = "1 1 2";
	c1 = clock();
	clusters = pattern5();
	c2 = clock();
	clusters_output(fp , clusters, c2-c1, "I2A1U1" );
	free_Clusters(clusters);
	//algorithm = "2 1 1";
	c1 = clock();
	clusters = pattern6();
	c2 = clock();
	clusters_output(fp , clusters, c2-c1, "I2A1U2" );
	free_Clusters(clusters);
	//algorithm = "2 1 2";
	c1 = clock();
	clusters = pattern7();
	c2 = clock();
	clusters_output(fp , clusters, c2-c1, "I2A2U1" );
	free_Clusters(clusters);
	//algorithm = "2 2 1";
	c1 = clock();
	clusters = pattern8();
	c2 = clock();
	clusters_output(fp , clusters, c2-c1, "I2A2U2" );
	free_Clusters(clusters);
	//algorithm = "2 2 2";
	
	cout<<"Program now terminates.............."<<endl; 
	fclose(fp); 

	for(int i = 0; i<arguments.L; i++){    
		delete hashtables[i];
	}

}
