//reads a pdb file and returns the minimal distance between sidechain atoms for each residue
//residues: variable which stores the 1st and last atom number for each residue
//atoms: xyz coordinates of all atoms
//the program assumes clean PDBS, i.e. one chain, no dual residue definitions, etc.
//max_dist defines the maximal distance between two residue in order to be defined as a contact


#include <string>
#include <iostream>
#include <fstream>
#include <new>
#include <vector>
#include <math.h>
#include <time.h>
#include <dirent.h>

using namespace std;

int contact_map(std::string pdb_filename){
    clock_t tic = clock();


//defines all variables
	string line;
	int atom_no=0;
	vector<string> pdb;
	int residue=0;
	int first_residue=10000;
	int last_residue=-100;
	//int sidechain=0;
	int start_atom=-1;
	int end_atom=-1;
	float dist;
	float best_dist;
	//maximal distance to be considered as a contact
	int max_dist=5;
	
	struct point{float x; float y; float z;};
	point mypoint;
	
	
	vector<string> atoms;
	
//opens pdb file
	std::string map_filename=pdb_filename.substr(0,pdb_filename.length()-4)+".txt";

	ifstream pdb_file;
	pdb_file.open(pdb_filename);
	

//reads through pdb file and write ATOM line into pdb vector
	while ( getline (pdb_file,line) ){
		if (line.substr(0,6)=="ATOM  "){
			atom_no++;
			pdb.push_back(line);
		}
		else if (line.substr(0,7)=="ENDMDL ")
			break;	
    }
    for (int i=0;i<20;i++)
    	std::cout <<pdb[i];
    pdb_file.clear();
    pdb_file.close();
    
    
	
	//pdb.clear();
	//mypoint={1,2,3};
	
    for (int j=0; j<atom_no; j++){
    	
    	//std::cout <<i<<pdb[i]<<std::endl;
		//mypoint={stof(pdb[i].substr(30,8)),stof(pdb[i].substr(38,8)),stof(pdb[i].substr(46,8))};	
		

	
		//std::cout <<std::to_string(mypoint.z)<<std::endl;
		std::cout<<j<<std::endl;
		atoms.push_back (std::to_string(1));
		
    }
    std::cout <<"bla";

ofstream map_file;
map_file.open(map_filename);    

//calculates distance between every atom
	for (int i=0;i<last_residue-first_residue-2;i++){
		for (int j=1;j<last_residue-first_residue-1;j++){
			if (j>i){
				dist=0;
				best_dist=max_dist;
				
				if (best_dist<max_dist){
					line=std::to_string(i+first_residue)+";"+std::to_string(j+first_residue)+";"+std::to_string(best_dist)+"\n";
					map_file <<line;
				}
			}
		}
	}
	std::cout <<"bla";

	
	
	
	
	//atoms.clear();	
	pdb.clear();
	map_file.close();
	clock_t toc = clock();
	std::cout <<pdb_filename+" finished in ";
	std::cout <<double(toc-tic);
	std::cout <<" msecs"<<std::endl;
	return 0;
}


int main()
{
	std::string temp_filename="";
	DIR *files;
	struct dirent *ent;
	std::string pdb_file="";
	std::string pdb_directory="D:/My documents/C++/chains/";

	

	contact_map("D:/My documents/C++/chains/1FXL_A.pdb");
	contact_map("D:/My documents/C++/chains/1B35_B.pdb");

	return 0;
		
	files=opendir(pdb_directory.c_str());
	while ((ent = readdir (files)) != NULL) {
		pdb_file=ent->d_name;
		std::cout <<pdb_file;
		if (pdb_file.length()>4 && pdb_file.substr(pdb_file.length()-4,4)==".pdb"){
		
	    	//printf ("%s\n", ent->d_name);
	    	std::cout << "     ";
	    	temp_filename=pdb_directory+pdb_file.substr(0,pdb_file.length()-4)+".txt";
	    	std::cout <<temp_filename<<std::endl;
	    	ifstream temp(temp_filename.c_str());
	    	if (temp.good())
	    		std::cout <<"already done"<<std::endl;
	    	else{
	    		contact_map(pdb_directory+ent->d_name);
	    		return 0;
	    	}
    	}
  	}

	closedir(files);
	//contact_map("D:/My documents/C++/chains/1FXL_A.pdb");
}
