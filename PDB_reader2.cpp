//reads a pdb file and saves its atoms in an object
//if multiple models are present, only the first one is read

#include <string>
#include <iostream>
#include <fstream>
#include <new>
#include <vector>
#include <math.h>
#include <time.h>
#include <dirent.h>

using namespace std;

int pdb_object(std::string pdb_filename){
    clock_t tic = clock();


//defines all variables
    string line;
    int atom_no=0;
    vector<string> pdb;
    int residue=0;
    int temp_residue=0;
    int first_residue=10000;
    int last_residue=-100;
    //int sidechain=0;
    int start_atom=-1;
    int end_atom=-1;
    float dist;
    float best_dist;
    int max_dist=5;
    int residue1_end=-1;
    int residue2_end=-1;
    int j=0;
    int jj=0;
    int ii=0;   
            
    struct point{float x; float y; float z;};
    struct atom{string record;string atomtype; string residuetype;string chain;int residue;point location;};
    
    point mypoint;
    atom myatom;
    
    
    vector<atom> atoms;
    
//opens pdb file
    ifstream pdb_file;
    pdb_file.open(pdb_filename);
    

//reads through pdb file and write ATOM line into pdb vector
    while ( getline (pdb_file,line) ){
        if (line.substr(0,6)=="ATOM  "||line.substr(0,6)=="HETATM"){
            atom_no++;
            myatom.record=line.substr(0,6);
            myatom.atomtype=line.substr(12,4).c_str();
            myatom.residuetype=line.substr(17,3).c_str();
            myatom.chain=line.substr(21,1).c_str();
            myatom.residue=stoi(line.substr(22,4));
            mypoint={stof(line.substr(30,8)),stof(line.substr(38,8)),stof(line.substr(46,8))};
            myatom.location=mypoint;
            atoms.push_back(myatom);
            pdb.push_back(line);
        }
        else if (line.substr(0,7)=="ENDMDL ")
            break;  
    }

    pdb_file.clear();
    pdb_file.close();
    
    
    //measure distance between all atoms of all residues, returns the shortest distance between a residue pair
    residue=atoms[0].residue-1;
    for (int i=0;i<atom_no;i++){
        residue2_end=-1;
        //if new residue starts, start calculating distances
        if (atoms[i].record=="ATOM  " &&residue!=atoms[i].residue){
            residue=atoms[i].residue;
            j=i;
            
            //finds the last atom which belongs to the first residue[i]
            while (j<atom_no){
                j++;
                if (residue!=atoms[j].residue && atoms[i].record=="ATOM  "){
                    residue1_end=j-1;
                    j=atom_no+1;                        
                }
            }
            
            //if not the last residue was taken, continue to find a partner to compare to
            if (j!=atom_no){
                residue=atoms[residue1_end+1].residue;
                ii=residue1_end+1;
                if (residue2_end==-1)
                    residue2_end=residue1_end+1;
                    
                
                //continues to scan all following residues
                while (ii<atom_no){
                    ii=residue2_end+1;
                    j=residue2_end+1;
                    
                    //tries to find next residue
                    while (j<atom_no){
                        residue2_end=j;
                        if (residue!=atoms[j].residue && atoms[j].record=="ATOM  "){
                            residue=atoms[j].residue;
                            j=atom_no+1;                        
                        }
                        j++;
                    }
                    if (j==atom_no)
                        ii=atom_no;
                    else{
                    
                        //std::cout <<residue2_end<<std::endl;
                        
                        best_dist=max_dist+1;
                        std::cout<<std::to_string(atoms[i].residue);
                        std::cout<<";";
                        std::cout<<std::to_string(atoms[residue2_end].residue);
                        std::cout<<std::endl;
                        j=i;
                        }
                        
                        residue=atoms[j].residue;   
                        
                            
                }
            
            residue=atoms[i].residue;   
            }
        }
    }
    
    return 0;
}


int main()
{
    std::string temp_filename="";
    DIR *files;
    struct dirent *ent;
    std::string pdb_file="";
    std::string pdb_directory="D:/My documents/C++/chains/";
    
    //contact_map("D:/My documents/C++/chains/1FXL_A.pdb");
    //contact_map("D:/My documents/C++/chains/1B35_B.pdb");
    pdb_object("D:/My documents/C++/chains/1PID_A.pdb");
    return 0;
        
    files=opendir(pdb_directory.c_str());
    while ((ent = readdir (files)) != NULL) {
        pdb_file=ent->d_name;
        std::cout <<pdb_file;
        if (pdb_file.length()>4 && pdb_file.substr(pdb_file.length()-4,4)==".pdb"){
        
            //printf ("%s\n", ent->d_name);
        
            
        }
    }

    closedir(files);
    
    //contact_map("D:/My documents/C++/chains/1FXL_A.pdb");
}
