//reads a pdb file and returns the minimal distance between sidechain atoms for each residue
//residues: variable which stores the 1st and last atom number for each residue
//atoms: xyz coordinates of all atoms
//the program assumes clean PDBs, i.e. one chain, no dual residue definitions, etc.
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


    vector<point> atoms;

//opens pdb file
    std::string map_filename=pdb_filename.substr(0,pdb_filename.length()-4)+".txt";

    ifstream pdb_file;
    pdb_file.open(pdb_filename);


//reads through pdb file and write ATOM line into pdb vector
    while ( getline (pdb_file,line) ){
        if (line.substr(0,6)=="ATOM  "||line.substr(0,6)=="HETATM"){
            atom_no++;
            pdb.push_back(line);
        }
        else if (line.substr(0,7)=="ENDMDL ")
            break;
    }

    pdb_file.clear();
    pdb_file.close();

//finds first and last residue in pdb vector
    for (int i=0;i<atom_no;i++){
        if (stoi(pdb[i].substr(22,4))>residue){
            if (first_residue>last_residue){
                first_residue=stoi(pdb[i].substr(22,4));
            }
            residue=stoi(pdb[i].substr(22,4));
            if (residue>last_residue)
                last_residue=residue;
        }
    }

//defines residues, even numbers: atom number of first atom of residue, odd numbers: last atom of residue, refers to atoms
    int * residues;
    residues=new int[(last_residue-first_residue+1)*2];

    residue=-1;
    for (int i=0; i<atom_no; i++) {
        if (stoi(pdb[i].substr(22,4))!=residue){
            residue=stoi(pdb[i].substr(22,4));
            start_atom=stoi(pdb[i].substr(6,5))-1;
            residues[(residue-first_residue)*2]=start_atom;
        }
        else
            residues[(residue-first_residue)*2+1]=stoi(pdb[i].substr(6,5))-1;
    }


    for (int i=0; i<atom_no; i++){
        mypoint={stof(pdb[i].substr(30,8)),stof(pdb[i].substr(38,8)),stof(pdb[i].substr(46,8))};
        atoms.push_back (mypoint);
    }

ofstream map_file;
map_file.open(map_filename);

//calculates distance between every atom
    for (int i=0;i<last_residue-first_residue-2;i++){
        for (int j=1;j<last_residue-first_residue-1;j++){
            if (j>i){
                if (pdb[residues[i*2]].substr(0,6)=="ATOM  "){

                    dist=0;
                    best_dist=max_dist;
                    for (int x=0;x<residues[i*2+1]-residues[i*2];x++){
                        for (int y=0;y<residues[j*2+1]-residues[j*2];y++){
                            dist=pow(atoms[residues[i*2]+x].x-atoms[residues[j*2]+y].x,2);
                            dist=dist+pow(atoms[residues[i*2]+x].y-atoms[residues[j*2]+y].y,2);
                            dist=dist+pow(atoms[residues[i*2]+x].z-atoms[residues[j*2]+y].z,2);
                            dist=sqrt(dist);

                            if (dist<best_dist){
                                best_dist=dist;
                            }
                        }
                    }
                }
                if (best_dist<max_dist){
                    line=std::to_string(i+first_residue)+";"+std::to_string(j+first_residue)+";"+std::to_string(best_dist)+"\n";
                    map_file <<line;
                }
            }
        }
    }

    delete [] residues;

    atoms.clear();
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

    //contact_map("D:/My documents/C++/chains/1FXL_A.pdb");
    //contact_map("D:/My documents/C++/chains/1B35_B.pdb");



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
            }
        }
    }

    closedir(files);
    return 0;
    //contact_map("D:/My documents/C++/chains/1FXL_A.pdb");
}
