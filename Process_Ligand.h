#ifndef PROCESS_LIGAND_H_
#define PROCESS_LIGAND_H_

#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <vector>

#define MIN(a, b)  (((a) < (b)) ? (a) : (b))

#define E  2.718281
#define PI 3.141592
#define Err 0.050      // 5.0% error estimate 
#define NPATHS 100     // maximum number of paths per atom (shortest path algo.)
#define PLANAR_CRIT 0.987 // planar criteria for aromatic
#define HYBRIDATION_LEVEL 3   // level of criterias to impose a bonding - 1: based on dist.
                              //                                          2: based on ang.
                              //                                          3: based on dist. and ang.

#define MAX_PATH 250
#define MAX_MAP 100000
#define MAX_CYCLE 250
#define MAX_HETATM 500
#define MAX_BONDS 10
#define ATOM_ORDER_LEN 10
#define MAXCONECT 8

#define STEREO_THRESHOLD 0.50
#define DEF_SIGMA_ANG 13
#define DEF_SIGMA_DIS 10
#define PROP_ANG 0.50
#define PROP_DIS 0.50
#define SP3_ANG 109.5
#define SP2_ANG 120
#define SP1_ANG 180

#define GAUSSIAN(x,u,o)  pow(E,-pow((x-u),2)/(2*o))


typedef struct struct_tarjan tarjan;
typedef struct atom_struct atom;
typedef struct bond_length_struct bond_length;
typedef struct struct_bond bond;
typedef struct struct_res residue;
typedef struct struct_subgraph subgraph;

struct struct_subgraph {
	int size;         // size of graph (fragment)
	int remain;       // number of nodes remaining to build
	int recsize;      // recursive size

	int id;           // unique id.
	bool anchor;      // can serve as anchor
	int state;        // state: 0 undone 1 wait 2 done

	atom* at;         // atom that connects to root
	subgraph* root;   // from which the fragments branches from
	subgraph* prev;   // linked-list graphs
};

struct struct_res {
	char  name[4];
	char  chain;
	int   number;
};

struct struct_tarjan {
	int   index;
	int   lowlink;
};

struct struct_bond {
	atom*  to;
	atom*  from;

	int    type;
	int    cyclic;
	int    flexible;
	float  dist;
	//char   sym_state;
};

struct atom_struct {
	char  name[5];
	char  element[3];
	int   number;
	int   index;
	char  type[6];       // sybyl atom type

	int   atomtype;      // atom type (property of atom)
	int   nonmetal;      // atom is a non-metal
	int   aromatic;      // atom is apart aromatic cycle
	int   n_bonds;       // number of bonds
	int   shortest;      // number of atoms in shortest path to GPA

	bond* conect;        // bond structure
	subgraph* graph;     // subgraph atom belongs to
	float coor[3];
	float charge;        // coulomb charge
	
	int   gpa;           // gpa flag 1/0
	int   build_state;   // builtlist state   1: built
	                     //                   0: ready
	                     //                  -1: failed

	atom* buildlist[3];  // first 3 neighbours
	int   sp_state;      // state -   0: none  - 1: ready
                             //           2: check - 3: done
	int     sp_n_paths;    // number of paths
	atom*** sp_paths;      // shortest path
	int     sp_paths_n;    // number of connections in paths
	
	atom* next_build;     // next atom built in IC

	bool marked;         // used in breadth frsirst search
	int  sym_index;      // seen in %d-th passage
	tarjan vertex;       // tarjan algo to find cycles

	// internal coordinates to neighbours
	float dis;
	float ang;
	float dih;

	int cycle5;
	int cycle6;
};


static const char SYBYLbondtypes[11][20] = { "none", "single", "double", "triple", "quadruple", "aromatic", "polymeric", "delocalised double", "pi", "amide", "unknown" };

static const char Types_GAUDREAULT[13][20] = { "None", "Strong_Doneptor", "Strong_Acceptor", "Strong_Donor", "Weak_Doneptor", "Weak_Acceptor",
					       "Halogen", "Hydrophobic", "Aromatic", "Neutral", "Positive", "Negative", "Electrophilic" };


static const char Types_SOBOLEV[11][20] = { "None", "Hydrophilic", "Acceptor", "Donor",
					    "Hydrophobic", "Aromatic", "Neutral", 
					    "Neutral-Donor", "Neutral-Acceptor" };


void print_command_line();
void parse_command_line(int argv, char** argc, char* filename, char* outname,int* verbose,
			int* hydro_flex, int* remove_hydro,int* force_gpa, float** force_pcg,
			int* atom_index, residue* res, char* extract, int* reference, 
			int* old_types, int* new_types, int* babel_types, int* convert_only, 
			int* process_only,int* gen3D, char* outformat, int* target);

void set_OutBase(char* filename,char* outname, char* basepath, char* informat);
residue* get_Extract_List(char* extract_string,int* n_extract, residue* res,char* informat);
int is_Mapped(int number, int* map);
int is_Extractable(residue* extract, int n_extract, char (*args)[25]);
void Replace_Hyphens(char* string);

//atom* read_PDB(char* filename, int* n_atoms, int* map, residue* extract, int n_extract, float* pcg_ori);
atom* read_MOL2(char* filename, int* n_atoms, int* map, residue* extract, int n_extract, float* pcg_ori, int atom_index);
int Copy_OriginalMOL2(char* oldfilename, char* error);
int Convert_2_MOL2(char* filename, const char* informat, const char* outformat, char* error,int convert_only,int gen3D);
int read_Type(char* type);
int get_Format(char* filename, char* informat);
int get_Target_Filename(char* filename, char* target_filename);
atom* get_Atom_from_coor(const float* coor, atom* atoms, int n_atoms);

int in_stack(atom* atomw, int* st, int n_st);
void strongconnect(atom* atomv, atom* atomf, atom* atoms, int n_atoms,int* st, int* n_st,
		   int* n_root, int* scc, int* n_scc, int* index_t);
int Tarjan(atom* atoms, int n_atoms, int* scc, int* n_scc);
void save_Shortest_to_GPA(atom* atoms, int n_atoms);

char* UPPER(char* string);
void get_Element_From_Hybridation(char* string, char* dest);
int is_Cyclic(atom* atomz, int* scc, int n_scc);
int is_Planar(atom* atom1, atom* atom2, atom* atom3, atom* atom4, atom* atom5, atom* atom6);
int is_NonMetal(char* element);
int is_Hydrogen(atom* atomz);
int get_Shortest_Path(atom* atomi, atom* atoms, int n_atoms);

// functions used with the old atom types
int bonds_Hydrophilic(atom* atomzero);
int bonds_Donor(atom* atomzero);
int bonds_Acceptor(atom* atomzero);
int bonds_Carbon(atom* atomzero);
int bonds_Carbocation(atom* atomzero);
int bonds_sp2_Carbon(atom* atomzero);

void set_AtomTypes(atom* atoms, int n_atoms, int old_types, int new_types, int babel_types, int verbose);
void set_AtomTypes_SOBOLEV(atom* atomzero, int verbose);
void set_AtomTypes_GAUDREAULT(atom* atomzero, int verbose);
void set_AtomTypes_SYBYL(atom* atomzero, int verbose);
void atomtype_by_charge(atom* atomzero);
int set_Flexible_Bonds(atom* atoms, int n_atoms);
void set_Cyclic_Bonds(atom* atoms, int n_atoms,int *scc, int n_scc);
int Bond_Exists(bond* b, bond* blist[], int nb);
void print_bond_status(bond* conect, int status);
void get_Flexible_Atoms(atom* atoms,int n_atoms,bond* flex,atom* atomlist[],int* sense,int* nlist, int remove_hydro);

int Get_NextConnection(char* buffer);
void Print_Connections(atom* atoms, int n_atoms);
void free_Paths(atom** atoms, int n_atoms);
int memAlloc_Paths(atom** atoms, int n_atoms);
void Copy_Paths(atom* atomdest, atom* atomsrc);
void Print_Paths(atom* atoms, int n_atoms);
int count_Hydrogens(atom* atomzero);
int count_Carbons(atom* atomzero);
int count_Heavy(atom* atomzero);
int count_Oxygens(atom* atomzero);
int count_Nitrogens(atom* atomzero);
int count_Phosphorus(atom* atomzero);
int count_Aromatic(atom* atomzero);
int count_Cyclic(atom* atomzero);

int is_Methyl(atom* atomzero);
int is_Guanidium(bond* conect);
int is_Imine(bond* conect);
int is_Amide(bond* conect);
int is_Triple(bond* conect);
int is_Planar_Amine(bond* conect);
int is_Amine(atom* atomzero);
int is_Aromatic_Amidine(bond* conect);
int is_Aromatic_Amine(bond* conect);
int is_Carbon_Amine(atom* atomzero);
int is_Carbon_Carboxylate(atom* atomzero);
int is_Aromatic_Sulfonate(bond* conect);
int is_Between_Perpendicular_Aromatic(bond* conect);
int is_Sulfonate(atom* atomzero);
int is_Aromatic_Nitro(bond* conect);
int is_Nitro(atom* atomzero);
int is_Aromatic_Carboxylate(bond* conect);
int is_Terminal(bond* conect);
int is_Meta(bond* conect);
int is_MetaAmine(bond* conect);
int has_MetaGroup(atom* atomz);

void roll_Cycle(const atom* atom1, atom** atom2, atom** atom3, bool inverse);

void Print_subGraph(subgraph* graph, atom* atoms, int n_atoms);
int validate_Graph(subgraph* graph);
subgraph* subGraph_Molecule(atom* atoms, int n_atoms);
atom* get_NextGraphAtom(atom* atoms, int n_atoms);
subgraph* get_LargestGraph(subgraph* graph);
void reset_Graph(subgraph* graph);
void connect_Graph(subgraph* graph, subgraph* anchor_graph, atom* atoms, int n_atoms);
void set_Recursive_Size(subgraph* graph);
subgraph* get_ChildestGraph(subgraph* graph, int* recsize);
int anchor_Graph(subgraph* build_graph, subgraph* anchor_graph, atom* atoms, int n_atoms);
subgraph* get_BuildableGraph(subgraph* graph);
subgraph* get_NextUnconnectedGraph(subgraph* graph);
subgraph* get_NextAnchorGraph(subgraph* graph);
int is_RigidPath(atom** sp_path, int sp_path_n, atom* dest);
int is_Flexible(atom* from, atom* to);

atom* get_Force_gpa(atom* atoms, int n_atoms, int force_gpa);
atom* get_Free_gpa(atom* atoms, int n_atoms);
void get_Ligand_Center_Geometry(atom* atoms, int n_atoms, float* lig_ori);
atom* get_Center_Geometry_Atom(atom* atoms, int n_atoms, float* lig_ori);
int get_2_Consecutive_Bonded_Atom(atom* from, bond* conect1, bond* conect2);
int get_2_Consecutive_Bonded_Graph_Atom(atom* from, bond* conect1, bond* conect2);
atom* get_2_Common_Bonded_Atom(atom* atoms, int n_atoms, atom** gpa2, atom** gpa3);
atom* get_GPA_from_AnchorGraph(subgraph* graph, subgraph* anchor_graph, atom* atoms, int n_atoms, atom* gpa, atom** gpa2, atom** gpa3);
subgraph* get_GPAGraph(subgraph* graph, atom* gpa);
int is_AngleGraph(subgraph* anchor_graph, atom* atoms, int n_atoms);
atom* get_GPA_from_AngleGraph(subgraph* anchor_graph, atom* atoms, int n_atoms);
atom* get_Free_gpa(atom* atoms, int n_atoms);
atom* get_Buildable(atom* atoms, int n_atoms,subgraph* graph, atom* gpa);
atom* BuildList(atom* atoms, int n_atoms, atom* build, atom* sequence);
int reset_BuildableGraph(subgraph* graph);
int validate_Graphs(subgraph* graph);
int is_Built(atom* atomb);
void Reset_Buildable(atom* atoms, int n_atoms);
void Print_BuildList(atom* atomz);
atom* get_NonBL_Connection(atom* atomz, atom* build);
void Translate(atom* pcg, float* force_pcg);

void Write_Target(char* original_filename, char* target_filename, atom* atoms, int n_atoms);
void Write_IC(char* filename, atom* atoms, int n_atoms, atom* gpa, int remove_hydro);
void Write_INP(char* filename, char* icfile, atom* atoms, int n_atoms, int remove_hydro, residue* res, residue* force_outres,atom* gpa, subgraph* graph);
void Write_REF(char* filename, atom* atoms, int n_atoms, int remove_hydro, residue* res, residue* force_outres);
void buildcc(atom* sequence);

// Geometry file
float angle(const float *a, const float *b, const float *c);
float bndang(const float* a,const float* b, const float* c);
float dihedral(const float *a1,const float *a2,const float *a3,const float *a4);
float dihang(const float *a1,const float *a2,const float *a3,const float *a4);
float vec_norm(const float *v);
void vec_sub(float *a, const float *b, const float *c);
float dot_prod(const float *v1, const float *v2);
float dist(const float *a, const float *b);
float *cross_prod(float *x1, const float *x2, const float *x3);

#endif
