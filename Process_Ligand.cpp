#include "Process_Ligand.h"
#include "openbabel/obconversion.h"

#define HYDROPHOBIC_THRESHOLD       0.0655    // based on long aliphatic chains with a Csp3 (  CH3-R   )
                                              // see TEST folder in src directory

//#define NEUTRAL_THRESHOLD           0.1250

using namespace std;
//using namespace OpenBabel;

int main(int argv, char* argc[])
{
	int i,j;           // dumb counters

	atom* atoms=NULL;    // atoms structure
	atom* gpa=NULL;      // reference atom to anchor ligand
	atom* gpa2=NULL,*gpa3=NULL;
	atom* build=NULL;    // next atom to build in buildlist
	atom* sequence=NULL;

	subgraph* graph=NULL, *build_graph;
	subgraph* anchor_graph=NULL;

	atom  pcg1,pcg2,pcg3;      //  ligand center geometry atom used in buildlist (number is 0)


    
	float  lig_pcg[3];
	float* ori_pcg=NULL;       //  protein center geometry (origin)
	float* force_pcg=NULL;     //forces a protein center geometry

	int* map_atom;            // maps atom number to index in atoms struct.

	int scc[MAX_CYCLE];   // tarjan algo: strongly connected components
	int n_scc;

	int n_atoms=0;       // number of PDB atoms
	int n_branch=0;      // number of branching atoms when the anchor graph is angular

	int remove_hydro=1;  // removes hydrogens in output files
	int n_hydro=0;       // number of hydrogen atoms in ligand

	int hydro_flex=0;    // includes hydrogen flexible bonds
	int force_gpa=-1;    // forces a reference atom

	int old_types=0;
	int new_types=0;
	int sybyl_types=1;
	
	int target=0;        // process the target molecule
	int atom_index=0;    // starts atom indexing at
	int reference=0;     // output reference PDB file

	int convert_only=0;  // convert only PDB, do not process;
	int process_only=0;  // convert only PDB, do not process;
	int gen3D=0;         // generate a 3d conformation of the ligand

	int verbose=0;       // outputing detailed informations

	char original_filename[MAX_PATH]; // original input file
	char filename[MAX_PATH]; // PDB ligand file
	char target_filename[MAX_PATH]; // inp pdb target file
	char basepath[MAX_PATH]; // base path of filename
	char outname[MAX_PATH]; // output BASE filenames
	char icfile[MAX_PATH]; // ligand ic file (.ic)
	char inpfile[MAX_PATH]; // ligand input file (.inp)
	char reffile[MAX_PATH]; // reference PDB file
	char error[MAX_PATH];   // error messages

	residue outres;
	residue force_outres;

	char informat[6];
	char outformat[6] = "mol2";

	char extract_string[1000];
	residue* extract;
	int n_extract = 0;

	int n=0;

	outname[0] = '\0';
	extract_string[0] = '\0';

	strcpy(force_outres.name,"   ");
	force_outres.chain = '-';
	force_outres.number = -1;

	strcpy(outres.name,"LIG");
	outres.chain = ' ';
	outres.number = 9999;

	// set default protein center geometry
	pcg1.number = 0;        //invariable
	pcg2.number = 0;        //invariable
	pcg3.number = 0;        //invariable

	for(i=0; i<3; i++){
		pcg1.coor[i] = 0.0f;
		pcg2.coor[i] = 0.0f;
		pcg3.coor[i] = 0.0f;
	}

	pcg1.coor[0] = 1.0f;
	pcg3.coor[1] = 1.0f;

	parse_command_line(argv,argc,filename,outname,&verbose,
			   &hydro_flex,&remove_hydro,&force_gpa,&force_pcg,
			   &atom_index,&force_outres,extract_string,
			   &reference,&old_types,&new_types,&sybyl_types,
			   &convert_only,&process_only,&gen3D,outformat,&target);


	map_atom = (int*)malloc(MAX_MAP*sizeof(int));
	if(map_atom == NULL){
		return 1;
	}else{
		memset(map_atom,-1,MAX_MAP*sizeof(int));
	}

	extract = get_Extract_List(extract_string,&n_extract,&outres,informat);

	if(n_extract > 0){
		ori_pcg = (float*)malloc(3*sizeof(float));

		if(ori_pcg == NULL){
			printf("ERROR: could not allocate memory for ori_pcg\n");
			return(2);
		}

		memset(ori_pcg,0.0f,3*sizeof(float));
	}

	if(!process_only)
	{

		strcpy(original_filename,filename);

		if(get_Format(filename,informat)){
			fprintf(stderr,"ERROR: unknown format for file %s\n", filename);
			return(2);
		}
		
		printf("File format is '%s'\n", informat);

		/*
		  if(convert_only){
		  printf("Output file format is '%s'\n", outformat);
		  }
		*/
        
//        cout << "informat is :" << informat << endl;
		if(strcmp(informat, "mol2")) // if format is different from 'mol2'
		{
			if( !(n=Convert_2_MOL2(filename,informat,outformat,error,convert_only,gen3D)) )
			{
				fprintf(stderr,"%s",error);
				return(2);
			}
			else
			{
				printf("OpenBabel :: converted %d molecule%s\n",n, n>1?"s":"");
			}

			if(convert_only)
			{
				printf("Done.\n");
				return 0;
			}
		}
		else //
		{
			cout << "No conversion needed for the file to the mol2 format." << endl;
		}
	}

	atoms = read_MOL2(filename,&n_atoms,map_atom,extract,n_extract,ori_pcg,atom_index);
	set_AtomTypes(atoms, n_atoms, old_types, new_types, sybyl_types, verbose);
	
	if(target){
		get_Target_Filename(filename,target_filename);
		Write_Target(original_filename,target_filename,atoms,n_atoms);
		return(0);
	}

	// no atoms matches the extraction list OR no atoms in PDB/MOL2 file.
	if(atoms == NULL || n_atoms == 0){ printf("ERROR: no atoms found\n"); return 1;}

	printf("read %d atoms\n", n_atoms);
	get_Ligand_Center_Geometry(atoms,n_atoms,lig_pcg);
	

	if(force_pcg != NULL){
		// Adds up the user-defined PCG
		printf("translation of force_pcg (%8.3f %8.3f %8.3f)\n",
		       force_pcg[0], force_pcg[1], force_pcg[2]);

		Translate(&pcg1,force_pcg);
		Translate(&pcg2,force_pcg);
		Translate(&pcg3,force_pcg);

	}else if(ori_pcg != NULL){
		// One residue was extracted
		printf("translation of ori_pcg (%8.3f %8.3f %8.3f)\n",
		       ori_pcg[0], ori_pcg[1], ori_pcg[2]);

		Translate(&pcg1,ori_pcg);
		Translate(&pcg2,ori_pcg);
		Translate(&pcg3,ori_pcg);
	}

 	// Force a GPA Atom
	if(force_gpa > -1){
		gpa = get_Force_gpa(atoms,n_atoms,force_gpa);
	}


	// calculate the total number of hydrogen atoms
	n_hydro = 0;
	for(i=0; i<n_atoms; i++){
		if(is_Hydrogen(&atoms[i])){
			n_hydro++;
		}
	}

	printf("number of Hydrogens is %d\n", n_hydro);
	// 3 global positionning atoms are needed to build a ligand
	if(n_atoms - n_hydro < 3){ printf("ERROR: the ligand needs at least 3 heavy atoms\n"); return 1; }


	//Print_Connections(atoms,n_atoms);
	if(Tarjan(atoms,n_atoms,scc,&n_scc) != n_atoms){
		printf("molecule is cyclic\n");
	}else{
		printf("molecule is acyclic\n");
	}
	
	// set bonds as flexible that meet criteria
	printf("Found %d flexible bond(s)\n", set_Flexible_Bonds(atoms,n_atoms));

	
	// divide molecule into sub graphs
	graph = subGraph_Molecule(atoms,n_atoms);
	//Print_subGraph(graph,atoms,n_atoms);

	if(gpa == NULL){
		printf("largest graph will serve as anchor\n");
		anchor_graph = get_LargestGraph(graph);
	}else{
		printf("gpa forced to %d: gpa graph will serve as anchor\n", gpa->number);
		anchor_graph = get_GPAGraph(graph,gpa);
	}

	printf("gpa graph has size %d\n", anchor_graph->size);
	anchor_graph->root = anchor_graph;

	connect_Graph(graph, anchor_graph, atoms, n_atoms);
	reset_Graph(graph);
	set_Recursive_Size(graph);

	// GPA processing
	gpa = get_GPA_from_AnchorGraph(graph, anchor_graph, atoms, n_atoms, gpa, &gpa2, &gpa3);

	if(gpa == NULL){
		printf("ERROR: could not assign gpa atoms.\n");
		return 3;
	}
	else
	{

		if((n_branch=is_AngleGraph(anchor_graph, atoms, n_atoms)))
		{
			printf("gpa is in angle graph: will validate gpa atom\n");

			if((gpa->n_bonds - count_Hydrogens(gpa)) != 1)
			{
				gpa = get_GPA_from_AngleGraph(anchor_graph, atoms, n_atoms);
			}

			gpa->build_state = 1;
			gpa->buildlist[0] = &pcg1;
			gpa->buildlist[1] = &pcg2;
			gpa->buildlist[2] = &pcg3;

			// Buildlist sequence
			sequence = gpa;

			for(i=1; i<=gpa->n_bonds; i++)
			{
				if(!is_Hydrogen(gpa->conect[i].to))
				{
					gpa2 = gpa->conect[i].to;
					break;
				}
			}

			gpa2->build_state = 1;
			gpa2->buildlist[0] = gpa;
			gpa2->buildlist[1] = &pcg2;
			gpa2->buildlist[2] = &pcg3;

			sequence->next_build = gpa2;
			sequence = sequence->next_build;

			for(j=0; j<n_branch; j++){
				for(i=1; i<=gpa2->n_bonds; i++){
					if(!is_Hydrogen(gpa2->conect[i].to) &&
					   gpa2->conect[i].to != gpa){

						sequence->next_build = gpa2->conect[i].to;
						sequence = sequence->next_build;

						sequence->build_state = 1;
						sequence->buildlist[0] = gpa2;
						sequence->buildlist[1] = gpa;
						sequence->buildlist[2] = &pcg3;

						if(j == 0){ gpa3 = sequence; }
					}
				}
			}

		}else{
			bond conect1,conect2;

			if(!get_2_Consecutive_Bonded_Graph_Atom(gpa,&conect1,&conect2)){
				gpa2 = conect1.to;
				gpa3 = conect2.to;
			}else{
				if(!get_2_Consecutive_Bonded_Atom(gpa,&conect1,&conect2)){
					gpa2 = conect1.to;
					gpa3 = conect2.to;
				}else{
					printf("ERROR: could not retrieve 3 consecutive atoms in graph\n");
					return 5;
				}
			}

			gpa->build_state = 1;
			gpa->buildlist[0] = &pcg1;
			gpa->buildlist[1] = &pcg2;
			gpa->buildlist[2] = &pcg3;

			// Buildlist sequence
			sequence = gpa;

			gpa2->build_state = 1;
			gpa2->buildlist[0] = gpa;
			gpa2->buildlist[1] = &pcg2;
			gpa2->buildlist[2] = &pcg3;

			sequence->next_build = gpa2;
			sequence = sequence->next_build;

			gpa3->build_state = 1;
			gpa3->buildlist[0] = gpa2;
			gpa3->buildlist[1] = gpa;
			gpa3->buildlist[2] = &pcg3;

			sequence->next_build = gpa3;
			sequence = sequence->next_build;

		}

	}


	gpa->gpa = 1;
	printf("reference gpa atom gpa is %d\n", gpa->number);

	get_Shortest_Path(gpa,atoms,n_atoms);
	save_Shortest_to_GPA(atoms,n_atoms);
	//Print_Paths(atoms,n_atoms);
	free_Paths(&atoms,n_atoms);

	// build atoms sequence
	reset_Graph(graph);

	build_graph = anchor_graph;

	gpa->graph->remain--;
	gpa2->graph->remain--;
	gpa3->graph->remain--;
	
	while((build=get_Buildable(atoms,n_atoms,build_graph,gpa)) != NULL)
	{
		get_Shortest_Path(build,atoms,n_atoms);
		//Print_Paths(atoms,n_atoms);
		sequence=BuildList(atoms,n_atoms,build,sequence);
		free_Paths(&atoms,n_atoms);
	}
		
	// Base name in output files
	set_OutBase(filename,outname,basepath,informat);
	printf("will output to prefix '%s'\n", outname);
	//printf("basepath is '%s'\n", basepath);

	// Write internal coordinates of neighbouring atoms
	strcpy(icfile,outname);
	strcat(icfile,".ic");
	Write_IC(icfile, atoms, n_atoms, gpa, remove_hydro);

	// Write neighbours and atom types and flexible bonds
	strcpy(inpfile,outname);
	strcat(inpfile,".inp");
	Write_INP(inpfile, icfile, atoms, n_atoms, remove_hydro, &outres, &force_outres, gpa, graph);

	// write reference PDB file
	if(reference){
		strcpy(reffile,outname);
		strcat(reffile,"_ref");
		strcat(reffile,".pdb");

		buildcc(gpa);
		Write_REF(reffile, atoms, n_atoms, remove_hydro, &outres, &force_outres);
	}

	// memory dealloc's
	if(force_pcg != NULL) { free(force_pcg); }
	if(ori_pcg != NULL) { free(ori_pcg); }
	if(extract != NULL) { free(extract); }

	if(graph != NULL){
		build_graph = graph->prev;
		free(graph);
		graph = build_graph;
	}

	// free up atoms structure in memory
	for(i=0; i<n_atoms; ++i){
		if(atoms[i].conect != NULL){
			free(atoms[i].conect);
		}
	}
	free(atoms);
	free(map_atom);


	printf("Done.\n");

	return 0;
}

void set_AtomTypes(atom* atoms, int n_atoms, int old_types, int new_types, int sybyl_types, int verbose){

	if(old_types){
		// all atoms except Carbons
		for(int j=0; j<n_atoms; ++j){
			if(strncmp(atoms[j].type,"C.",2)){
				set_AtomTypes_SOBOLEV(&atoms[j],verbose);
			}
		}

		// Carbons only
		for(int j=0; j<n_atoms; ++j){
			if(!strncmp(atoms[j].type,"C.",2)){
				set_AtomTypes_SOBOLEV(&atoms[j],verbose);
			}
		}

	}else if(sybyl_types){
		for(int j=0; j<n_atoms; ++j){
			set_AtomTypes_SYBYL(&atoms[j],verbose);
		}

	}else{ // if(new_types){
		for(int j=0; j<n_atoms; ++j){
			set_AtomTypes_GAUDREAULT(&atoms[j],verbose);
		}

	}
}

int Copy_OriginalMOL2(char* oldfilename, char* error){

	char newfilename[255];

	strcpy(newfilename,oldfilename);
	strcat(oldfilename,".old");

	//cout << "will copy to " << oldfilename << endl;
	ofstream oldfile(oldfilename);
	ifstream newfile(newfilename);
	string line;

	if(!oldfile || !newfile){
		strcpy(error,"ERROR: could not open filename to copy original MOL2 file\n");
		return 1;
	}

	while (newfile.good()){
		getline(newfile, line);
		oldfile << line << endl;
	}

	newfile.close();
	oldfile.close();

	return 0;

}

int get_Target_Filename(char* filename, char* target_filename)
{
	char* pch;

	strcpy(target_filename,filename);

	pch = strstr(target_filename,".mol2.tmp");
	if(pch != NULL){
		strcpy(pch,".inp.pdb");
	}else{
		return 1;
	}

	return 0;
}

int Convert_2_MOL2(char* filename, const char* informat, const char* outformat, char* error,int convert_only,int gen3D){

	ifstream ifs_exist;
	char suffix[10];
	//int i=2;

	//cout << "will open " << filename << " for input" << endl;

	ifstream ifs(filename);
	if(!ifs){
		strcpy(error,"ERROR: could not open input file for Convert\n");
		return 0;
	}

	char outfilename[MAX_PATH];
	strcpy(outfilename,filename);

	char* pch = strrchr(outfilename,'.');
	strcpy(suffix,".mol2.tmp");
	strcpy(pch,suffix);

	/*
	  strcpy(suffix,"_tmp.");
	  if(convert_only){
	  strcat(suffix,outformat);
	  }else{
	  strcat(suffix,"mol2");
	  }
	  strcpy(pch,suffix);
	*/

	//cout << "will open " << outfilename << " for output" << endl;

	ofstream ofs(outfilename);
	if(!ofs){
		strcpy(error,"ERROR: could not open output file for Convert\n");
		return 0;
	}

	OpenBabel::OBConversion conv(&ifs,&ofs);

	if(!conv.SetInAndOutFormats(informat,outformat)){
		strcpy(error,"ERROR: OpenBabel: Formats not available\n");
		return 0;
	}

	// default error msg
	strcpy(error,"ERROR: no molecules converted\n");

	// add Hydrogens
//	conv.AddOption("h",OpenBabel::OBConversion::GENOPTIONS);
	
	// physiological pH
	conv.AddOption("p",OpenBabel::OBConversion::GENOPTIONS,"7.0");
	
	// generate 3D conf.
	if(gen3D){ conv.AddOption("gen3D",OpenBabel::OBConversion::GENOPTIONS); }

	int n = conv.Convert();
	
	ofs.close();
	ifs.close();
	
	strcpy(filename,outfilename);
	
	return(n);
}


int get_Format(char* filename, char* informat){

	char* pch = strrchr(filename,'.');

	if(pch != NULL)
		strcpy(informat,++pch);
	else return 1;

	return 0;
}

char* UPPER(char* string){

	static char STRING[MAX_PATH]={'\0'};
	strcpy(STRING,string);

	// Upper-str function
	for(int i=0; i<(int)strlen(STRING); i++)
		STRING[i] = toupper(STRING[i]);

	return STRING;
}

residue* get_Extract_List(char* extract_string,int* n_extract, residue* outres, char* informat){

	residue* res_list=NULL;
	residue* r;
	char* pch;
	char temp[12];
	int len;

	*n_extract = 0;

	if(!strcmp(extract_string,"")) { return res_list; }

	pch = strtok (extract_string,",");
	while (pch != NULL)
	{
		//printf ("%s\n",pch);
		(*n_extract)++;

		if(*n_extract == 1){
			res_list = (residue*)malloc(*n_extract*sizeof(residue));
			if(res_list == NULL){
				fprintf(stderr,"could not allocate memory for res_list\n");
				exit(2);
			}
		}else{
			res_list = (residue*)realloc(res_list,*n_extract*sizeof(residue));
			if(res_list == NULL){
				fprintf(stderr,"could not reallocate memory for res_list\n");
				exit(2);
			}
		}

		r = &res_list[*n_extract-1];

		if(!strcmp(UPPER(informat),"MOL2")){
			//file is MOL2
			r->number = atoi(pch);
			printf("will extract residue number %d\n", r->number);

			if((*n_extract-1) == 0)
				outres->number = r->number;
			else outres->number = 9999;

		}else if(!strcmp(UPPER(informat),"PDB")){
			//file is PDB
			strncpy(r->name,pch,3);
			r->name[3]='\0';

			r->chain = pch[strlen(pch)-1];

			len = strlen(pch)-4;
			strncpy(temp,&pch[3],len);
			temp[len]='\0';

			r->number = atoi(temp);

			Replace_Hyphens(r->name);
			if(r->chain == '-'){ r->chain = ' '; }

			printf("will extract %s %d %c\n", r->name, r->number, r->chain);

			if((*n_extract-1) == 0){
				outres->chain = r->chain;
				outres->number = r->number;
				strcpy(outres->name,r->name);
			}else{
				outres->chain = ' ';
				outres->number = 9999;
				strcpy(outres->name,"LIG");
			}
		}


		pch = strtok (NULL, ",");
	}

	return res_list;

}

void Replace_Hyphens(char* string){

	int i;

	for(i=0; i<(int)strlen(string); ++i){
		string[i] = string[i] == '-' ? ' ' : string[i];
	}

}

void set_OutBase(char* filename, char* outname, char* basepath,char* informat){

	char* pch, *pch2;

	// default outbase is the filename with .ext removed
	// e.g. ../BTN.mol2 -> ../BTN

	//printf("outname: '%s'\n", outname);
	if(!strcmp(outname,"")) {

		if((pch=strstr(filename,".mol2")) != NULL){

#ifdef _WIN32
			if((pch2=strrchr(filename,'\\')) != NULL && pch2 < pch){
				strncpy(basepath,filename,pch2-filename+1);
				basepath[pch2-filename+1]='\0';
			}else{
				strcpy(basepath,"");
			}

#else
			if((pch2=strrchr(filename,'/')) != NULL && pch2 < pch) {
				strncpy(basepath,filename,pch2-filename+1);
				basepath[pch2-filename+1]='\0';
			}else{
				strcpy(basepath,"");
			}
#endif

			strncpy(outname,filename,pch-filename);
			outname[pch-filename]='\0';

		}else{
			strcpy(outname,filename);
			strcpy(basepath,"");
		}
	}

}


void Print_BuildList(atom* atomz){

	int i;

	printf("buildlist of %6d: ", atomz->number);
	for(i=0; i<3; i++){
		printf("%6d", atomz->buildlist[i]->number);
	}
	printf("\n");

}

void Translate(atom* pcg, float* translate_pcg){

	pcg->coor[0] += translate_pcg[0];
	pcg->coor[1] += translate_pcg[1];
	pcg->coor[2] += translate_pcg[2];

}

subgraph* subGraph_Molecule(atom* atoms, int n_atoms){

	subgraph* newgraph, *graph = NULL;
	atom* root = NULL;

	while((root=get_NextGraphAtom(atoms,n_atoms)) != NULL){

		//printf("using number %d as graph root\n", root->number);

		newgraph = (subgraph*)malloc(sizeof(subgraph));
		if(newgraph == NULL){
			printf("memory allocation error for newgraph\n");
			return NULL;
		}

		root->graph = newgraph;
		newgraph->size = 1;
		newgraph->at = NULL;
		newgraph->root = NULL;
		newgraph->anchor = true;
		newgraph->state = 0;

		if(graph == NULL){
			newgraph->id = 1;
		}else{
			newgraph->id = graph->id + 1;
		}

		get_Shortest_Path(root,atoms,n_atoms);

		for(int i=0; i<n_atoms; i++)
		{
			if(root != &atoms[i] &&
			   atoms[i].graph == NULL &&
			   /*!is_Hydrogen(&atoms[i]) &&*/
			   is_RigidPath(atoms[i].sp_paths[0],atoms[i].sp_paths_n,&atoms[i])){
//				printf("ADDED %d TO NODE LIST\n", atoms[i].number);
				atoms[i].graph = newgraph;
				newgraph->size++;
			}
			//getchar();
		}

		free_Paths(&atoms,n_atoms);

		newgraph->remain = newgraph->size;

		// graphs linked list
		newgraph->prev = graph;
		graph = newgraph;
	}

	return graph;
}

int is_RigidPath(atom** sp_path, int sp_path_n, atom* dest){

	//printf("sp_path_n = %d\n", sp_path_n);
	for(int i=0; i<sp_path_n; i++){
		if(i == (sp_path_n-1)){
			if(is_Flexible(sp_path[i],dest) == 1){
				return 0;
			}
		}else{
			if(is_Flexible(sp_path[i],sp_path[i+1]) == 1){
				return 0;
			}
		}
	}

	return 1;

}

int is_Flexible(atom* from, atom* to){

	//printf("analyzing from atom[%d] to atom[%d]\n", from->number, to->number);
	for(int i=1; i<=from->n_bonds; i++){
		if(from->conect[i].to == to){
			//printf("atoms are connected!\n");
			if(from->conect[i].flexible){
				// atoms are connect and bond is flexible
				return 1;
			}

			// atoms are connected but bond is rigid
			return 0;
		}
	}

	// atoms are not connected
	return -1;
}

atom* get_NextGraphAtom(atom* atoms, int n_atoms){

	for(int i=0; i<n_atoms; i++){
		if(!is_Hydrogen(&atoms[i]) &&
		   atoms[i].graph == NULL) {
			return &atoms[i];
		}
	}

	return NULL;

}

subgraph* get_NextAnchorGraph(subgraph* graph){

	while(graph != NULL){
		if(graph->root != NULL &&
		   graph->anchor){
			return graph;
		}
		graph = graph->prev;
	}

	return NULL;
}

subgraph* get_NextUnconnectedGraph(subgraph* graph){

	while(graph != NULL){
		if(graph->root == NULL &&
		   !graph->state){
			return graph;
		}

		graph = graph->prev;
	}

	return NULL;
}

subgraph* get_GPAGraph(subgraph* graph, atom* gpa){

	while(graph != NULL){
		if(gpa->graph == graph){
			return graph;
		}

		graph = graph->prev;
	}

	return NULL;
}

subgraph* get_LargestGraph(subgraph* graph){

	int max = 0;
	subgraph* maxgraph = NULL;

	while(graph != NULL){
		if(graph->size > max){
			max = graph->size;
			maxgraph = graph;
		}

		graph = graph->prev;
	}

	return maxgraph;
}

void set_Recursive_Size(subgraph* graph){

	int recsize = 0;
	subgraph* child_graph;

	while((child_graph=get_ChildestGraph(graph,&recsize)) != NULL){
		child_graph->recsize = recsize + child_graph->size;
		child_graph->state = 2;
	}

}

subgraph* get_ChildestGraph(subgraph* graph, int* recsize){

	subgraph* root, *default_graph = NULL;
	subgraph* graph1 = graph,* graph2 = graph;

	// retrieve the graph none graph roots to and not already flagged as done
	while(graph1 != NULL){

		if(graph1->state == 2){ graph1 = graph1->prev; continue; }

		root = NULL;
		*recsize = 0;

		while(graph2 != NULL){

			// is the graph rooted from another graph
			if(graph2->root == graph1){
				(*recsize) += graph2->recsize;

				/*
				  printf("graph[%d] roots graph[%d] (%s)\n", graph1->id, graph2->id,
				  graph2->done?"DONE":"NOT DONE");
				*/

				if(!graph2->state){ root = graph1; }
			}

			graph2 = graph2->prev;
		}

		if(root == NULL){ return graph1; }

		default_graph = graph1;

		graph1 = graph1->prev;
		graph2 = graph;

	}

	return default_graph;
}

void connect_Graph(subgraph* graph, subgraph* anchor_graph, atom* atoms, int n_atoms){

	subgraph* build_graph;

	do{
		//printf("anchor set to graph[%d]\n", anchor_graph->id);

		while((build_graph = get_NextUnconnectedGraph(graph)) != NULL){

			//printf("next graph to be built is graph[%d]\n", build_graph->id);

			anchor_Graph(build_graph, anchor_graph, atoms, n_atoms);
			build_graph->state = 2;
		}

		anchor_graph->anchor = false;
		reset_Graph(graph);

	}while((anchor_graph=get_NextAnchorGraph(graph)) != NULL);

}

void reset_Graph(subgraph* graph){

	while(graph != NULL){
		graph->state = 0;
		graph = graph->prev;
	}

}

int anchor_Graph(subgraph* build_graph, subgraph* anchor_graph, atom* atoms, int n_atoms){

	//printf("will try to anchor graph[%d] onto graph[%d]\n", build_graph->id, anchor_graph->id);

	for(int i=0; i<n_atoms; i++){
		if(atoms[i].graph == build_graph){
			//printf("for atom[%d] in graph[%d]\n", atoms[i].number, atoms[i].graph->id);

			for(int j=0; j<n_atoms; j++){
				if(atoms[j].graph == anchor_graph){
					/*
					  printf("check connection with atom[%d] from graph[%d]\n",
					  atoms[j].number, atoms[j].graph->id);
					*/

					if(is_Flexible(&atoms[i],&atoms[j]) == 1){
						/*
						  printf("graph[%d] anchored onto graph[%d]\n",
						  build_graph->id, anchor_graph->id);
						*/

						build_graph->root = anchor_graph;

						return 1;
					}
				}
			}
		}
	}

	return 0;
}

void Reset_Buildable(atom* atoms, int n_atoms){

	int i;

	for(i=0; i<n_atoms; i++){
		if(atoms[i].build_state == -1) { atoms[i].build_state = 0; }
	}

}

int validate_Graphs(subgraph* graph)
{
	while(graph != NULL){
		if(graph->state != 2){
			return(1);
		}
		graph = graph->prev;
	}
	return(0);
}

int reset_BuildableGraph(subgraph* graph)
{
	static int i = -1;
	int j = 0;

	// resets wait graphs to undone
	while(graph != NULL){
		if(graph->state == 1){
			graph->state = 0;
			j++;
		}
		graph = graph->prev;
	}

	if(i != -1 && i == j)
		return(1);

	i = j;
	return 0;
}

subgraph* get_BuildableGraph(subgraph* graph){

	subgraph* graph_start = graph;
	
	while(graph != NULL){
		if(!graph->state){
			//graph->root->state == 2){
			return graph;
		}
		graph = graph->prev;
	}

	// no undone graphs, are there any wait graphs
	if(reset_BuildableGraph(graph_start))
		return NULL;

	graph = graph_start;
	while(graph != NULL){
		if(!graph->state){
			//graph->root->state == 2){
			return graph;
		}
		graph = graph->prev;
	}

	return NULL;
}

atom* get_Buildable(atom* atoms, int n_atoms, subgraph* build_graph, atom* gpa){

	atom* mincon_atom = NULL;
	int mincon = 1000;

	for(int i=0; i<n_atoms; i++){
		if(atoms[i].build_state == 0){
			if(atoms[i].shortest < mincon){
				mincon = atoms[i].shortest;
				mincon_atom = &atoms[i];
			}
		}
	}
	
	return mincon_atom;
}

int is_Built(atom* atomb){

	//if atom is not built (0) or failed (-1)
	//the returned value is false

	/*
	  printf("is %d built: %s\n", atomb->number,
	  atomb->build_state == 1 ? "YES" : "NO");
	*/

	if (atomb == NULL)
		return 0;

	return (atomb->build_state == 1 ? 1 : 0);

}

void save_Shortest_to_GPA(atom* atoms, int n_atoms){

	for(int i=0; i<n_atoms; i++){
		atoms[i].shortest = atoms[i].sp_paths_n;
		//printf("atoms[%d].shortest = %d\n", atoms[i].number, atoms[i].shortest);
	}

}

atom* BuildList(atom* atoms, int n_atoms, atom* build, atom* sequence)
{
	int i,j,k,n;
	int path_fail=0;
	atom* atomc = NULL;
	
	//printf("==================\n");
	//printf("in BuildList for %d\n", build->number);

	for(i=0; i<n_atoms; i++)
	{
		if(build->build_state == 1){break;}
		
		if(atoms[i].sp_paths_n == 0 || is_Hydrogen(&atoms[i])){continue;}

		//printf("----------------\n");
		//printf("from atom %d\n", atoms[i].number);
		
		n=atoms[i].sp_paths_n >=4 ? 4 : atoms[i].sp_paths_n;
		
		// loop through all paths of atom[i] to atom 'build'
		for(j=0; j<atoms[i].sp_n_paths; j++)
		{
			path_fail = 0;

			build->buildlist[0] = NULL;
			build->buildlist[1] = NULL;
			build->buildlist[2] = NULL;
			
			for(k=1; k<n; k++)
			{
				if(!is_Built(atoms[i].sp_paths[j][k]))
				{
					path_fail = 1;
					break;
				}
				else
				{
					build->buildlist[k-1] = atoms[i].sp_paths[j][k];
				}
			}

			
			if(n < 4 && !path_fail){
				if(!is_Built(&atoms[i])){
					path_fail = 1;
				}else{
					build->buildlist[n-1] = &atoms[i];

					if((n + 1) < 4){
						if(!is_Built(atomc=get_NonBL_Connection(&atoms[i],build))){
							path_fail = 1;
						}else{
							build->buildlist[n] = atomc;
							
							if((n + 2) < 4){
								if(!is_Built(atomc=get_NonBL_Connection(atomc,build))){
									path_fail = 1;
								}else{
									build->buildlist[n+1] = atomc;
								}
							}
						}
					}
				}
			}

			if(path_fail)
			{
//				printf("%d failed\n", atoms[i].number);
				continue;
			}
			else
			{
				build->build_state = 1;
				break;
			}

		}
	}

	if(build->build_state == 1)
	{
		//Print_Paths(atoms,n_atoms);
		sequence->next_build = build;
		sequence = sequence->next_build;

		Reset_Buildable(atoms,n_atoms);
		//Print_BuildList(build);

		build->graph->remain--;
	}
	else
	{
//		printf("%d FLAGGED AS FAILED\n",build->number);
		build->build_state = -1;
	}

	return sequence;

}

atom* get_NonBL_Connection(atom* atomz, atom* build){

	int i,j,in;

	for(i=1; i<=atomz->n_bonds; i++){
		if(atomz->conect[i].to == build ||
		   is_Hydrogen(atomz->conect[i].to))
		{
			continue;
		}


		in=0;
		for(j=0; j<3; j++){
			if(build->buildlist[j] != NULL &&
			   atomz->conect[i].to == build->buildlist[j])
			{
				in=1;
				break;
			}
		}

		if(!in)
			return atomz->conect[i].to;
	}

	return NULL;
}

atom* get_Free_gpa(atom* atoms, int n_atoms){

	int i;

	for(i=0; i<n_atoms; i++){
		if(is_Hydrogen(&atoms[i])) { continue; }

		if(atoms[i].gpa == 0){
			atoms[i].gpa = 1;
			return &atoms[i];
		}
	}

	return NULL;
}

void get_Ligand_Center_Geometry(atom* atoms, int n_atoms, float* lig_ori){

	int i,tot=0;

	lig_ori[0] = 0.0f;
	lig_ori[1] = 0.0f;
	lig_ori[2] = 0.0f;

	for(i=0; i<n_atoms; i++){
		if(is_Hydrogen(&atoms[i])) { continue; }

		lig_ori[0] += atoms[i].coor[0];
		lig_ori[1] += atoms[i].coor[1];
		lig_ori[2] += atoms[i].coor[2];

		tot++;
	}

	lig_ori[0] /= (float)tot;
	lig_ori[1] /= (float)tot;
	lig_ori[2] /= (float)tot;

}

atom* get_Center_Geometry_Atom(atom* atoms, int n_atoms, float* lig_ori){

	int i;
	float d,dmin=9999.999f;
	atom* atommin=NULL;

	for(i=0; i<n_atoms; i++){
		if(is_Hydrogen(&atoms[i])) { continue; }

		if((d=dist(lig_ori,atoms[i].coor)) < dmin){
			dmin = d;
			atommin = &atoms[i];
		}
	}

	return atommin;
}

atom* get_GPA_from_AngleGraph(subgraph* anchor_graph, atom* atoms, int n_atoms){

	for(int i=0; i<n_atoms; i++){
		if(atoms[i].graph == anchor_graph){
			int heavy_conect = atoms[i].n_bonds - count_Hydrogens(&atoms[i]);
			if(heavy_conect == 1){
				return &atoms[i];
			}
		}
	}

	return NULL;

}


int is_AngleGraph(subgraph* anchor_graph, atom* atoms, int n_atoms){

	int n_branch = 0, n_root = 0;

	for(int i=0; i<n_atoms; i++){
		if(atoms[i].graph == anchor_graph){

			int heavy_conect = atoms[i].n_bonds - count_Hydrogens(&atoms[i]);
			if(heavy_conect == 1){
				n_branch++;
			}else if(heavy_conect > 1){
				n_root++;
			}
		}
	}

	return (n_root == 1 && n_branch > 1) ? n_branch : 0;
}

atom* get_GPA_from_AnchorGraph(subgraph* graph, subgraph* anchor_graph, atom* atoms, int n_atoms,
			       atom* gpa,atom** gpa2, atom** gpa3){

	int recsize = 0;
	subgraph *lgraph = NULL, *_graph = graph;

	if(gpa == NULL){

		// find the 2nd largest graph
		while(_graph != NULL){
			if(_graph->root == anchor_graph &&
			   _graph != anchor_graph &&
			   _graph->recsize > recsize){
				recsize = _graph->recsize;
				lgraph = _graph;
			}

			_graph = _graph->prev;
		}

		if(lgraph == NULL){
			for(int i=0; i<n_atoms; i++){
				if(atoms[i].graph == anchor_graph){
					return &atoms[i];
				}
			}
		}


		//printf("2nd largest recursive graph is %d\n", lgraph->id);

		// find the atom that connects to anchor_graph from 2nd largest graph
		for(int i=0; i<n_atoms; i++){
			if(atoms[i].graph == lgraph){
				for(int j=1; j<=atoms[i].n_bonds; j++){
					if(atoms[i].conect[j].to->graph == anchor_graph &&
					   atoms[i].conect[j].flexible){
						gpa = atoms[i].conect[j].to;
						break;
					}
				}

				if(gpa != NULL){ break; }
			}
		}
	}

	return gpa;

}

int get_2_Consecutive_Bonded_Graph_Atom(atom* from, bond* conect1, bond* conect2){

	for(int i=1; i<=from->n_bonds; i++){
		atom* to = from->conect[i].to;

		if(to->graph == from->graph &&
		   !is_Hydrogen(to)){

			*conect1 = from->conect[i];

			for(int j=1; j<=to->n_bonds; j++){
				atom* to2 = to->conect[j].to;

				if(to2->graph == from->graph &&
				   !is_Hydrogen(to2) &&
				   to2 != from){

					*conect2 = to->conect[j];

					return 0;
				}
			}
		}
	}

	return 1;
}

int get_2_Consecutive_Bonded_Atom(atom* from, bond* conect1, bond* conect2){

	for(int i=1; i<=from->n_bonds; i++){
		atom* to = from->conect[i].to;

		if(!is_Hydrogen(to)){

			*conect1 = from->conect[i];

			for(int j=1; j<=to->n_bonds; j++){
				atom* to2 = to->conect[j].to;

				if(!is_Hydrogen(to2) &&
				   to2 != from){

					*conect2 = to->conect[j];

					return 0;
				}
			}
		}
	}

	return 1;
}


atom* get_2_Common_Bonded_Atom(atom* atoms, int n_atoms, atom** gpa2, atom** gpa3){

	int i,j;
	int k,l;
	atom* _atom = NULL, *__atom = NULL;

	for(i=0; i<n_atoms; i++){

		_atom = &atoms[i];
		if(is_Hydrogen(_atom)) { continue; }

		if((_atom->n_bonds-count_Hydrogens(_atom)) != 1){ continue; }

		for(j=0; j<n_atoms; j++){
			if(i == j){ continue; }

			__atom = &atoms[j];
			if(is_Hydrogen(__atom)) { continue; }

			if((__atom->n_bonds-count_Hydrogens(__atom)) != 1){ continue; }
			(*gpa3) = __atom;

			for(k=1; k<=_atom->n_bonds; k++){
				if(is_Hydrogen(_atom->conect[k].to)){ continue; }

				for(l=1; l<=__atom->n_bonds; l++){
					if(is_Hydrogen(__atom->conect[l].to)){ continue; }

					if(_atom->conect[k].to == __atom->conect[l].to &&
					   _atom->conect[k].to->n_bonds == 3){

						(*gpa2) = _atom->conect[k].to;
						return _atom;
					}
				}
			}
		}
	}

	(*gpa2) = NULL;
	(*gpa3) = NULL;
	_atom = NULL;
	__atom = NULL;

	return _atom;
}

atom* get_Force_gpa(atom* atoms, int n_atoms, int force_gpa){

	int i;

	for(i=0; i<n_atoms; i++){
		if(atoms[i].number == force_gpa &&
		   !is_Hydrogen(&atoms[i])) {
			return &atoms[i];
		}
	}

	return NULL;
}

void print_bond_status(bond* conect, int status){

	if(status == 1){
		printf("bond[%d][%d] is Terminal\n", conect->to->number, conect->from->number);
	}else if(status == 2){
		printf("bond[%d][%d] is Imine\n", conect->to->number, conect->from->number);
	}else if(status == 3){
		printf("bond[%d][%d] is Double_or_Amide\n", conect->to->number, conect->from->number);
	}else if(status == 4){
		printf("bond[%d][%d] is Triple\n", conect->to->number, conect->from->number);
	}else if(status == 5){
		printf("bond[%d][%d] is Planar_Amine\n", conect->to->number, conect->from->number);
	}else if(status == 6){
		printf("bond[%d][%d] is Aromatic_Amidine\n", conect->to->number, conect->from->number);
	}else if(status == 7){
		printf("bond[%d][%d] is Aromatic_Sulfonate\n", conect->to->number, conect->from->number);
	}else if(status == 8){
		printf("bond[%d][%d] is Aromatic_Nitro\n", conect->to->number, conect->from->number);
	}else if(status == 9){
		printf("bond[%d][%d] is Aromatic_Carboxylate\n", conect->to->number, conect->from->number);
	}else if(status == 10){
		printf("bond[%d][%d] is Aromatic_Amide\n", conect->to->number, conect->from->number);
	}
	
}

int set_Flexible_Bonds(atom* atoms,int n_atoms){

	int i,j;
	int flex_counter = 0;
	atom* atom1_ptr, *atom2_ptr;

	for(i=0; i<n_atoms; i++)
	{
		atom1_ptr=&atoms[i];

		if(is_Hydrogen(atom1_ptr)) { continue; }
		else if(is_Methyl(atom1_ptr)) { continue; }

		for(j=1; j<=atom1_ptr->n_bonds; j++)
		{
			atom2_ptr = atom1_ptr->conect[j].to;
			if(is_Hydrogen(atom2_ptr)) { continue; }
			else if(atom1_ptr->conect[j].type != 1 || atom1_ptr->conect[j].cyclic) { continue; }

			else if(is_Terminal(&atom1_ptr->conect[j])){ print_bond_status(&atom1_ptr->conect[j],1); continue; }
			else if(is_Imine(&atom1_ptr->conect[j])){ print_bond_status(&atom1_ptr->conect[j],2); continue; }
			else if(is_Amide(&atom1_ptr->conect[j])){ print_bond_status(&atom1_ptr->conect[j],3); continue; }
			else if(is_Triple(&atom1_ptr->conect[j])){ print_bond_status(&atom1_ptr->conect[j],4); continue; }
			else if(is_Planar_Amine(&atom1_ptr->conect[j])){ print_bond_status(&atom1_ptr->conect[j],5); continue; }
			else if(is_Aromatic_Amidine(&atom1_ptr->conect[j])){ print_bond_status(&atom1_ptr->conect[j],6); continue; }
			else if(is_Aromatic_Sulfonate(&atom1_ptr->conect[j])){ print_bond_status(&atom1_ptr->conect[j],7); continue; }
			else if(is_Aromatic_Nitro(&atom1_ptr->conect[j])){ print_bond_status(&atom1_ptr->conect[j],8); continue; }
			else if(is_Aromatic_Carboxylate(&atom1_ptr->conect[j])){ print_bond_status(&atom1_ptr->conect[j],9); continue; }
			else if(is_Aromatic_Amide(&atom1_ptr->conect[j])){ print_bond_status(&atom1_ptr->conect[j],10); continue; }
			
			atom1_ptr->conect[j].flexible = 1;
			flex_counter++;
		}
	}

	return flex_counter/2;
}

void set_Cyclic_Bonds(atom* atoms,int n_atoms, int* scc,int n_scc){

	int i,j;
	atom* atomz1,*atomz2;

	for(i=0; i<n_atoms; i++)
	{
		atomz1=&atoms[i];

		if(is_Cyclic(atomz1,scc,n_scc))
		{
			for(j=1; j<=atomz1->n_bonds; j++)
			{
				atomz2=atomz1->conect[j].to;

				if(is_Cyclic(atomz2,scc,n_scc))
				{
					atomz1->conect[j].cyclic = 1;
					//printf("bond[%d][%d] is cyclic\n", atomz1->number, atomz2->number);
				}
			}
		}
	}

}


int in_stack(atom* atomw, int* st, int n_st){

	int i;

	for(i=0; i<n_st; i++){
		if(st[i] == atomw->number){
			return 1;
		}
	}

	return 0;
}

void strongconnect(atom* atomv, atom* atomf, atom* atoms, int n_atoms, int* st, int* n_st,
		   int* n_root, int* scc, int* n_scc, int* index_t){

	atom* atomw=NULL;
	int i,t;

	atomv->vertex.index = *index_t;
	atomv->vertex.lowlink = *index_t;
	(*index_t)++;

	st[(*n_st)++] = atomv->number;

	//printf("strongconnect: strongconnect of atom %d\n", atomv->number);

	for(i=1; i<=atomv->n_bonds; i++){
		atomw = atomv->conect[i].to;

		if(atomf != NULL && atomw->number == atomf->number) { continue; }

		if(atomw->vertex.index == -1){
			strongconnect(atomv->conect[i].to,atomv,atoms,n_atoms,st,n_st,n_root,scc,n_scc,index_t);
			atomv->vertex.lowlink = MIN(atomv->vertex.lowlink,atomw->vertex.lowlink);
		}else if(in_stack(atomv->conect[i].to,st,*n_st)){
			atomv->vertex.lowlink = MIN(atomv->vertex.lowlink,atomw->vertex.index);
		}
	}

	if(atomv->vertex.lowlink == atomv->vertex.index){

		(*n_root)++;
		*n_scc = 0;

		t=0;

		do{
			scc[(*n_scc)++] = st[--(*n_st)];
			t++;
		}while(scc[*n_scc-1] != atomv->number);

		// remove individual nodes
		if(t == 1)
		{
			(*n_scc)--;
		}else{
			/*
			  printf("Tarjan(%d): ", *n_scc);
			  for(i=0; i<*n_scc; i++){
			  printf("%6d",scc[i]);
			  }
			  printf("\n");
			*/

			set_Cyclic_Bonds(atoms,n_atoms,scc,*n_scc);

		}

	}

}

int Tarjan(atom* atoms, int n_atoms, int* scc, int* n_scc){

	int i;

	int st[MAX_HETATM];
	int n_st=0;
	int n_root=0; // returns the number of root nodes found

	int index_t=0;

	for(i=0; i<n_atoms; i++){
		if(atoms[i].vertex.index == -1){
			//printf("Tarjan: strongconnect of atom %d\n", atoms[i].number);
			strongconnect(&atoms[i],NULL,atoms,n_atoms,st,&n_st,&n_root,scc,n_scc,&index_t);
		}
	}

	return n_root;
}

int count_Hydrogens(atom* atomzero){

	int i;
	int tot=0;

	for(i=1; i<=atomzero->n_bonds; i++){
		if(!strcmp(atomzero->conect[i].to->type,"H"))
			tot++;
	}

	return tot;
}

int count_Heavy(atom* atomzero){

	int i;
	int tot=0;

	for(i=1; i<=atomzero->n_bonds; i++){
		if(strcmp(atomzero->conect[i].to->type,"H"))
			tot++;
	}

	return tot;
}

int count_Oxygens(atom* atomzero){

	int i;
	int tot=0;

	for(i=1; i<=atomzero->n_bonds; i++){
		if(!strncmp(atomzero->conect[i].to->type,"O.",2))
			tot++;
	}

	return tot;
}

int count_Nitrogens(atom* atomzero){

	int i;
	int tot=0;

	for(i=1; i<=atomzero->n_bonds; i++){
		if(!strncmp(atomzero->conect[i].to->type,"N.",2))
			tot++;
	}

	return tot;
}

int is_Methyl(atom* atomzero){

	int n_hydro = 0;

	if(!strncmp(atomzero->type,"C.3",3)){

		for(int i=1; i<=atomzero->n_bonds; i++){
			if(is_Hydrogen(atomzero->conect[i].to)){
				n_hydro++;
			}
		}

		if(n_hydro > 2) { return 1; }
	}

	return 0;

}

int is_Imine(bond* conect){

	if(!strncmp(conect->from->type,"N.2",3) &&
	   !strncmp(conect->to->type,"N.2",3)){
		return 1;
	}

	return 0;
}

int is_Amide(bond* conect){

	if(!strncmp(conect->to->type,"C.2",3) &&
	   (!strncmp(conect->from->type,"N.2",3) ||
	    !strncmp(conect->from->type,"N.AM",4) ||
	    !strncmp(conect->from->type,"N.PL3",5))) {
		return 1;

	}else if(!strncmp(conect->from->type,"C.2",3) &&
		 (!strncmp(conect->to->type,"N.2",3) ||
		  !strncmp(conect->to->type,"N.AM",4) ||
		  !strncmp(conect->to->type,"N.PL3",5))) {
		return 1;
	}

	return 0;
}

int is_Triple(bond* conect){

	if(!strncmp(conect->to->type,"C.1",3) &&
	   !strncmp(conect->from->type,"C.1",3)){
		return 1;
	}

	return 0;
}

int is_Planar_Amine(bond* conect){

	if(!strncmp(conect->to->type,"C.",2) &&
	   is_Amine(conect->from)){
		return 1;

	}else if(!strncmp(conect->from->type,"C.",2) &&
		 is_Amine(conect->to)){
		return 1;
	}

	return 0;
}

int is_Amine(atom* atomzero){

	int n_hydro = 0;

	if(!strncmp(atomzero->type,"N.2",3) || !strncmp(atomzero->type,"N.PL3",5)){
		for(int i=1; i<=atomzero->n_bonds; i++){
			if(is_Hydrogen(atomzero->conect[i].to)){
				n_hydro++;
			}
		}

		if(n_hydro == 2) { return 1; }
	}

	return 0;
}

int is_Aromatic_Amidine(bond* conect){

	if(!strncmp(conect->from->type,"C.AR",4) &&
	   is_Carbon_Amine(conect->to)){
		return 1;

	}else if(!strncmp(conect->to->type,"C.AR",4) &&
		 is_Carbon_Amine(conect->from)){

		return 1;
	}

	return 0;

}

int is_Aromatic_Amide(bond* conect){

	if(!strncmp(conect->from->type,"C.AR",4) &&
	   !strncmp(conect->to->type,"N.AM",4)){
		return 1;

	}else if(!strncmp(conect->to->type,"C.AR",4) &&
		 !strncmp(conect->from->type,"N.AM",4)){

		return 1;
	}
	
	return 0;

}

int is_Aromatic_Sulfonate(bond* conect){

	if(!strncmp(conect->from->type,"C.AR",4) &&
	   is_Sulfonate(conect->to)){
		return 1;

	}else if(!strncmp(conect->to->type,"C.AR",4) &&
		 is_Sulfonate(conect->from)){

		return 1;
	}

	return 0;

}

int is_Sulfonate(atom* atomzero){

	int n_oco2 = 0;
	if(!strncmp(atomzero->type,"S.3",3)){
		for(int i=1; i<=atomzero->n_bonds; i++){
			if(!strncmp(atomzero->conect[i].to->type,"O.CO2",5)){
				n_oco2++;
			}
		}

		if(n_oco2 == 3){ return 1; }
	}

	return 0;
}

int is_Aromatic_Nitro(bond* conect){

	if(!strncmp(conect->from->type,"C.AR",4) &&
	   is_Nitro(conect->to)){
		return 1;

	}else if(!strncmp(conect->to->type,"C.AR",4) &&
		 is_Nitro(conect->from)){

		return 1;
	}

	return 0;

}

int is_Nitro(atom* atomzero){

	if(!strncmp(atomzero->type,"N.",2)){
		for(int i=1; i<=atomzero->n_bonds; i++){
			if(!strncmp(atomzero->conect[i].to->type,"O.2",3)){
				return 1;
			}
		}
	}

	return 0;
}

int is_Aromatic_Carboxylate(bond* conect){

	if(!strncmp(conect->from->type,"C.AR",4) &&
	   is_Carbon_Carboxylate(conect->to)){
		return 1;

	}else if(!strncmp(conect->to->type,"C.AR",4) &&
		 is_Carbon_Carboxylate(conect->from)){

		return 1;
	}

	return 0;

}

int is_Carbon_Amine(atom* atomzero){

	int n_amine = 0;

	if(!strncmp(atomzero->type,"C.CAT",5) || !strncmp(atomzero->type,"C.2",3)){
		for(int i=1; i<=atomzero->n_bonds; i++){
			if(!strncmp(atomzero->conect[i].to->type,"N.PL3",5) ||
			   !strncmp(atomzero->conect[i].to->type,"N.2",3)){
				n_amine++;
			}
		}

		if(n_amine == 2) { return 1; }
	}

	return 0;
}

int is_Carbon_Carboxylate(atom* atomzero){

	int n_oco2 = 0;

	if(!strncmp(atomzero->type,"C.CAT",5) || !strncmp(atomzero->type,"C.2",2)){
		for(int i=1; i<=atomzero->n_bonds; i++){
			if(!strncmp(atomzero->conect[i].to->type,"O.CO2",5)){
				n_oco2++;
			}
		}

		if(n_oco2 == 2) { return 1; }
	}

	return 0;
}

int is_Guanidium(atom* atomzero){

	int i,j,dbl,sng;

	for(i=1; i<=atomzero->n_bonds; i++){
		if(!strncmp(atomzero->conect[i].to->type,"C.",2)){

			sng = 0;
			dbl = 0;
			if(count_Nitrogens(atomzero->conect[i].to) == 3 &&
			   count_Heavy(atomzero->conect[i].to) == 3){

				for(j=1; j<=atomzero->conect[i].to->n_bonds; j++){
					if(count_Heavy(atomzero->conect[i].to->conect[j].to) == 2){
						dbl = 1;

					}else if(count_Heavy(atomzero->conect[i].to->conect[j].to) == 1){
						sng++;
					}
				}

			}

			if(dbl == 1 && sng == 2) { return 1; }
		}
	}

	return 0;
}

int count_Carbons(atom* atomzero){

	int i;
	int tot=0;

	for(i=1; i<=atomzero->n_bonds; i++){
		if(!strncmp(atomzero->conect[i].to->type,"C.",2))
			tot++;
	}

	return tot;
}

int count_Phosphorus(atom* atomzero){

	int i;
	int tot=0;

	for(i=1; i<=atomzero->n_bonds; i++){
		if(!strncmp(atomzero->conect[i].to->type,"P.",2))
			tot++;
	}

	return tot;
}

int count_Cyclic(atom* atomzero){

	int i;
	int tot=0;

	for(i=1; i<=atomzero->n_bonds; i++){
		if(atomzero->conect[i].cyclic)
			tot++;
	}

	return tot;
}

int count_Aromatic(atom* atomzero){

	int i;
	int tot=0;

	for(i=1; i<=atomzero->n_bonds; i++){
		if(!strncmp(atomzero->conect[i].to->type,"C.AR",4))
			tot++;
	}

	return tot;
}

int bonds_Carbon(atom* atomzero){

	int nc=0;
	for(int i=1; i<=atomzero->n_bonds; i++){
		if(!strncmp(atomzero->conect[i].to->type,"C.",2)){ nc++; }
	}

	return nc;
}

int bonds_Carbocation(atom* atomzero){

	int ncat=0;
	for(int i=1; i<=atomzero->n_bonds; i++){
		if(!strncmp(atomzero->conect[i].to->type,"C.CAT",5)){ ncat++; }
	}

	return ncat;
}

int bonds_sp2_Carbon(atom* atomzero){

	int nc=0;
	for(int i=1; i<=atomzero->n_bonds; i++){
		if(!strncmp(atomzero->conect[i].to->type,"C.2",3)){ nc++; }
	}

	return nc;
}

int bonds_Donor(atom* atomzero){

	int nd=0;
	for(int i=1; i<=atomzero->n_bonds; i++){
		if(atomzero->conect[i].to->atomtype == 3){ nd++; }
	}

	return nd;
}

int bonds_Acceptor(atom* atomzero){

	int na=0;
	for(int i=1; i<=atomzero->n_bonds; i++){
		if(atomzero->conect[i].to->atomtype == 2){ na++; }
	}

	return na;
}

int bonds_Hydrophilic(atom* atomzero){

	int nh=0;
	for(int i=1; i<=atomzero->n_bonds; i++){
		if(atomzero->conect[i].to->atomtype == 1 ||
		   atomzero->conect[i].to->atomtype == 2 ||
		   atomzero->conect[i].to->atomtype == 3 ){ nh++; }
	}

	return nh;
}


void set_AtomTypes_SOBOLEV(atom* atomzero, int verbose){

	// Carbon
	if(!strncmp(atomzero->type,"C.",2)){
		if(!strncmp(atomzero->type,"C.AR",4)){
			atomzero->atomtype = 5;

		}else{
			if(!bonds_Hydrophilic(atomzero)){
				atomzero->atomtype = 4;
			}else if(bonds_Donor(atomzero) == 1){
				atomzero->atomtype = 7;
			}else if(bonds_Acceptor(atomzero) == 1){
				atomzero->atomtype = 8;
			}else{
				atomzero->atomtype = 6;
			}
		}

		//Oxygen
	}else if(!strncmp(atomzero->type,"O.",2)){

		if(!strncmp(atomzero->type,"O.3",3)){
			// Hydroxyl group
			if(count_Hydrogens(atomzero) > 0){
				atomzero->atomtype = 1;
			}else{
				atomzero->atomtype = 2;
			}
		}else if(!strncmp(atomzero->type,"O.2",3)){
			// Carbonyl group
			atomzero->atomtype = 2;
		}else if(!strncmp(atomzero->type,"O.CO2",5)){
			// Carboxyl group
			atomzero->atomtype = 2;
		}

		// Nitrogen
	}else if(!strncmp(atomzero->type,"N.",2)){

		if(bonds_Carbon(atomzero) == 3){
			atomzero->atomtype = 6;
		}else{

			if(!strncmp(atomzero->type,"N.4",3)){
				atomzero->atomtype = 3;
			}else if(!strncmp(atomzero->type,"N.PL3",5)){
				atomzero->atomtype = 3;
			}else if(!strncmp(atomzero->type,"N.AM",4)){
				// Amide group (resonance)
				if(count_Hydrogens(atomzero) > 0){
					atomzero->atomtype = 3;
				}else{
					atomzero->atomtype = 2;
				}
			}else if(!strncmp(atomzero->type,"N.AR",4)){
				if(count_Hydrogens(atomzero) > 0){
					atomzero->atomtype = 3;
				}else{
					atomzero->atomtype = 2;
				}
			}else{
				if(count_Hydrogens(atomzero) > 0){
					atomzero->atomtype = 1;
				}else{
					atomzero->atomtype = 2;
				}
			}
		}

		// Sulphur
	}else if(!strncmp(atomzero->type,"S.",2)){
		atomzero->atomtype = 6;
	}else if(!strncmp(atomzero->type,"F",1)){ // FLuoride
		atomzero->atomtype = 6;


	}else if(!strncmp(atomzero->type,"CL",2)){ // Chloride
		atomzero->atomtype = 4;
	}else if(!strncmp(atomzero->type,"BR",2)){ // Bromide
		atomzero->atomtype = 4;
	}else if(!strncmp(atomzero->type,"I",1)){ // Iodine
		atomzero->atomtype = 4;


	}else{
		// neutral (default)
		atomzero->atomtype = 6;
	}

	if(verbose >= 1){
		printf("atom(%s).charge=%6.3f\t.type(%s)=%s\n",
		       atomzero->name, atomzero->charge,
		       atomzero->type,
		       Types_SOBOLEV[atomzero->atomtype]);
	}

}

void set_AtomTypes_SYBYL(atom* atomzero, int verbose){

	if(!strncmp(atomzero->type,"C.1",3)){
		atomzero->atomtype = 1;
	}else if(!strncmp(atomzero->type,"C.2",3)){
		atomzero->atomtype = 2;
	}else if(!strncmp(atomzero->type,"C.3",3)){
		atomzero->atomtype = 3;
	}else if(!strncmp(atomzero->type,"C.AR",4)){
		atomzero->atomtype = 4;
	}else if(!strncmp(atomzero->type,"C.CAT",5)){
		atomzero->atomtype = 5;
	}else if(!strncmp(atomzero->type,"N.1",3)){
		atomzero->atomtype = 6;
	}else if(!strncmp(atomzero->type,"N.2",3)){
		atomzero->atomtype = 7;
	}else if(!strncmp(atomzero->type,"N.3",3)){
		atomzero->atomtype = 8;
	}else if(!strncmp(atomzero->type,"N.4",3)){
		atomzero->atomtype = 9;
	}else if(!strncmp(atomzero->type,"N.AR",4)){
		atomzero->atomtype = 10;
	}else if(!strncmp(atomzero->type,"N.AM",4)){
		atomzero->atomtype = 11;
	}else if(!strncmp(atomzero->type,"N.PL3",5)){
		atomzero->atomtype = 12;
	}else if(!strncmp(atomzero->type,"O.2",3)){
		atomzero->atomtype = 13;
	}else if(!strncmp(atomzero->type,"O.3",3)){
		// important distinction to make
		/*
		  if(count_Hydrogens(atomzero) == 1){
		  atomzero->atomtype = 14; // alcohol (O.AL)
		  }else{
		  atomzero->atomtype = 15; // ether (O.ET)
		  }
		*/
		atomzero->atomtype = 14;
	}else if(!strncmp(atomzero->type,"O.CO2",5)){
		//atomzero->atomtype = 16;
		atomzero->atomtype = 15;
	}else if(!strncmp(atomzero->type,"O.AR",4)){
		//atomzero->atomtype = 17;
		atomzero->atomtype = 16;
	}else if(!strncmp(atomzero->type,"S.2",3)){
		//atomzero->atomtype = 18;
		atomzero->atomtype = 17;
	}else if(!strncmp(atomzero->type,"S.3",3)){
		// important distinction to make (alcohol vs ether)
		/*
		  if(count_Hydrogens(atomzero) == 1){
		  atomzero->atomtype = 19; // alcohol
		  }else{
		  atomzero->atomtype = 20; // ether
		  }
		*/
		atomzero->atomtype = 18; // alcohol
	}else if(!strncmp(atomzero->type,"S.O2",4)){
		atomzero->atomtype = 22;
		atomzero->atomtype = 19;
	}else if(!strncmp(atomzero->type,"S.O",3)){
		//atomzero->atomtype = 21;
		atomzero->atomtype = 20;
	}else if(!strncmp(atomzero->type,"S.AR",4)){
		//atomzero->atomtype = 23;
		atomzero->atomtype = 21;
	}else if(!strncmp(atomzero->type,"P.3",3)){
		//atomzero->atomtype = 24;
		atomzero->atomtype = 22;
	}else if(!strncmp(atomzero->type,"CL",2)){
		//atomzero->atomtype = 26;
		atomzero->atomtype = 24;
	}else if(!strncmp(atomzero->type,"BR",2)){
		//atomzero->atomtype = 27;
		atomzero->atomtype = 25;
	}else if(!strncmp(atomzero->type,"SE",2)){
		//atomzero->atomtype = 29;
		atomzero->atomtype = 27;
	}else if(!strncmp(atomzero->type,"MG",2)){
		//atomzero->atomtype = 30;
		atomzero->atomtype = 28;
	}else if(!strncmp(atomzero->type,"SR",2)){
		//atomzero->atomtype = 31;
		atomzero->atomtype = 29;
	}else if(!strncmp(atomzero->type,"CU",2)){
		//atomzero->atomtype = 32;
		atomzero->atomtype = 30;
	}else if(!strncmp(atomzero->type,"MN",2)){
		//atomzero->atomtype = 33;
		atomzero->atomtype = 31;
	}else if(!strncmp(atomzero->type,"HG",2)){
		//atomzero->atomtype = 34;
		atomzero->atomtype = 32;
	}else if(!strncmp(atomzero->type,"CD",2)){
		//atomzero->atomtype = 35;
		atomzero->atomtype = 33;
	}else if(!strncmp(atomzero->type,"NI",2)){
		//atomzero->atomtype = 36;
		atomzero->atomtype = 34;
	}else if(!strncmp(atomzero->type,"ZN",2)){
		//atomzero->atomtype = 37;
		atomzero->atomtype = 35;
	}else if(!strncmp(atomzero->type,"CA",2)){
		//atomzero->atomtype = 38;
		atomzero->atomtype = 36;
	}else if(!strncmp(atomzero->type,"FE",2)){
		//atomzero->atomtype = 39;
		atomzero->atomtype = 37;
	}else if(!strncmp(atomzero->type,"F",1)){
		//atomzero->atomtype = 25;
		atomzero->atomtype = 23;
	}else if(!strncmp(atomzero->type,"I",1)){
		//atomzero->atomtype = 28;
		atomzero->atomtype = 26;
	}else if(!strncmp(atomzero->type,"CO.OH",5)){
		//atomzero->atomtype = 40;
		atomzero->atomtype = 38;
	}else{
		// dummy atom type
		//atomzero->atomtype = 41;
		atomzero->atomtype = 39;
	}
}

void set_AtomTypes_GAUDREAULT(atom* atomzero, int verbose){

	/*
	  I                  II                 III             IV               V
	  "Strong_Doneptor", "Strong_Acceptor", "Strong_Donor", "Weak_Doneptor", "Weak_Acceptor",
	  IV          VII          VIII        IX          X          XI             XII
	  "Halogen", "Hydrophobic", "Aromatic", "Neutral", "Positive", "Negative", "Electrophilic"
	*/

	// Carbon
	if(!strncmp(atomzero->type,"C.",2)){
		if(!strncmp(atomzero->type,"C.AR",4)){
			atomzero->atomtype = 8;

		}else if(!strncmp(atomzero->type,"C.CAT",5)){
			// guanadinium cation group
			// positively charged carbon
			atomzero->atomtype = 9;
		}else{
			atomtype_by_charge(atomzero);
		}

		//Oxygen
	}else if(!strncmp(atomzero->type,"O.",2)){

		// oxygen from phosphate
		if(count_Heavy(atomzero) == 1 && count_Phosphorus(atomzero) == 1){
			atomzero->atomtype = 1;

		}else if(count_Cyclic(atomzero) == 2 && count_Aromatic(atomzero) == 2){
			//printf("%d O is aromatic\n", atomzero->number);

			// in aromatic cycle
			atomtype_by_charge(atomzero);
			if(atomzero->atomtype == 9){
				atomzero->atomtype = 12;
			}else{
				atomzero->atomtype = 8;
			}

		}else if(count_Aromatic(atomzero) == 2){
			//printf("%d O is weak arom-O-arom\n", atomzero->number);

			// in between 2 aromatic cycles
			atomzero->atomtype = 5;

		}else if(count_Aromatic(atomzero) == 1 && count_Hydrogens(atomzero) == 0 &&
			 count_Heavy(atomzero) == 2){
			//printf("%d O is weak arom-O-aliph\n", atomzero->number);

			// in between an aromatic and aliphatic (ether)
			atomzero->atomtype = 5;

		}else if(!strncmp(atomzero->type,"O.2",3)){
			// Carbonyl group
			atomzero->atomtype = 2;

		}else if(!strncmp(atomzero->type,"O.CO2",5)){
			// Carboxyl group
			atomzero->atomtype = 11;

		}else if(bonds_sp2_Carbon(atomzero)){
			//printf("%d O is weak ester\n", atomzero->number);

			// binds a double bonded carbon
			// e.g. sp3 oxygen of an ester
			atomzero->atomtype = 5;

		}else if(!strncmp(atomzero->type,"O.3",3)){
			if(count_Hydrogens(atomzero) > 0){
				// Hydroxyl group
				atomzero->atomtype = 1;
			}else{
				// ether - aliph-O-aliph
				atomzero->atomtype = 2;
			}

		}else{
			printf("Unknown type for O %d\n", atomzero->number);
			atomzero->atomtype = 9;
		}

		// Nitrogen
	}else if(!strncmp(atomzero->type,"N.",2)){

		if(!strncmp(atomzero->type,"N.4",3)){
			// charged amide (ammomium primary)
			atomzero->atomtype = 10;

		}else if(!strncmp(atomzero->type,"N.PL3",5)){

			if(is_Guanidium(atomzero)){ //bonds_Carbocation(atomzero)){
				atomzero->atomtype = 10;

			}else if(count_Hydrogens(atomzero) > 0){
				atomzero->atomtype = 3;

			}else{
				// bound to 3 heavy atoms
				atomzero->atomtype = 12;
			}

		}else if(!strncmp(atomzero->type,"N.AM",4)){
			// Amide group (resonance)
			if(count_Hydrogens(atomzero) > 0){
				atomzero->atomtype = 3;
			}else{
				// Amide with 3 heavy atoms
				atomtype_by_charge(atomzero);
				if(atomzero->atomtype == 9){
					atomzero->atomtype = 12;
				}

			}

		}else if(!strncmp(atomzero->type,"N.AR",4)){
			if(count_Hydrogens(atomzero) > 0){
				atomzero->atomtype = 3;
			}else{

				if(count_Cyclic(atomzero) == 2 && count_Aromatic(atomzero) == 2 &&
				   count_Heavy(atomzero) == 2){
					//printf("%d N is aromatic\n", atomzero->number);
					atomzero->atomtype = 2;

				}else{
					atomtype_by_charge(atomzero);
					if(atomzero->atomtype == 9){
						atomzero->atomtype = 12;
					}else{
						atomzero->atomtype = 8;
					}
				}
			}
		}else{
			if(count_Hydrogens(atomzero) > 0){
				atomzero->atomtype = 1;
			}else{
				atomzero->atomtype = 2;
			}
		}

		// Sulphur
	}else if(!strncmp(atomzero->type,"S.",2)){


		if(count_Heavy(atomzero) == 2 && count_Cyclic(atomzero) == 2 &&
		   count_Aromatic(atomzero) == 2){

			// aromatic sulfur
			// if hydrophobic, assign as aromatic
			atomtype_by_charge(atomzero);
			if(atomzero->atomtype != 9){
				atomzero->atomtype = 8;
			}

		}else if(!strncmp(atomzero->type,"S.O2",4)){
			atomtype_by_charge(atomzero);

		}else if(!strncmp(atomzero->type,"S.O",3)){
			atomtype_by_charge(atomzero);

		}else{
			if(atomzero->n_bonds >= 4){
				atomtype_by_charge(atomzero);

			}else{
				if(count_Hydrogens(atomzero) > 0){
					atomzero->atomtype = 4;
				}else{
					atomzero->atomtype = 5;
				}
			}
		}

	}else if(!strncmp(atomzero->type,"P.",2)){

                // phophate group ?
		if(count_Oxygens(atomzero) < 3){
			// No, P.3 not in phosphate
			atomzero->atomtype = 5;

		}else{
			atomtype_by_charge(atomzero);
		}

	}else if(!strncmp(atomzero->type,"BR",2) ||
		 !strncmp(atomzero->type,"CL",2) ||
		 !strncmp(atomzero->type,"F",1) ||
		 !strncmp(atomzero->type,"I",1)) {

		// Halogen
		atomzero->atomtype = 6;

		//atomtype_by_charge(atomzero);
		/*
		  if(atomzero->atomtype == 9){
		  atomzero->atomtype = 12;
		  }
		*/
	}else{

		atomtype_by_charge(atomzero);
	}


	if(verbose >= 1){
		printf("atom(%s).charge=%6.3f\t.type(%s)=%s\n",
		       atomzero->name, atomzero->charge,
		       atomzero->type,
		       Types_GAUDREAULT[atomzero->atomtype]);
	}
}

void atomtype_by_charge(atom* atomzero){

	// neutral if over
	// otherwise, hydrophobic
	atomzero->atomtype = fabs(atomzero->charge) > HYDROPHOBIC_THRESHOLD ? 9 : 7;

}

void print_command_line(){

	printf("\n");
	printf("Obligatory argument(s):\n");
	printf("\t%-50s%-50s\n", "-f", "sets the ligand file");
	printf("\n");

	printf("Facultative argument(s):\n");
	printf("\t%-50s%-50s\n", "-v <INT>", "sets verbose level");
	printf("\t%-50s%-50s\n", "-o <STR>", "sets the output base filename");
	printf("\t%-50s%-50s\n", "-e <STR>", "sets a residue to extract in the file (only works for .mol2/.pdb files)");
	printf("\t%-50s%-50s\n", "-c <STR>", "converts a molecule only without processing it to a specified format\n");
	printf("\t%-50s%-50s\n", "-hf", "includes hydrogen flexible bonds");
	printf("\t%-50s%-50s\n", "-wh", "adds hydrogen atoms in output files\n");
	printf("\t%-50s%-50s\n", "-ref", "outputs the final PDB file from the IC\n");
	printf("\t%-50s%-50s\n", "-target", "parses a target (append atomtypes in columns[76-77] - PDB format only)\n");
	printf("\t%-50s%-50s\n", "--atom_index <INT>", "starts atom indexing at");
	printf("\t%-50s%-50s\n", "--res_name <STR>", "sets the 3-chars code for the ligand");
	printf("\t%-50s%-50s\n", "--res_chain <CHAR>", "sets the ligand chain");
	printf("\t%-50s%-50s\n", "--res_number <INT>", "sets the ligand number\n");
	printf("\t%-50s%-50s\n", "--force_gpa <INT>", "forces a reference atom (not yet working when converting is necessary)");
	printf("\t%-50s%-50s\n", "--force_pcg <FLOAT FLOAT FLOAT>", "forces a protein center of geometry\n");
	printf("\t%-50s%-50s\n", "--gen3D", "generates a 3D conformation of the ligand\n");
	//printf("\t%-50s%-50s\n", "--old_types", "uses the old atom types (I-VIII)");
	//printf("\t%-50s%-50s\n", "--new_types", "uses the new atom types (I-XII)");
	//printf("\t%-50s%-50s\n", "--sybyl_types", "uses the sybyl atom types (I-XXVI)\n");
	printf("\t%-50s%-50s\n", "--help", "prints this help menu");
	printf("\n");

}

void parse_command_line(int argv, char** argc, char* filename, char* outname,int* verbose,
			int* hydro_flex, int* remove_hydro, int* force_gpa, float** force_pcg,
			int* atom_index, residue* force_outres, char* extract_string, int* reference,
			int* old_types, int* new_types, int* sybyl_types, int* convert_only,
			int* process_only, int* gen3D, char* outformat, int* target){

	int i;

	if(argv == 1){
		print_command_line();
		exit(0);
	}

	for(i=0; i<argv; ++i)
	{
		if(!strcmp(argc[i],"-f")){
			strcpy(filename,argc[++i]);
		}else if(!strcmp(argc[i],"-o")){
			strcpy(outname,argc[++i]);
		}else if(!strcmp(argc[i],"-e")){
			strcpy(extract_string,argc[++i]);
		}else if(!strcmp(argc[i],"-c")){
			*convert_only=1;
			//strcpy(outformat,argc[++i]);
		}else if(!strcmp(argc[i],"-p")){
			*process_only=1;
		}else if(!strcmp(argc[i],"--gen3D")){
			*gen3D=1;
		}else if(!strcmp(argc[i],"-v")){
			*verbose=atoi(argc[++i]);
			//*verbose=1;
		}else if(!strcmp(argc[i],"--help") || !strcmp(argc[i],"-h")){
			print_command_line();
			exit(0);
		}
		else if(!strcmp(argc[i],"-hf")){
			*hydro_flex = 1;
			*remove_hydro = 0;
	  	}
	  	else if(!strcmp(argc[i],"-wh")){
			*remove_hydro = 0;
			*hydro_flex = 1;
		}else if(!strcmp(argc[i],"-ref")){
			*reference = 1;
		}else if(!strcmp(argc[i],"-target")){
			*target = 1;
		}else if(!strcmp(argc[i],"--force_gpa")){
			*force_gpa = atoi(argc[++i]);
		}else if(!strcmp(argc[i],"--force_pcg")){
			if((*force_pcg) == NULL){
				(*force_pcg) = (float*)malloc(3*sizeof(float));

				if((*force_pcg) != NULL){
					(*force_pcg)[0] = atof(argc[++i]);
					(*force_pcg)[1] = atof(argc[++i]);
					(*force_pcg)[2] = atof(argc[++i]);
				}else{
					printf("ERROR: could not allocate memory for force_pcg\n");
					exit(2);
				}
			}
		}else if(!strcmp(argc[i],"--atom_index")){
			*atom_index = atoi(argc[++i]);

			/*else if(!strcmp(argc[i],"--old_types")){
			 *old_types = 1;
			 }else if(!strcmp(argc[i],"--new_types")){
			 *new_types = 1;
			 }else if(!strcmp(argc[i],"--sybyl_types")){
			 *sybyl_types = 1;
			 }*/

		}else if(!strcmp(argc[i],"--res_name")){
			strncpy(force_outres->name,argc[++i],3);
			force_outres->name[3] = '\0';
		}else if(!strcmp(argc[i],"--res_number")){
			force_outres->number = atoi(argc[++i]);
		}else if(!strcmp(argc[i],"--res_chain")){
			force_outres->chain = argc[++i][0];
		}
	}

}

void Print_subGraph(subgraph* graph, atom* atoms, int n_atoms){

	while(graph != NULL){

		printf("graph[%d]: ", graph->id);
		for(int i=0; i<n_atoms; i++){
			if(atoms[i].graph == graph){
				printf("%5d", atoms[i].number);
			}
		}
		printf("\n");

		graph = graph->prev;
	}

}


void Print_Connections(atom* atoms, int n_atoms){

	int i,j;

	for(i=0; i<n_atoms; ++i){
		printf("Connections for atom %d:",atoms[i].number);
		for(j=1; j<=atoms[i].n_bonds; ++j){
			printf(" %d(%d)",atoms[i].conect[j].to->number, atoms[i].conect[j].type);
		}
		printf("\n\n");
	}

}

int is_Terminal(bond* conect){

	int z1,z2;

	z1 = conect->to->n_bonds-count_Hydrogens(conect->to);
	z2 = conect->from->n_bonds-count_Hydrogens(conect->from);

	return (z1 == 1 || z2 == 1) ? 1 : 0;
}

int is_Hydrogen(atom* atomz){

	int i=0;

	if(atomz->type[i] == ' ')
		while(atomz->type[i] == ' ' && atomz->type[i] != '\0')
			i++;

	if(atomz->type[i] == 'H')
		return 1;
	else
		return 0;

}

int memAlloc_Paths(atom** atoms, int n_atoms){

	int i,j;

	for(i=0; i<n_atoms; ++i){

		(*atoms)[i].sp_n_paths = 0;
		(*atoms)[i].sp_paths_n = 0;
		(*atoms)[i].sp_state = 0;

		(*atoms)[i].sp_paths = (atom***)malloc(NPATHS*sizeof(atom**));
		if((*atoms)[i].sp_paths == NULL){
			fprintf(stderr,"could not allocate memory for sp_paths.\n");
			return 0;
		}
		//printf("[-1]: allocated %p\n", (*atoms)[i].sp_paths);

		for(j=0; j<NPATHS; ++j){
			(*atoms)[i].sp_paths[j] = (atom**)malloc(n_atoms*sizeof(atom*));
			if((*atoms)[i].sp_paths[j] == NULL){
				fprintf(stderr,"could not allocate memory for sp_paths.\n");
				return 0;
			}
			//printf("[%2d]: allocated %p\n", j, (*atoms)[i].sp_paths[j]);
			memset((*atoms)[i].sp_paths[j],0,n_atoms*sizeof(atom*));
		}
	}

	return 1;
}

void free_Paths(atom** atoms, int n_atoms){

	int i,j;

	for(i=0; i<n_atoms; ++i){

		for(j=0; j<NPATHS; ++j){

			free((*atoms)[i].sp_paths[j]);
		}
		free((*atoms)[i].sp_paths);
	}

}

void Copy_Paths(atom* atomdest, atom* atomsrc){

	int i,j;

	for(i=0; i<atomsrc->sp_n_paths; ++i){

		for(j=0; j<atomsrc->sp_paths_n; ++j){

			atomdest->sp_paths[atomdest->sp_n_paths][j] =
				atomsrc->sp_paths[i][j];

		}

		atomdest->sp_paths[atomdest->sp_n_paths][atomsrc->sp_paths_n] = atomsrc;
		atomdest->sp_paths_n = atomsrc->sp_paths_n + 1;

		atomdest->sp_n_paths++;
	}

}

int is_Cyclic(atom* atomz, int* scc, int n_scc){

	int i;

	for(i=0; i<n_scc; i++){
		if(scc[i] == atomz->number)
			return 1;
	}

	return 0;
}

int get_Shortest_Path(atom* atomi,atom* atoms, int n_atoms){

	int i,j;
	int counter;
	int ndone = 0;

	if(!memAlloc_Paths(&atoms,n_atoms)){ return 1; }

	// Find Shortest Path to atomi

	atomi->sp_state = 2;
	atomi->sp_n_paths = 1;
	//printf("shortest path from %d\n", atomi->number);

	while(ndone < n_atoms){

		counter = 0;
		for(i=0; i<n_atoms; ++i){

			// Atom is READY-CHECK
			if(atoms[i].sp_state == 2){

				if(atoms[i].n_bonds == 0 ||
				   atoms[i].n_bonds == count_Hydrogens(&atoms[i])){
					return 0;
				}

				for(j=1; j<=atoms[i].n_bonds; ++j){

					// Atom is not set to READY
					if(atoms[i].conect[j].to->sp_state == 0 ||
					   atoms[i].conect[j].to->sp_state == 1){

						Copy_Paths(atoms[i].conect[j].to,&atoms[i]);

						atoms[i].conect[j].to->sp_state = 1;

					}

				}

				counter++;
			}

		}

		for(i=0; i<n_atoms; ++i){

			if(atoms[i].sp_state == 1){
				atoms[i].sp_state = 2;
			}else if(atoms[i].sp_state == 2){
				atoms[i].sp_state = 3;
				ndone++;
			}
		}

		if(counter == 0) { return 0; }
	}

	return 1;

}

void Print_Paths(atom* atoms, int n_atoms){

	int i,j,k;

	for(i=0; i<n_atoms; ++i){

		printf("Paths for atom %d (%d):\n", atoms[i].number, atoms[i].sp_paths_n);
		for(j=0; j<atoms[i].sp_n_paths; ++j){
			printf("[%d]:",j);
			for(k=0; k<atoms[i].sp_paths_n; ++k){
				printf(" %d",atoms[i].sp_paths[j][k]->number);
			}
			printf("\n");
		}
		printf("\n\n");
	}
}

int is_NonMetal(char* element) {

	if(!strcmp(element,"H") || !strcmp(element,"D") || !strcmp(element,"B") || !strcmp(element,"C") ||
	   !strcmp(element,"N") || !strcmp(element,"O") || !strcmp(element,"F") || !strcmp(element,"Si") ||
	   !strcmp(element,"P") || !strcmp(element,"S") || !strcmp(element,"Cl") || !strcmp(element,"As") ||
	   !strcmp(element,"Se") || !strcmp(element,"Br") || !strcmp(element,"Te") || !strcmp(element,"I") ||
	   !strcmp(element,"At") || !strcmp(element,"He") || !strcmp(element,"Ne") || !strcmp(element,"Ar") ||
	   !strcmp(element,"Kr") || !strcmp(element,"Xe") || !strcmp(element,"Rn")){

		return 1;
	}

	return 0;
}


int is_Extractable(residue* extract, int n_extract, char (*args)[25]){

	int i;

	if(extract != NULL){

		for(i=0; i<n_extract; i++){
			if(extract[i].number == atoi(args[6]))
				return 1;
		}

	}else{
		return 1;
	}

	return 0;
}


atom* read_MOL2(char* filename, int* n_atoms, int* map_atom, residue *extract, int n_extract, float* ori_pcg, int atom_index){

	//@<TRIPOS>ATOM
	//@<TRIPOS>BOND
	//1  N         19.0030    9.2780   15.3460 N.3    63  ALA63      -0.3202

	int i=0;

	atom* MOL2 = NULL;
	FILE* infile_ptr = NULL;
	char buffer[100];
	char* pch;
	char fields[9][25];

	int atom_flag = 0;
	int bond_flag = 0;
	int n_ori = 0;

	int   n_malloc = 25;

	infile_ptr = fopen(filename,"r");
	if(infile_ptr == NULL){
		fprintf(stderr, "could not open MOL2 file: %s\n", filename);
		return MOL2;
	}

	MOL2 = (atom*)malloc(n_malloc*sizeof(atom));

	if(MOL2 == NULL){
		fprintf(stderr, "could not allocate memory for atoms\n");
		return MOL2;
	}

	while(fgets(buffer,sizeof(buffer),infile_ptr) != NULL){

		if(!strncmp("@<TRIPOS>ATOM",buffer,13)){
			atom_flag = 1;
			bond_flag = 0;
			continue;
		}else if(!strncmp("@<TRIPOS>BOND",buffer,13)){
			atom_flag = 0;
			bond_flag = 1;
			continue;
		}else if(!strncmp("@<TRIPOS>",buffer,9)){
			atom_flag = 0;
			bond_flag = 0;
			continue;			
		}

		if(atom_flag){

			//0         1         2         3         4         5         6         7
			//01234567890123456789012345678901234567890123456789012345678901234567890123456789
			//     11  NAX       -8.1990  -30.9250   23.7700 N.pl32070  ZZP2070    -0.2965
			//      3  O2G       41.2980   20.9880   25.9490 O.co2 501  ANP501     -0.5861
			//      7  O2B       45.7340   19.3320   23.7330 O.3   501  ANP501     -0.3119

			strncpy(fields[0],&buffer[1],6);            // atom number
			fields[0][6]='\0';
			strncpy(fields[1],&buffer[8],4);            // atom name
			fields[1][4]='\0';
			strncpy(fields[2],&buffer[16],10);          // x-coord
			fields[2][10]='\0';
			strncpy(fields[3],&buffer[26],10);          // y-coord
			fields[3][10]='\0';
			strncpy(fields[4],&buffer[36],10);          // z-coord
			fields[4][10]='\0';
			strncpy(fields[5],&buffer[47],5);           // hybridation
			fields[5][5]='\0';
			strncpy(fields[6],&buffer[52],4);           // residue number
			fields[6][4]='\0';
			strncpy(fields[7],&buffer[58],7);           // residue name/number
			fields[7][7]='\0';
			strncpy(fields[8],&buffer[66],10);          // charge
			fields[8][10]='\0';

			//     ---- FIELDS ----
			/*
			  36
			  H
			  15.4394
			  2.2921
			  20.8702
			  H
			  65
			  ALA65
			  0.0252
			*/

			/*
			  for(i=0; i<9; i++)
			  printf("fields[%d]=%s\n",i,fields[i]);
			  printf("\n");
			*/

			if(ori_pcg != NULL){
				for(i=0; i<3; i++)
					ori_pcg[i] += atof(fields[2+i]);
				n_ori++;
			}


			if(!is_Extractable(extract,n_extract,fields)){
				continue;
			}


			if(*n_atoms == n_malloc){
				n_malloc += 5;

				MOL2 = (atom*)realloc(MOL2,n_malloc*sizeof(atom));

				if(MOL2 == NULL){
					fprintf(stderr, "could not reallocate memory for atoms\n");
					return MOL2;
				}
			}


			MOL2[*n_atoms].vertex.index = -1;
			MOL2[*n_atoms].vertex.lowlink = -1;

			MOL2[*n_atoms].number = atoi(fields[0]);
			MOL2[*n_atoms].index = atom_index > 0 ? atom_index + *n_atoms : MOL2[*n_atoms].number;
			//printf("old index: %d - new index: %d\n", MOL2[*n_atoms].number, MOL2[*n_atoms].index);

			strcpy(MOL2[*n_atoms].name,fields[1]);

			MOL2[*n_atoms].coor[0] = atof(fields[2]);
			MOL2[*n_atoms].coor[1] = atof(fields[3]);
			MOL2[*n_atoms].coor[2] = atof(fields[4]);

			MOL2[*n_atoms].conect = (bond*)malloc(MAX_BONDS*sizeof(bond));
			if(MOL2[*n_atoms].conect == NULL){
				fprintf(stderr, "could not allocate memory for bond\n");
				exit(2);
			}

			MOL2[*n_atoms].n_bonds=0;

			MOL2[*n_atoms].gpa=0;
			MOL2[*n_atoms].build_state=0;
			MOL2[*n_atoms].shortest=1000;

			MOL2[*n_atoms].graph = NULL;

			MOL2[*n_atoms].next_build = NULL;
			MOL2[*n_atoms].sp_paths = NULL;
			MOL2[*n_atoms].aromatic = 0;

			sscanf(fields[5],"%s",MOL2[*n_atoms].type);
			for(int j=0; j<6; j++){
				MOL2[*n_atoms].type[j] = toupper(MOL2[*n_atoms].type[j]);
			}

			get_Element_From_Hybridation(MOL2[*n_atoms].type,MOL2[*n_atoms].element);
			MOL2[*n_atoms].nonmetal = is_NonMetal(MOL2[*n_atoms].element);

			//set_AtomTypes_SYBYL(&MOL2[*n_atoms],0);
			
			MOL2[*n_atoms].charge = atof(fields[8]);

			map_atom[MOL2[*n_atoms].number] = *n_atoms;

			(*n_atoms)++;

		}else if(bond_flag){

			i=0;

			pch = strtok (buffer," ");
			while (pch != NULL)
			{
				//printf ("%s\n",pch);
				strcpy(fields[i++],pch);
				pch = strtok (NULL, " ");
			}

			//     ---- FIELDS ----
			//     1    11    12    1
			//     2    11    10    1

			if(is_Mapped(atoi(fields[1]),map_atom) && is_Mapped(atoi(fields[2]),map_atom)){

				MOL2[map_atom[atoi(fields[1])]].n_bonds++;
				MOL2[map_atom[atoi(fields[2])]].n_bonds++;

				MOL2[map_atom[atoi(fields[1])]].conect[MOL2[map_atom[atoi(fields[1])]].n_bonds].to = &MOL2[map_atom[atoi(fields[2])]];
				MOL2[map_atom[atoi(fields[1])]].conect[MOL2[map_atom[atoi(fields[1])]].n_bonds].from = &MOL2[map_atom[atoi(fields[1])]];

				MOL2[map_atom[atoi(fields[2])]].conect[MOL2[map_atom[atoi(fields[2])]].n_bonds].to = &MOL2[map_atom[atoi(fields[1])]];
				MOL2[map_atom[atoi(fields[2])]].conect[MOL2[map_atom[atoi(fields[2])]].n_bonds].from = &MOL2[map_atom[atoi(fields[2])]];

				MOL2[map_atom[atoi(fields[1])]].conect[MOL2[map_atom[atoi(fields[1])]].n_bonds].type = read_Type(fields[3]);
				MOL2[map_atom[atoi(fields[2])]].conect[MOL2[map_atom[atoi(fields[2])]].n_bonds].type = read_Type(fields[3]);

				MOL2[map_atom[atoi(fields[1])]].conect[MOL2[map_atom[atoi(fields[1])]].n_bonds].cyclic = 0;
				MOL2[map_atom[atoi(fields[2])]].conect[MOL2[map_atom[atoi(fields[2])]].n_bonds].cyclic = 0;

				MOL2[map_atom[atoi(fields[1])]].conect[MOL2[map_atom[atoi(fields[1])]].n_bonds].flexible = 0;
				MOL2[map_atom[atoi(fields[2])]].conect[MOL2[map_atom[atoi(fields[2])]].n_bonds].flexible = 0;

				MOL2[map_atom[atoi(fields[1])]].conect[MOL2[map_atom[atoi(fields[1])]].n_bonds].dist =
					dist(MOL2[map_atom[atoi(fields[1])]].coor,
					     MOL2[map_atom[atoi(fields[1])]].conect[MOL2[map_atom[atoi(fields[1])]].n_bonds].to->coor);
				MOL2[map_atom[atoi(fields[2])]].conect[MOL2[map_atom[atoi(fields[2])]].n_bonds].dist =
					dist(MOL2[map_atom[atoi(fields[2])]].coor,
					     MOL2[map_atom[atoi(fields[2])]].conect[MOL2[map_atom[atoi(fields[2])]].n_bonds].to->coor);
			}

		}else{
			// do nothing
		}

	}

	if(n_ori > 0)
		for(i=0; i<3; i++)
			ori_pcg[i] /= (float)n_ori;


	return MOL2;

}

void get_Element_From_Hybridation(char* string, char* dest){

	if(!strncmp(string,"O.CO2",5)  || !strncmp(string,"O.AR",4) || !strncmp(string,"O.2",3) || !strncmp(string,"O.3",3)){
		strcpy(dest,"O");
	}else if(!strncmp(string,"P.3",3)){
		strcpy(dest,"P");
	}else if(!strncmp(string,"S.O2",4) || !strncmp(string,"S.2",3) || !strncmp(string,"S.3",3) ||
		 !strncmp(string,"S.O",3)  || !strncmp(string,"S.AR",4)){
		strcpy(dest,"S");
	}else if(!strncmp(string,"N.AM",4) || !strncmp(string,"N.AR",4) || !strncmp(string,"N.PL3",5) ||
		 !strncmp(string,"N.1",3) || !strncmp(string,"N.2",3) || !strncmp(string,"N.3",3) || !strncmp(string,"N.4",3)){
		strcpy(dest,"N");
	}else if(!strncmp(string,"C.CAT",5) || !strncmp(string,"C.AR",4) || !strncmp(string,"C.1",3) ||
		 !strncmp(string,"C.2",3) || !strncmp(string,"C.3",3)){
		strcpy(dest,"C");
	}else if(!strncmp(string,"CL",2)){
		strcpy(dest,"Cl");
	}else if(!strncmp(string,"BR",2)){
		strcpy(dest,"Br");
	}else if(!strncmp(string,"SE",2)){
		strcpy(dest,"Se");
	}else if(!strncmp(string,"MG",2)){
		strcpy(dest,"Mg");
	}else if(!strncmp(string,"SR",2)){
		strcpy(dest,"Sr");
	}else if(!strncmp(string,"CU",2)){
		strcpy(dest,"Cu");
	}else if(!strncmp(string,"MN",2)){
		strcpy(dest,"Mn");
	}else if(!strncmp(string,"MG",2)){
		strcpy(dest,"Mg");
	}else if(!strncmp(string,"HG",2)){
		strcpy(dest,"Hg");
	}else if(!strncmp(string,"CD",2)){
		strcpy(dest,"Cd");
	}else if(!strncmp(string,"NI",2)){
		strcpy(dest,"Ni");
	}else if(!strncmp(string,"ZN",2)){
		strcpy(dest,"Zn");
	}else if(!strncmp(string,"CA",2)){
		strcpy(dest,"Ca");
	}else if(!strncmp(string,"FE",2)){
		strcpy(dest,"Fe");
	}else if(!strncmp(string,"F",1)){
		strcpy(dest,"F");
	}else if(!strncmp(string,"I",1)){
		strcpy(dest,"I");
	}else if(!strncmp(string,"H",1)){
		strcpy(dest,"H");
	}else if(!strncmp(string,"CO.OH",5)){
		strcpy(dest,"Co");
	}else{
		printf("WARNING: unknown element for type %s\n", string);
		strcpy(dest,"  ");
	}

}

int read_Type(char* type){

	if(!strncmp(type,"1",1)) { return 1; }
	else if(!strncmp(type,"2",1)) { return 2; }
	else if(!strncmp(type,"3",1)) { return 3; }
	else if(!strncmp(type,"4",1)) { return 4; }
	else if(!strncmp(type,"ar",2)) { return 5; }
	else if(!strncmp(type,"am",2)) { return 9; }
	else {
		return 10;
	}

}


int is_Mapped(int number, int* map_atom){

	return ( map_atom[number] != -1 ? 1 : 0 );

}

int Get_NextConnection(char* buffer){

	int   i;
	int   num;
	char  bufnul[6];

	for(i=0; i<5; ++i)
	{
		if(buffer[i] != ' '){
			strncpy(bufnul,buffer,5);
			bufnul[5]='\0';

			num = atoi(bufnul);
			//printf("connected to %d\n",num);
			return(num);
		}
	}

	return 0;

}

int Bond_Exists(bond* b, bond* blist[], int nb){

	int i;

	for(i=0; i<nb; i++){
		if((blist[i]->from->number == b->from->number &&
		    blist[i]->to->number == b->to->number) ||
		   (blist[i]->from->number == b->to->number &&
		    blist[i]->to->number == b->from->number))
		{
			return 1;
		}
	}

	return 0;
}

void get_Flexible_Atoms(atom* atoms,int n_atoms,bond* flex,atom* atomlist[], int* sense, int *nlist, int remove_hydro){

	int i;
	//int mobile;

	for(i=0; i<n_atoms; i++){

		if(remove_hydro && is_Hydrogen(&atoms[i]))
			continue;

		if((atoms[i].buildlist[0]->number == flex->to->number &&
		    atoms[i].buildlist[1]->number == flex->from->number) ||
		   (atoms[i].buildlist[0]->number == flex->from->number &&
		    atoms[i].buildlist[1]->number == flex->to->number))
		{

			if(atoms[i].buildlist[0]->number == 0 ||
			   atoms[i].buildlist[1]->number == 0 ||
			   atoms[i].buildlist[2]->number == 0)
				continue;

			sense[(*nlist)] = atoms[i].buildlist[0]->number == flex->to->number ? 1 : -1;
			atomlist[(*nlist)] = &atoms[i];
			printf("added atom %d to atomlist\n", atoms[i].index);

			(*nlist)++;
		}
	}

}

atom* get_Atom_from_coor(const float* coor, atom* atoms, int n_atoms){

	for(int i=0; i<n_atoms; i++){
		int match = 1;
		for(int j=0; j<3; j++){
			if(fabs(atoms[i].coor[j]-coor[j]) > 0.0001){
				match = 0;
				break;
			}
		}
		if(match){return &atoms[i];}
	}

	return NULL;
}

void Write_Target(char* original_filename, char* target_filename,
		  atom* atoms, int n_atoms){

	FILE* infile_ptr = fopen(original_filename,"r");
	FILE* outfile_ptr = fopen(target_filename,"w");

	char buffer[100];

	if(infile_ptr == NULL){
		fprintf(stderr,"could not open file %s for reading\n", original_filename);
		return;
	}

	if(outfile_ptr == NULL){
		fprintf(stderr,"could not open file %s for writing\n", target_filename);
		return;
	}

	while(fgets(buffer,sizeof(buffer),infile_ptr) != NULL){
		if(!strncmp(buffer,"ATOM  ",6) || !strncmp(buffer,"HETATM",6)){
			float coor[3];
			for(int i=0; i<3; i++){
				coor[i] = atof(&buffer[30+i*8]);
			}

			atom* atomz = get_Atom_from_coor(coor,atoms,n_atoms);
			if(atomz == NULL){
				fprintf(stderr, "ERROR: No matching atoms for coordinates.\n");
			}else{
				//0         1         2         3         4         5         6         7
				//01234567890123456789012345678901234567890123456789012345678901234567890123456789
				//HETATM  709 MG    MG     2       8.551  -2.525  11.085  1.00  8.13          MG
				buffer[76] = '\0';
				fprintf(outfile_ptr, "%s%2d\n", buffer, atomz->atomtype);
			}
		}else{
			fprintf(outfile_ptr, "%s", buffer);
		}
	}

	fclose(outfile_ptr);
}

void Write_IC(char* filename, atom* atoms, int n_atoms, atom* gpa, int remove_hydro){

	int i;
	FILE* infile_ptr;
	atom* atomw;

	infile_ptr = fopen(filename,"w");

	if(infile_ptr == NULL)
	{
		fprintf(stderr,"could not open file %s for writing\n", filename);
		return;
	}


	//    1:    1.462  112.662   82.309
	for(i=0; i<n_atoms; i++)
	{
		atomw = &atoms[i];

		if(remove_hydro && is_Hydrogen(atomw)) { continue; }

		fprintf(infile_ptr, "%5d: %8.3f %8.3f %8.3f\n",
			atomw->index,
			atomw->dis=dist(atomw->coor,atomw->buildlist[0]->coor),
			atomw->ang=bndang(atomw->coor,atomw->buildlist[0]->coor,atomw->buildlist[1]->coor),
			atomw->dih=dihedral(atomw->coor,atomw->buildlist[0]->coor,atomw->buildlist[1]->coor,atomw->buildlist[2]->coor));
	}

	fprintf(infile_ptr, "REFPCG %8.3f %8.3f %8.3f\n",
		gpa->buildlist[1]->coor[0],
		gpa->buildlist[1]->coor[1],
		gpa->buildlist[1]->coor[2]);

	fclose(infile_ptr);

}


void Write_INP(char* filename, char* icfile, atom* atoms, int n_atoms, int remove_hydro,
	       residue* outres, residue* force_outres, atom* gpa, subgraph* graph){

	int i,j;
	FILE* infile_ptr;
	atom* atomw;

	bond* flex[MAX_CYCLE];
	int n_flex=0;
	int validated=1;

	atom* atomlist[MAXCONECT];
	int sense[MAXCONECT];
	int nlist=0;

	infile_ptr = fopen(filename,"w");

	if(infile_ptr == NULL){
		fprintf(stderr,"could not open file %s for writing\n", filename);
		return;
	}

	fprintf(infile_ptr, "RESIDU %3s %c %d\n",
		strcmp(force_outres->name,"   ") ? force_outres->name : outres->name,
		force_outres->chain != '-' ? force_outres->chain : outres->chain,
		force_outres->number != -1 ? force_outres->number : outres->number);

	for(i=0; i<n_atoms; i++){
		atomw = &atoms[i];

		if(remove_hydro && is_Hydrogen(atomw)) { continue; }

		fprintf(infile_ptr, "HETTYP%5d%2d %4s %c %5d%5d%5d%5d\n",
			atomw->index,
			atomw->atomtype,
			atomw->name,
			'm',
			atomw->buildlist[0]->number == 0 ? 0 : atomw->buildlist[0]->index,
			atomw->buildlist[1]->number == 0 ? 0 : atomw->buildlist[1]->index,
			atomw->buildlist[2]->number == 0 ? 0 : atomw->buildlist[2]->index,
			atomw->graph->id
			);
	}

	fprintf(infile_ptr, "GPATOM %d %d %d\n",
		gpa->index,
		gpa->next_build->index,
		gpa->next_build->next_build->index);


	for(i=0; i<n_atoms; i++){
		atomw = &atoms[i];

		if(remove_hydro && is_Hydrogen(atomw)) { continue; }

		for(j=1; j<=atomw->n_bonds; j++){
			if(atomw->conect[j].flexible){
				if(!Bond_Exists(&atomw->conect[j],flex,n_flex)){
					{
						flex[n_flex++] = &atomw->conect[j];
					}
				}
			}
		}
	}


	//printf("number of unique flexible bonds: %d\n", n_flex);
	int n_flex_counter = 1;
	for(i=0; i<n_flex; i++){
		nlist=0;
		get_Flexible_Atoms(atoms,n_atoms,flex[i],atomlist,sense,&nlist,remove_hydro);
		
		if(nlist == 0){
			printf("the flexible bond [%d][%d] could not be set\n",
			       flex[i]->to->number, flex[i]->from->number);
		}else{
			validated = sense[0];
			for(j=1; j<nlist; j++){
				if(validated != sense[j]) {
					validated = 0;
				}
			}

			if(validated){
				fprintf(infile_ptr,"FLEDIH%3d ", n_flex_counter++);
				for(j=0; j<nlist; j++){ fprintf(infile_ptr,"%5d", atomlist[j]->index); }
				fprintf(infile_ptr,"\n");
			}else{
				printf("could not define flexible bond [%d-%d]\n",
				       flex[i]->to->number, flex[i]->from->number);
			}
		}
	}


	for(i=0; i<n_atoms; i++){
		atomw = &atoms[i];

		if(remove_hydro && is_Hydrogen(atomw)) { continue; }

		fprintf(infile_ptr, "CONECT%5d", atomw->index);
		for(j=1; j<=atomw->n_bonds; j++){
			if(remove_hydro && is_Hydrogen(atomw->conect[j].to)) { continue; }

			fprintf(infile_ptr, "%5d", atomw->conect[j].to->index);
		}
		fprintf(infile_ptr, "\n");
	}
	
	fprintf(infile_ptr,"ENDINP\n");

	fclose(infile_ptr);

}

void Write_REF(char* filename, atom* atoms, int n_atoms, int remove_hydro, residue* outres, residue* force_outres){

	FILE* infile_ptr = NULL;

	infile_ptr = fopen(filename,"w");

	if(infile_ptr == NULL){
		fprintf(stderr,"could not open file %s for writing\n", filename);
		return;
	}

	for(int i=0; i<n_atoms; i++){

		if(remove_hydro && is_Hydrogen(&atoms[i])) { continue; }

		fprintf(infile_ptr,"HETATM%5d %4s %3s %c%4d    %8.3f%8.3f%8.3f  1.00  1.00          %2s\n",
			atoms[i].index,
			atoms[i].name,
			strcmp(force_outres->name,"   ") ? force_outres->name : outres->name,
			force_outres->chain != '-' ? force_outres->chain : outres->chain,
			force_outres->number != -1 ? force_outres->number : outres->number,
			atoms[i].coor[0],atoms[i].coor[1],atoms[i].coor[2],
			atoms[i].element
			);
	}

	for(int i=0; i<n_atoms; i++){

		if(remove_hydro && is_Hydrogen(&atoms[i])) { continue; }
		fprintf(infile_ptr, "CONECT%5d", atoms[i].index);

		for(int j=1; j<=atoms[i].n_bonds; j++){

			if(remove_hydro && is_Hydrogen(atoms[i].conect[j].to)) { continue; }
			fprintf(infile_ptr, "%5d", atoms[i].conect[j].to->index);
		}

		fprintf(infile_ptr, "\n");
	}


	fprintf(infile_ptr, "END\n");

	fclose(infile_ptr);

}

void buildcc(atom* sequence){

	float x[4],y[4],z[4];
	float a,b,c,op,cx,cy,cz,d,xn,yn,zn,ct,st,xk,yk,zk,angPI,dihPI;

	do{
		//printf("buildcc of %d\n", sequence->number);

		for(int i=1; i<=3; i++){
			x[i] = sequence->buildlist[i-1]->coor[0];
			y[i] = sequence->buildlist[i-1]->coor[1];
			z[i] = sequence->buildlist[i-1]->coor[2];

			//perturb atom coordinates
			x[i] += 1e-10;
			y[i] += 1e-10;
			z[i] += 1e-10;
		}

		a=y[1]*(z[2]-z[3])+y[2]*(z[3]-z[1])+y[3]*(z[1]-z[2]);
		b=z[1]*(x[2]-x[3])+z[2]*(x[3]-x[1])+z[3]*(x[1]-x[2]);
		c=x[1]*(y[2]-y[3])+x[2]*(y[3]-y[1])+x[3]*(y[1]-y[2]);

		op=sqrtf(a*a+b*b+c*c);

		cx=a/op;
		cy=b/op;
		cz=c/op;

		a=x[2]-x[1];
		b=y[2]-y[1];
		c=z[2]-z[1];

		d=1.0f/sqrtf(a*a+b*b+c*c);
		op=sequence->dis*d;
		xn=a*op;
		yn=b*op;
		zn=c*op;

		a=cx*cx;
		b=cy*cy;
		c=cz*cz;

		angPI = sequence->ang*PI/180.0f;
		ct=cos(angPI);
		st=-sin(angPI);

		op=1.0f-ct;

		xk=(cx*cz*op-cy*st)*zn+((1.0f-a)*ct+a)*xn+(cx*cy*op+cz*st)*yn;
		yk=(cy*cx*op-cz*st)*xn+((1.0f-b)*ct+b)*yn+(cy*cz*op+cx*st)*zn;
		zk=(cz*cy*op-cx*st)*yn+((1.0f-c)*ct+c)*zn+(cz*cx*op+cy*st)*xn;

		dihPI = sequence->dih*PI/180.0f;
		ct=cos(dihPI);
		st=sin(dihPI);

		op=1.0f-ct;

		cx=(x[2]-x[1])*d;
		cy=(y[2]-y[1])*d;
		cz=(z[2]-z[1])*d;

		a=cx*cx;
		b=cy*cy;
		c=cz*cz;

		x[0]=(cx*cz*op-cy*st)*zk+((1.0f-a)*ct+a)*xk+(cx*cy*op+cz*st)*yk+x[1];
		y[0]=(cy*cx*op-cz*st)*xk+((1.0f-b)*ct+b)*yk+(cy*cz*op+cx*st)*zk+y[1];
		z[0]=(cz*cy*op-cx*st)*yk+((1.0f-c)*ct+c)*zk+(cz*cx*op+cy*st)*xk+z[1];

		sequence->coor[0]=x[0];
		sequence->coor[1]=y[0];
		sequence->coor[2]=z[0];

		//printf("build sequence atom number %d\n", sequence->number);
		sequence = sequence->next_build;

	}while(sequence != NULL);

}
