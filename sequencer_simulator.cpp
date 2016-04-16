/* This code is shared with no explicit or implicit guarantee about 
 * functionalities. It's part of a thesis. For more information:
 * http://pikkolamakkia.com/thesis
 * 
 * Author: Carlo Masaia
 * Version: 1.0-thesis
 * Last update: 2016.Feb.12
 * Email: karurochari@mail.com
 * 
 * The rights of the authors on this code are enforced by Italian and 
 * European laws on copyright. Every kind of usage, alteration, copy, 
 * printing of what it follows is strictly forbidden without the explicit
 * permission of the authors.
*/

#include <iostream>
#include <fstream>
#include <cstdlib> 
#include <ctime>
#include <cstring>

using namespace std;

//DNA sequencer simulator.

inline char to_base(int i){if(i==0)return 'A';else if(i==1)return 'T';else if(i==2)return 'C';else return 'G';}

int main(int argc, char *argv[]){
	if(argc!=6){
		cout<<"Syntax error! Good syntaxes are the following.\n\n"
			<<"Read DNA from file [WIP]:\n"
			<<argv[0]<<" -f {filename} {reads} {reads-len} {output}\n"
			<<"Generate DNA by internal source:\n"
			<<argv[0]<<" -g {dna-length} {reads} {reads-len} {output}\n";
			
		return -1;
	}
	
	srand (time(NULL));
	
	int reads_n=atoi(argv[3]);
	int reads_len=atoi(argv[4]);
	int DNA_len=0;

	
	char *DNA=nullptr;
	
	if(strcmp(argv[1],"-g")==0){
		DNA_len=atoi(argv[2]);
		if(reads_len>DNA_len)throw "Reads bigger than DNA";
		DNA = new char[DNA_len+1];
		DNA[DNA_len]=0;
	
		for(int i=0;i<DNA_len;i++)DNA[i]=to_base(rand()%4);
	}
	
	else if(strcmp(argv[1],"-f")==0){
		fstream infile(argv[2],ios::in);
		if(!infile.is_open()){cout<<"File not found: "<<argv[2]<<"\n";return 0;}
		infile>>DNA_len;

		if(reads_len>DNA_len)throw "Reads bigger than DNA";

		DNA = new char[DNA_len+1];
		DNA[DNA_len]=0;	
		
		int j=0;
		for(int i=0;infile.good();i++){char c;infile.get(c);if(c!='\n'){DNA[j]=c;j++;}}	
	}
	
	char *reads = new char[(reads_len+1)*reads_n];
	
	for(int i=0;i<reads_n;i++){
		memcpy(&reads[i*(reads_len+1)],&DNA[rand()%(DNA_len-reads_len)],reads_len);
		reads[(i+1)*(reads_len+1)-1]='\n';
	}
	reads[(reads_len+1)*(reads_n+1)-1]=0;
	
	string fname=argv[5];
	
	{
	fstream f(fname+".reads.txt",ios::out);
	f<<reads_n<<" "<<reads_len<<"\n";
	f<<reads;
	f.close();
	}
	
	{
	fstream f(fname+".dna.txt",ios::out);
	f<<DNA_len<<"\n";
	f<<DNA;
	f.close();
	}
	
	delete[] DNA;
	delete[] reads;
	
	return 0;
}
