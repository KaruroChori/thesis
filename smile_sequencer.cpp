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
#include <cstring>
#include <fstream>
#include <vector>
#include <algorithm>
#include <thread>
#include <atomic>

#include <unistd.h>

#include "headers/dna_read.h"
#include "headers/metric_matrix.h"
#include "headers/linear_optimizer.h"
#include "headers/precompressor.h"


using namespace std;

//Import from the other file. It may be a built in one or user defined
struct query{
    static string GREEDY, GREEDY_STAR, GREEDY_N, GEN_AND_TEST, PARTICLE;
};

string query::GREEDY("select MIN(0,1) take 1 with RAND");
string query::GREEDY_STAR("select MIN(0,1) take ALL with ANY");
string query::GREEDY_N("select MIN(0,n) take ALL with ASC_ORDER");		//DELETE IT!
string query::GEN_AND_TEST("select IDENTITY take ALL with ASC_ORDER");
string query::PARTICLE("select IDENTITY take 1 with RAND_P_INV");

/*
Greedy      select MIN(0,1) take 1    order by RAND
Greedy*     select MIN(0,1) take ALL order by ANY
Greedy(n)   select MIN(0,n) take ALL order by ASC_ORDER
Gen&Test    select IDENTITY take ALL order by ASC_ORDER
Particle    select IDENTITY take 1 order by RAND_P_INV
*/

struct mytime{
    float time;
    enum{us,ms,s,m,h,d,ANY_TIME} suffix;
   
    mytime(){}
    mytime(const char* str){
        if(strcmp("any", str)==0){suffix=ANY_TIME;time=0;}
        else{
			char number[16];
			char unit[3];
			int i=0;
			for(;((str[i]>='0' and str[i]<='9') or str[i]=='.') and i<15;i++)number[i]=str[i];
			number[i]=0;
			int j=i;
			for(;str[j]!='\0' and j<2+i;j++){unit[j-i]=str[j];}
			unit[j-i]=0;
			time=atof(number);
			
			if(strcmp(unit,"us")==0)suffix=us;
			else if(strcmp(unit,"ms")==0)suffix=ms;
			else if(strcmp(unit,"s")==0)suffix=s;
			else if(strcmp(unit,"m")==0)suffix=m;
			else if(strcmp(unit,"h")==0)suffix=h;
			else if(strcmp(unit,"d")==0)suffix=d;
			else if(strcmp(unit,"")==0)suffix=s;
			else throw "Time suffix not valid!";
			
        }
    }
    friend ostream& operator<<(ostream& out, mytime& t){
		if(t.suffix==ANY_TIME){out<<"no-limit";return out;}
		out<<t.time;
		if(t.suffix==us)out<<"us";
		else if(t.suffix==ms)out<<"ms";
		else if(t.suffix==s)out<<"s";
		else if(t.suffix==m)out<<"m";
		else if(t.suffix==h)out<<"h";
		else if(t.suffix==d)out<<"d";
		return out;
	}   
	
	uint64_t to_us(){
		if(suffix==ANY_TIME){return ~(uint64_t)0;}
		else if(suffix==us){return time;}
		else if(suffix==ms){return time*1000;}
		else if(suffix==s){return time*1000*1000;}
		else if(suffix==m){return time*1000*1000*60;}
		else if(suffix==h){return time*1000*1000*60*60;}
		else if(suffix==d){return time*1000*1000*60*60*24;}		
		else return 0;
	}
};

struct program_config{
    //Input file;
    string input_file;
    enum{SMILE,FASTA,FASTQ} input_format;
   
    //Clustering
    enum{NONE,BARCODE} clustering_clustering;
    enum{CHI_SQUARED,CHI_SQUARED_C} clustering_test;
    float clustering_threshold;
   
    //Metric matrix
    enum{DELTA,HEAD_TAIL,K_HEAD_TAIL,SHIFT_HAMMING} matrix_metric;
    string matrix_metric_class;
    enum{TIME,SPACE} matrix_opt_profile;
   
    //Optimizer
    string optimizer_query;
    int optimizer_instances;
    mytime optimizer_timeout;
   
    //Output
    bool output_keep_order;
    bool output_generate_matrix;
    enum{ASCII,BINARY,HTML} output_output_format;
   
    //Final compression (NONE=0) already declared in contex
    enum compression_method{DEFLATE=1}        compression_section[6];
                       
    bool already_set[20];
   
    program_config(){
        for(int i=0;i!=20;i++)already_set[i]=false;
       
        input_file="";
       
        input_format=SMILE;
        clustering_clustering=NONE;
        clustering_test=CHI_SQUARED_C;
        clustering_threshold=0.005f;
        matrix_metric=HEAD_TAIL;
        matrix_metric_class="";
        matrix_opt_profile=SPACE;
        optimizer_query=query::GREEDY;
        optimizer_instances=1;
        optimizer_timeout=mytime("5s");
        output_keep_order=false;
        output_generate_matrix=false;
        output_output_format=BINARY;
        for(int i=0;i!=6;i++)compression_section[i]=DEFLATE;
    }
   
    static const char* VARS[20];
   
    static const char* IN_FORMAT[3];
    static const char* CLUSTERING_CLUSTERING[2];
    static const char* CLUSTERING_TEST[2];
    static const char* MATRIX_METRIC[4];
    static const char* MATRIX_OPT_PROFILE[2];
    static const char* BOOLEAN[2];
    static const char* OUTPUT_OUTPUT_FORMAT[3];
    static const char* COMPRESSION_SECTION[2];

   
    int test(const char* var, const char* val, const char** space, int n){
        for(int i=0;i!=n;i++){
            if(strcmp(val,space[i])==0)return i;
        }
       
        cerr<<"Error: variable {"<<var<<"} cannot accept {"<<val<<"} as input.\n";
        cerr<<"       Accepted values are";
        for(int i=0;i!=n;i++)cerr<<" {"<<space[i]<<"}";
        cerr<<"\n";
       
        return -1;
       
    }

    int set(const char* var, const char* val){
        for(int i=0;i<20;i++){
            if(strcmp(var,VARS[i])==0){
                if(already_set[i]==true)cerr<<"Warning: variable {"<<var<<"} already set. It will be overwritten!\n";
               
                if(i==0){
					if(val[0]=='{'){
						if(val[1]=='}'){input_file="";return 0;}
						input_file.resize(strlen(val));
						int k=1;
						for(int p=1;val[k]!=0;k++){
							if(val[k]=='{')p++;
							else if(val[k]=='}')p--;
							if(p==0)break;
							input_file[k-1]=val[k];
						}
						input_file[k-1]=0;
					}
                    else input_file=val;
                }
                else if(i==1){
                    int t=test(var,val,IN_FORMAT,3);
                    if(t>=0)input_format=(decltype(input_format))t;
                    else return -1;
                }
                else if(i==2){
                    int t=test(var,val,CLUSTERING_CLUSTERING,2);
                    if(t>=0)clustering_clustering=(decltype(clustering_clustering))t;
                    else return -1;
                }
                else if(i==3){
                    int t=test(var,val,CLUSTERING_TEST,2);
                    if(t>=0)clustering_test=(decltype(clustering_test))t;
                    else return -1;
                }
                else if(i==4){
                    float t=atof(val);
                    //error?
                    clustering_threshold=t;
                }
                else if(i==5){
                    int t=test(var,val,MATRIX_METRIC,4);
                    if(t>=0)matrix_metric=(decltype(matrix_metric))t;
                    else return -1;
                }
                else if(i==6){
					if(val[0]=='{'){
						if(val[1]=='}'){matrix_metric_class="";return 0;}
						matrix_metric_class.resize(strlen(val));
						int k=1;
						for(int p=1;val[k]!=0;k++){
							if(val[k]=='{')p++;
							else if(val[k]=='}')p--;
							if(p==0)break;
							matrix_metric_class[k-1]=val[k];
						}
						matrix_metric_class[k-1]=0;
					}
                    else matrix_metric_class=val;
                }
                else if(i==7){
                    int t=test(var,val,MATRIX_OPT_PROFILE,2);
                    if(t>=0)matrix_opt_profile=(decltype(matrix_opt_profile))t;
                    else return -1;
                }
               
                else if(i==8){
					if(val[0]=='{'){
						if(val[1]=='}'){optimizer_query="";return 0;}
						optimizer_query.resize(strlen(val));
						int k=1;
						for(int p=1;val[k]!=0;k++){
							if(val[k]=='{')p++;
							else if(val[k]=='}')p--;
							if(p==0)break;
							optimizer_query[k-1]=val[k];
						}
						optimizer_query[k-1]=0;
					}
                    else optimizer_query=val;
                }
                else if(i==9){
                    int t=atoi(val);
                    //error?
                    optimizer_instances=t;
                }
                else if(i==10){
                    optimizer_timeout=mytime(val);
                    //error?
                }

                else if(i==11){
                    int t=test(var,val,BOOLEAN,2);
                    if(t>=0)output_keep_order=(decltype(output_keep_order))t;
                    else return -1;
                }
                else if(i==12){
                    int t=test(var,val,BOOLEAN,2);
                    if(t>=0)output_generate_matrix=(decltype(output_generate_matrix))t;
                    else return -1;
                }
                else if(i==13){
                    int t=test(var,val,OUTPUT_OUTPUT_FORMAT,3);
                    if(t>=0)output_output_format=(decltype(output_output_format))t;
                    else return -1;
                }
                else if(i>=14 and i<14+6){
                    int t=test(var,val,COMPRESSION_SECTION,2);
                    if(t>=0)compression_section[i-13]=(compression_method)t;
                    else return -1;
                }
                               
                //Set as already set. Bad joke.
                already_set[i]=true;
               
                return 0;
            }
            else continue;
        }
       
        cerr<<"Error: variable {"<<var<<"} is not supported by this version or may be just misspelled\n";
        return -1;
    }
       
    void dump(ostream& out, bool mode){
        if(mode==true){    //Cmdline
            out<<"--in-file="<<input_file.c_str()<<" ";
            out<<"--in-format="<<IN_FORMAT[(int)input_format]<<" ";
            out<<"--clustering-mode="<<CLUSTERING_CLUSTERING[(int)clustering_clustering]<<" ";
            out<<"--clustering-test="<<CLUSTERING_TEST[(int)clustering_test]<<" ";
            out<<"--clustering-thres="<<clustering_threshold<<" ";
            out<<"--metric="<<MATRIX_METRIC[(int)matrix_metric]<<" ";
            out<<"--metric-class=\""<<matrix_metric_class.c_str()<<"\" ";
            out<<"--matrix-profile="<<MATRIX_OPT_PROFILE[(int)matrix_opt_profile]<<" ";
            out<<"--opt-query=\""<<optimizer_query.c_str()<<"\" ";
            out<<"--opt-instances="<<optimizer_instances<<" ";
            out<<"--opt-timeout="<<optimizer_timeout<<" ";
            out<<"--out-keep-order="<<BOOLEAN[output_keep_order]<<" ";
            out<<"--out-gen-matrix="<<BOOLEAN[output_generate_matrix]<<" ";
            out<<"--out-format="<<OUTPUT_OUTPUT_FORMAT[(int)output_output_format]<<" ";
            for(int i=0;i!=6;i++){
                out<<"--section["<<i<<"]="<<COMPRESSION_SECTION[(int)compression_section[i]]<<" ";
            }
        }
        else{ //File
            out<<"in-file:\t{"<<input_file.c_str()<<"}\n";
            out<<"in-format:\t"<<IN_FORMAT[(int)input_format]<<"\n";
            out<<"clustering-mode:\t"<<CLUSTERING_CLUSTERING[(int)clustering_clustering]<<"\n";
            out<<"clustering-test:\t"<<CLUSTERING_TEST[(int)clustering_test]<<"\n";
            out<<"clustering-thres:\t"<<clustering_threshold<<"\n";
            out<<"metric:\t"<<MATRIX_METRIC[(int)matrix_metric]<<"\n";
            out<<"metric-class:\t{"<<matrix_metric_class.c_str()<<"}\n";
            out<<"matrix-profile:\t"<<MATRIX_OPT_PROFILE[(int)matrix_opt_profile]<<"\n";
            out<<"opt-query:\t{"<<optimizer_query.c_str()<<"}\n";
            out<<"opt-instances:\t"<<optimizer_instances<<"\n";
            out<<"opt-timeout:\t"<<optimizer_timeout<<"\n";
            out<<"out-keep-order:\t"<<BOOLEAN[output_keep_order]<<"\n";
            out<<"out-gen-matrix:\t"<<BOOLEAN[output_generate_matrix]<<"\n";
            out<<"out-format:\t"<<OUTPUT_OUTPUT_FORMAT[(int)output_output_format]<<"\n";
            for(int i=0;i!=6;i++){
                out<<"section["<<i<<"]:\t"<<COMPRESSION_SECTION[(int)compression_section[i]]<<"\n";
            }

        }
    }
};


const char* program_config::VARS[20]={
		"in-file", "in-format", "clustering-mode",
		"clustering-test", "clustering-thres", "metric",
		"metric-class", "matrix-profile", "opt-query",
		"opt-instances", "opt-timeout", "out-keep-order",
		"out-gen-matrix", "out-format",
		"section[0]", "section[1]", "section[2]", "section[3]",
		"section[4]", "section[5]"
};

const char* program_config::IN_FORMAT[3]={"smile","fasta","fastq"};
const char* program_config::CLUSTERING_CLUSTERING[2]={"none","barcode"};
const char* program_config::CLUSTERING_TEST[2]={"chi-squared","chi-squared-c"};
const char* program_config::MATRIX_METRIC[4]={"delta","head-tail","k-head-tail","shift-hamming"};
const char* program_config::MATRIX_OPT_PROFILE[2]={"time","space"};
const char* program_config::BOOLEAN[2]={"false","true"};
const char* program_config::OUTPUT_OUTPUT_FORMAT[3]={"ascii","binary","html"};
const char* program_config::COMPRESSION_SECTION[2]={"none","deflate"};

int main(int argc, char* argv[]){
	cout<<"\n";
    program_config myconfig;
    
    if(argc==1){cerr<<"Syntax error. Use {-f filename} to load a configuration file or pass parameters by shell.\nUse {-?} to get help.\n";return -1;}
    else if(strcmp(argv[1],"-?")==0){
		cout<<"Sorry, I'm a lazy person.\nDo you seriously thought I would have spent precious time writing how to use this program? :P\n";
		return 0;
	}
    else if(strcmp(argv[1],"-f")==0){
        if(argc>=2){
            char buffer[360];
            char var[180], val[180];
            //Load config file
            cout<<"Reading file {"<<argv[2]<<"}\n";
           
            ifstream cfg(argv[2], ios::in);
            if(!cfg.is_open()){
            cerr<<"Warning: unable to open file {"<<argv[2]<<"}.\n";   
            }
            while(cfg.good()){
                cfg.getline(buffer,360);
               
                //Skip empty line!
                if(buffer[0]=='\0')continue;
               
                int j=0;
                for(;buffer[j]==' ' or buffer[j]=='\t';j++);
                for(;buffer[j]!='\0' && buffer[j]!=' ' && buffer[j]!='\t' && buffer[j]!=':' && j<180-1 ;j++){
                    var[j]=buffer[j];
                }
                var[j]=0;
                       
                //If parameter is inside continue!
                for(;buffer[j]==' ' or buffer[j]=='\t';j++);

                if(buffer[j]==':')j++;
               
                for(;buffer[j]==' ' or buffer[j]=='\t';j++);
                       
                int j2=0;
                for(;buffer[j]!='\0' && j2<180-1 ;j++,j2++){
                    val[j2]=buffer[j];
                }       
                val[j2]=0;
                       
                        //cout<<var<<": "<<val<<"\n";
                        myconfig.set(var,val);
               
            }
           
            cfg.close();
            //Skip visited options and read from command line!
            argc-=2;
            argv+=2;
        }
        else{
            cerr<<"Syntax error: name required after file request!";
            return -1;
        }
    }
   
    char var[180], val[180];
   
    //Continue reading from command line
    for(int i=1;i<argc;i++){
        //Set variable in myconfig

       
        if(argv[i][0]!='-' and argv[i][1]!='-'){
            cerr<<"Warning: parameter "<<i<<" has a syntax error! {"<<argv[i]<<"}\n";
            continue;
        }
        int j=2;
        for(;argv[i][j]!='\0' && argv[i][j]!='=' && j<180+2-1 ;j++){
            var[j-2]=argv[i][j];
        }
        var[j-2]=0;
       
        //If parameter is inside continue!
        if(argv[i][j]=='=')j++;
       
        int j2=0;
        for(;argv[i][j]!='\0' && j2<180+2-1 ;j++,j2++){
            val[j2]=argv[i][j];
        }       
        val[j2]=0;
       
        //cout<<var<<": "<<val<<"\n";
        myconfig.set(var,val);
    }
   
    //Now I can start with the main program!
    
    cout<<"~~~ Smile will use this final configuration:\n";
    myconfig.dump(cout,0);
    
    fstream fprofile(myconfig.input_file+"-profile.txt",ios::out|ios::binary);
    myconfig.dump(fprofile,0);
    fprofile.close();

	int NR=0;
	int LR=0;
	char *reads=nullptr;
	size_t *opt=nullptr;
	size_t *opt_best=nullptr;
	int best_len=0;

	fstream a0,a1,a2,a3,a4;
	
	try{
	
    cout<<"\n~~~ Opening reads from {"<<myconfig.input_file+".reads.txt"<<"}\n";
	fstream f(myconfig.input_file+".reads.txt",ios::in);
	if(!f.is_open()){
		cerr<<"Unable to open the file {"<<myconfig.input_file+".reads.txt"<<"}. Smile will be closed.\n";
		return 1;
	}
	
	if(myconfig.input_format!=program_config::SMILE){
		cerr<<"Sorry but the input format {"<<program_config::IN_FORMAT[myconfig.input_format]<<"} is still not supported!\n";
		cerr<<"Smile will be closed.\n";
		return 1;
	};
	
	if(myconfig.input_format==program_config::SMILE){
		f>>NR;
		f>>LR;
		
	    cout<<"\n~~~ Header found!\n";
	    cout<<"    Reads number:\t"<<NR<<"\n";
	    cout<<"    Reads length:\t"<<LR<<"\n";	
	    cout<<"\n~~~ Reserving memory\n";
	    reads=new char[NR*(LR+1)];
	    cout<<"\n~~~ Reading reads\n    ";
		for(int i=0;i<NR;i++){f>>(char*)&reads[i*(LR+1)];reads[(i+1)*(LR+1)-1]=0;if(i%100==0)cout<<".";};
	    cout<<"\n~~~ All the reads have been copied\n";

	}	
	
	f.close();
	
	if(myconfig.clustering_clustering!=program_config::NONE){
		cerr<<"Sorry but the clustering based on {"<<program_config::CLUSTERING_CLUSTERING[myconfig.input_format]<<"} is still not supported!\n";
		cerr<<"Smile will be continue working skipping this rule.\n";
	};	

	if(myconfig.matrix_metric==program_config::HEAD_TAIL){
		fast_dna_read<uint64_t> mask(LR);
		mask.create_mask();
		fast_dna_read<uint64_t>::current_mask=&mask;
		
		typedef distance_matrix<fast_dna_read<uint64_t>,fast_dna_read<uint64_t>::head_tail_struct>  dmatrix;
		dmatrix	MATRICE(NR, LR, (const uint8_t*)reads, fast_dna_read<uint64_t>::head_tail, false);
		
		cout<<"\n~~~ Parsing query\n    "<<myconfig.optimizer_query<<"\n";
		linear_optimizer<dmatrix, int> optimizer;
		
		if(myconfig.optimizer_query=="GREEDY")optimizer.QUERY.parse(query::GREEDY.c_str());
		else if(myconfig.optimizer_query=="GREEDY_STAR")optimizer.QUERY.parse(query::GREEDY_STAR.c_str());
		else if(myconfig.optimizer_query=="GREEDY_N")optimizer.QUERY.parse(query::GREEDY_N.c_str());
		else if(myconfig.optimizer_query=="GEN_AND_TEST")optimizer.QUERY.parse(query::GEN_AND_TEST.c_str());
		else if(myconfig.optimizer_query=="PARTICLE")optimizer.QUERY.parse(query::PARTICLE.c_str());
		else if(optimizer.QUERY.parse(myconfig.optimizer_query.c_str())!=0){
			cerr<<"Parsing error in the query.\nSmile will be closed.\n";
			return 1;
		}
				
		opt=new size_t[NR];
		opt_best=new size_t[NR];
		best_len=NR*LR;
		
		for(int i=0;i!=NR;i++)opt_best[i]=i;

		cout<<"\n~~~ Evaluating distance table. It may take a while\n";
		MATRICE.eval_matrix();
		cout<<"\n~~~ Evaluation completed!\n";

		cout<<"\n~~~ Time limited magic begins: "<<myconfig.optimizer_timeout<<" per task\n";
		for(int magic_index=0;magic_index<myconfig.optimizer_instances;magic_index++){
			cout<<"\n~~~ Magic ("<<magic_index+1<<"|"<<myconfig.optimizer_instances<<") is happening. It may take a long!\n";
			
			uint64_t ms=myconfig.optimizer_timeout.to_us()/1000;
			atomic<bool> closed(1);
			
			try{
				std::thread th([&](){try{optimizer(closed, MATRICE,opt,0,best_len);}catch(const char* e){cerr<<"Error: "<<e;}closed=0;});	
				th.detach();
				
				for(uint64_t i=0;i<ms;i++){usleep(1000);if(closed==0)break;}
				closed=0;
				usleep(5000);
				
			}
			catch(...){}
			
			if(closed!=0){cout<<"\n    Thread fault!\n";}
			else{cout<<"\n    Thread has correctly ended\n";}

			
			if(best_len>optimizer.result()){
				cout<<"    New better result: "<<optimizer.result()<<"\n";
				best_len=optimizer.result();
				memcpy(opt_best,opt,NR*sizeof(size_t));
			}
		}
		
		delete[] opt;
		
		cout<<"\n~~~ Magic ended. Best result reached is: "<<best_len<<"\n";

		a0.open(myconfig.input_file+"-a0.txt",ios::out|ios::binary);
		a1.open(myconfig.input_file+"-a1.txt",ios::out|ios::binary);
		a2.open(myconfig.input_file+"-a2.txt",ios::out|ios::binary);
		a3.open(myconfig.input_file+"-a3.txt",ios::out|ios::binary);

		if(!a0.is_open()){cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-a0.txt"<<"}\nSmile will be closed.\n";return 1;}		
		if(!a1.is_open()){cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-a1.txt"<<"}\nSmile will be closed.\n";return 1;}
		if(!a2.is_open()){cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-a2.txt"<<"}\nSmile will be closed.\n";return 1;}
		if(!a3.is_open()){cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-a3.txt"<<"}\nSmile will be closed.\n";return 1;}

		intersect_compression<dmatrix> compressor(MATRICE,opt_best,a1,a2,a3);
		
		cout<<"\n~~~ Precompression started. It may take a while\n";
		compressor.compress();
		
		//ADD a0 information!
		{
			bit_ofstream<fstream> header(a0);
			header.write(NR,64);
			header.write(LR,64);
			header.write(0 ,64);
		}
		
		cout<<"\n~~~ Precompression ended\n";
		
		a0.close();
		a1.close();
		a2.close();
		a3.close();

		if(myconfig.output_keep_order==true){
			cout<<"\n~~~ Exporting the original order vector\n";
			
			a4.open(myconfig.input_file+"-a4.txt",ios::out|ios::binary);
			if(!a4.is_open()){cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-a4.txt"<<"}\nSmile will be closed.\n";return 1;}
			
			bit_ofstream<fstream> swappings(a4);
			int bits=log2(NR);
			for(int i=0;i!=NR;i++)swappings.write(opt_best[i],bits);	//Consider max 2^32 reads
			a4.close();

			cout<<"\n~~~ Order vector exported\n";
		}
		
		delete[] opt_best;
	
		if(myconfig.output_generate_matrix==true){
			cout<<"\n~~~ Exporting the distance table\n";

			cerr<<"Warning, the generation of a human-readable table is still incomplete!\n";
			fstream ftable(myconfig.input_file+"-table.html",ios::out);
			if(ftable.is_open())MATRICE.dump(ftable, (dmatrix::dump_format)myconfig.output_output_format);
			else{cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-table.html"<<"}\n";}
			ftable.close();
			
			cout<<"\n~~~ Distance table exported\n";
		}

	}
	else if(myconfig.matrix_metric==program_config::SHIFT_HAMMING){
		typedef distance_matrix<fast_dna_read<uint64_t>,fast_dna_read<uint64_t>::shift_hamming_struct>  dmatrix;
		dmatrix	MATRICE(NR, LR, (const uint8_t*)reads, fast_dna_read<uint64_t>::shift_hamming, false);
		
		cout<<"\n~~~ Parsing query\n    "<<myconfig.optimizer_query<<"\n";
		linear_optimizer<dmatrix, int> optimizer;
		
		if(myconfig.optimizer_query=="GREEDY")optimizer.QUERY.parse(query::GREEDY.c_str());
		else if(myconfig.optimizer_query=="GREEDY_STAR")optimizer.QUERY.parse(query::GREEDY_STAR.c_str());
		else if(myconfig.optimizer_query=="GREEDY_N")optimizer.QUERY.parse(query::GREEDY_N.c_str());
		else if(myconfig.optimizer_query=="GEN_AND_TEST")optimizer.QUERY.parse(query::GEN_AND_TEST.c_str());
		else if(myconfig.optimizer_query=="PARTICLE")optimizer.QUERY.parse(query::PARTICLE.c_str());
		else if(optimizer.QUERY.parse(myconfig.optimizer_query.c_str())!=0){
			cerr<<"Parsing error in the query.\nSmile will be closed.\n";
			return 1;
		}
				
		opt=new size_t[NR];
		opt_best=new size_t[NR];
		best_len=NR*LR;
		
		for(int i=0;i!=NR;i++)opt_best[i]=i;

		cout<<"\n~~~ Evaluating distance table. It may take a while\n";
		MATRICE.eval_matrix();
		cout<<"\n~~~ Evaluation completed!\n";

		cout<<"\n~~~ Time limited magic begins: "<<myconfig.optimizer_timeout<<" per task\n";
		for(int magic_index=0;magic_index<myconfig.optimizer_instances;magic_index++){
			cout<<"\n~~~ Magic ("<<magic_index+1<<"|"<<myconfig.optimizer_instances<<") is happening. It may take a long!\n";
			
			uint64_t ms=myconfig.optimizer_timeout.to_us()/1000;
			atomic<bool> closed(1);
			
			try{
				std::thread th([&](){try{optimizer(closed, MATRICE,opt,0,best_len);}catch(const char* e){cerr<<"Error: "<<e;}closed=0;});	
				th.detach();
				
				for(uint64_t i=0;i<ms;i++){usleep(1000);if(closed==0)break;}
				closed=0;
				usleep(5000);
				
			}
			catch(...){}
			
			if(closed!=0){cout<<"\n    Thread fault!\n";}
			else{cout<<"\n    Thread has correctly ended\n";}

			
			if(best_len>optimizer.result()){
				cout<<"    New better result: "<<optimizer.result()<<"\n";
				best_len=optimizer.result();
				memcpy(opt_best,opt,NR*sizeof(size_t));
			}
		}
		
		delete[] opt;
		
		cout<<"\n~~~ Magic ended. Best result reached is: "<<best_len<<"\n";

		a0.open(myconfig.input_file+"-a0.txt",ios::out|ios::binary);
		a1.open(myconfig.input_file+"-a1.txt",ios::out|ios::binary);
		a2.open(myconfig.input_file+"-a2.txt",ios::out|ios::binary);
		a3.open(myconfig.input_file+"-a3.txt",ios::out|ios::binary);

		if(!a0.is_open()){cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-a0.txt"<<"}\nSmile will be closed.\n";return 1;}		
		if(!a1.is_open()){cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-a1.txt"<<"}\nSmile will be closed.\n";return 1;}
		if(!a2.is_open()){cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-a2.txt"<<"}\nSmile will be closed.\n";return 1;}
		if(!a3.is_open()){cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-a3.txt"<<"}\nSmile will be closed.\n";return 1;}

		strip_compression<dmatrix> compressor(MATRICE,opt_best,a1,a2,a3);
		
		cout<<"\n~~~ Precompression started. It may take a while\n";
		compressor.compress();
		
		//ADD a0 information!
		{
			bit_ofstream<fstream> header(a0);
			header.write(NR,64);
			header.write(LR,64);
			header.write(1 ,64);
		}
		
		cout<<"\n~~~ Precompression ended\n";
		
		a0.close();
		a1.close();
		a2.close();
		a3.close();

		if(myconfig.output_keep_order==true){
			cout<<"\n~~~ Exporting the original order vector\n";
			
			a4.open(myconfig.input_file+"-a4.txt",ios::out|ios::binary);
			if(!a4.is_open()){cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-a4.txt"<<"}\nSmile will be closed.\n";return 1;}
			
			bit_ofstream<fstream> swappings(a4);
			int bits=log2(NR);
			for(int i=0;i!=NR;i++)swappings.write(opt_best[i],bits);	//Consider max 2^32 reads
			a4.close();

			cout<<"\n~~~ Order vector exported\n";
		}
		
		delete[] opt_best;
	
		if(myconfig.output_generate_matrix==true){
			cout<<"\n~~~ Exporting the distance table\n";

			cerr<<"Warning, the generation of a human-readable table is still incomplete!\n";
			fstream ftable(myconfig.input_file+"-table.html",ios::out);
			if(ftable.is_open())MATRICE.dump(ftable, (dmatrix::dump_format)myconfig.output_output_format);
			else{cerr<<"Warning, unable to use file {"<<myconfig.input_file+"-table.html"<<"}\n";}
			ftable.close();
			
			cout<<"\n~~~ Distance table exported\n";
		}		
	}
	else{
		cerr<<"Sorry but the metric {"<<program_config::MATRIX_METRIC[myconfig.matrix_metric]<<"} is still not supported!\n";
		cerr<<"Smile will be closed.\n";
		return 1;
	}

	cout<<"\n~~~ Operations completed. Smile will stop!\n";
	}
	catch(const char* str){
		cerr<<"Excepion: {"<<str<<"}\nSmile will be closed.";
		delete[] reads;
		delete[] opt_best;
		delete[] opt;

		if(a0.is_open()){a0.close();}		
		if(a1.is_open()){a1.close();}		
		if(a2.is_open()){a2.close();}		
		if(a3.is_open()){a3.close();}		
		if(a4.is_open()){a4.close();}		

		return 1;
	}
	
	delete[] reads;

    return 0;
}
