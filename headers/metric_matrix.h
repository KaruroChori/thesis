/* This code is shared with no explicit or implicit guarantee about 
 * functionalities. It's part of a thesis. For more information:
 * http://pikkolamakkia.com/thesis
 * 
 * Author: Carlo Masaia
 * Version: 1.0-thesis
 * Last update: 2016.Feb.12
 * 
 * The rights of the authors on this code are enforced by Italian and 
 * European laws on copyright. Every kind of usage, alteration, copy, 
 * printing of what it follows is strictly forbidden without the explicit
 * permission of the authors.
*/

#pragma once

#include <iostream>
#include <cstring>
#include <fstream>
#include <vector>
#include <algorithm>
#include <random>
#include <chrono>

#include "generic_distance.h"
#include "set.h"
#include "swapper.h"
#include "sorting.h"

using namespace std;

template<typename READ, typename DISTANCE, typename SCALAR=int>
struct distance_matrix{
	public:
		
		//Build from file or memory text array
		
		distance_matrix(size_t n_reads, size_t len_reads, ostream& data, DISTANCE& dt, bool binary=false):dist_type(dt){
			
		}
		
		distance_matrix(size_t n_reads, size_t len_reads, const uint8_t* data, DISTANCE& dt, bool binary=false)
						:mask(n_reads),dist_type(dt){
				
			if(binary==true){	//Binary
				
				this->n_reads=n_reads;
				this->len_reads=len_reads;
				
				void* raw_memory = operator new[]( n_reads * sizeof( READ ) );
				READ* ptr = static_cast<READ*>( raw_memory );
				for( size_t i = 0; i < n_reads; ++i ) {
					new( &ptr[i] )READ( len_reads, &data[READ::kappas(len_reads)*i], true);
				}
				read_data=raw_memory;
				mask=-mask;
			}
			else{				//Text

				this->n_reads=n_reads;
				this->len_reads=len_reads;
				
				void* raw_memory = operator new[]( n_reads * sizeof( READ ) );
				READ* ptr = static_cast<READ*>( raw_memory );
				for( size_t i = 0; i < n_reads; ++i ) {
					new( &ptr[i] )READ( len_reads, &data[(len_reads+1)*i], false);
				}
				read_data=raw_memory;
				mask=-mask;
			}
			
			parent=this;
			
			distances=new typename DISTANCE::naked[n_reads*n_reads];
			//Init swappers!
			{
			void* raw_memory = operator new[]( n_reads * sizeof( swapper<size_t> ) );
			swapper<size_t>* ptr = static_cast<swapper<size_t>*>( raw_memory );
			for( size_t i = 0; i < n_reads; ++i ) {
				new( &ptr[i] )swapper<size_t>(n_reads);
			}
			rows_swp=(swapper<size_t>*)raw_memory;
			}			
		}
		
		//Reference to an other table:
		
		distance_matrix(const distance_matrix& ref):mask(ref.mask),dist_type(ref.dist_type){
			parent=&ref;
			distances=ref.distances;
			rows_swp=ref.rows_swp;
			read_data=ref.read_data;
			
			n_reads=ref.n_reads;
			len_reads=ref.len_reads;
		}
		
		/*USELESS
		distance_matrix(const distance_matrix& ref, const set<size_t>& rows):dist_type(ref.dist_type){
			parent=&ref;
			distances=ref.distances;
			rows_swp=ref.rows_swp;
			read_data=ref.read_data;
			
			n_reads=ref.n_reads;
			len_reads=ref.len_reads;
			
			if(rows<ref.mask)mask=rows;
			else throw "Row doesn't exist";			
		}*/
		
		~distance_matrix(){
						
			//Delete the table only if I'm the real owner!
			if(parent==this){
				delete[] distances;
				
				{
				swapper<size_t>* ptr = rows_swp;

				for( int i = n_reads - 1; i >= 0; --i ) {
					ptr[i].~swapper<size_t>();
				}
				operator delete[]( rows_swp );		
				}	
				
				{
				READ* ptr = static_cast<READ*>( read_data );

				for( int i = n_reads - 1; i >= 0; --i ) {
					ptr[i].~READ();
				}
				operator delete[]( read_data );		
				}	
			}
			
			distances=nullptr;
			rows_swp=nullptr;
			parent=nullptr;
			read_data=nullptr;
		}
		
		enum dump_format{RAW_TXT,BINARY,HTML,RTF};
		void dump(ostream& out, dump_format format=BINARY){
		
		//TODO: Check for symmetry!
		/*if(format==HTML){
			out<<"<table>";
			for(size_t i=0;i!=n_reads;i++){
				out<<"<tr>";
				
				out<<"<td colspan=\""<<i+1<<"\" align=\"right\">"<<distance(i,i)<<"</td>";
				for(size_t j=i+1;j<n_reads;j++){
					out<<"<td>"<<distance(i,j)<<"</td>";
				}
				out<<"</tr>";
			}
			out<<"</table>";
			
		}*/
		
		if(format==HTML){
			out<<"<table>";
			for(size_t i=0;i!=n_reads;i++){
				out<<"<tr>";
				
				for(size_t j=0;j<n_reads;j++){
					out<<"<td>"<<distance(i,j)<<"</td>";
				}
				out<<"</tr>";
			}
			out<<"</table>";
			
		}
		
		else if(format==RAW_TXT){
			out<<"#distance_table\n";
			for(size_t i=0;i!=n_reads;i++){
				for(size_t j=0;j!=n_reads;j++){
					cout<<distance(i,j)<<"\t";
				}
				out<<"\n";
			}
			out<<"\n";
			
		}

		else throw "Format still not implemented!";
			
		}
		
		struct query_struct{
			enum SET{IDENTITY,MIN,MAX}set;
			int parameter_1, parameter_2, threshold;
			int limit;
			enum SELECTOR{ANY,ASC_ORDER,DESC_ORDER,RAND,RAND_P,RAND_P_INV}selector;
			
			int parse(const char* str){
				char token[64];

				const char* s=str;
				for(;(*s==' ' or *s=='{' or *s=='\t' or *s=='\n') and *s!='0';s++);//Skip spaces
				if(memcmp(s,"select",strlen("select"))!=0)return -1;	//Ask for "select"
				s+=strlen("select");									//Skip "select"
				for(;(*s==' ' or *s=='\t' or *s=='\n') and *s!='0';s++);//Skip spaces
				
				int i=0;
				for(;s[i]!=' '  and
					 s[i]!='\t' and
					 s[i]!='\n' and
					 s[i]!='\0' and
					 s[i]!='('  and i<63;i++)token[i]=s[i];				//Copy untill spaces	
				token[i]=0;

				if(strcmp(token,"IDENTITY")==0){
					s+=strlen("IDENTITY");
					set=IDENTITY;
				}
				else if(strcmp(token,"MIN")==0){
					s+=strlen("MIN");
					set=MIN;
					parameter_1=0;parameter_2=1;
				}
				else if(strcmp(token,"MAX")==0){
					s+=strlen("MAX");
					set=MAX;
					parameter_1=0;parameter_2=1;
				}
				else return -1;
				
				for(;(*s==' ' or *s=='\t' or *s=='\n') and *s!='0';s++);//Skip spaces
				if(*s=='('){											//Save parameters
					s++;
					
					for(;(*s==' ' or *s=='\t' or *s=='\n') and *s!='0';s++);//Skip spaces
					for(i=0;s[i]>='0' and s[i]<='9';i++)token[i]=s[i];
					token[i]=0;
					s+=i;
					parameter_1=atoi(token);
					
					for(;(*s==' ' or *s=='\t' or *s=='\n' or *s==',') and *s!='0';s++);//Skip spaces
					for(i=0;s[i]>='0' and s[i]<='9';i++)token[i]=s[i];
					token[i]=0;
					s+=i;
					parameter_2=atoi(token);
					
					for(;(*s==' ' or *s=='\t' or *s=='\n') and *s!='0';s++);//Skip spaces
					
					if(*s!=')')return -1;
					else s++;		
				}
				
				for(;(*s==' ' or *s=='\t' or *s=='\n') and *s!='0';s++);//Skip spaces
				//Continue with "take"
				if(memcmp(s,"take",strlen("take"))!=0)return -1;		//Ask for "select"
				s+=strlen("take");									//Skip "select"
				for(;(*s==' ' or *s=='\t' or *s=='\n') and *s!='0';s++);//Skip spaces
				
				if(memcmp(s,"ALL",strlen("ALL"))==0){limit=-1;s+=strlen("ALL");}
				else{
					for(i=0;s[i]>='0' and s[i]<='9';i++)token[i]=s[i];
					token[i]=0;
					s+=i;
					limit=atoi(token);					
				}
				
				for(;(*s==' ' or *s=='\t' or *s=='\n') and *s!='0';s++);//Skip spaces
				
				if(memcmp(s,"with",strlen("with"))!=0)return -1;		//Ask for "select"
				s+=strlen("with");									//Skip "select"	

				for(;(*s==' ' or *s=='\t' or *s=='\n') and *s!='0';s++);//Skip spaces
				i=0;
				for(;s[i]!=' '  and
					 s[i]!='\t' and
					 s[i]!='\n' and
					 s[i]!='}' and
					 s[i]!='\0' and
					 s[i]!='('  and i<63;i++)token[i]=s[i];				//Copy untill spaces	
				token[i]=0;

				if(strcmp(token,"ANY")==0){
					s+=strlen("ANY");
					selector=ANY;
				}
				else if(strcmp(token,"ASC_ORDER")==0){
					s+=strlen("ASC_ORDER");
					selector=ASC_ORDER;
				}
				else if(strcmp(token,"DESC_ORDER")==0){
					s+=strlen("DESC_ORDER");
					selector=DESC_ORDER;
				}
				else if(strcmp(token,"RAND")==0){
					s+=strlen("RAND");
					selector=RAND;
				}
				else if(strcmp(token,"RAND_P")==0){
					s+=strlen("RAND_P");
					selector=RAND_P;
				}
				else if(strcmp(token,"RAND_P_INV")==0){
					s+=strlen("RAND_P_INV");
					selector=RAND_P_INV;
				}
				else return -1;	
				
				return 0;						
			}
		};
		
		//TODO: test MIN, copy back to MAX
		vector<size_t> query(query_struct& q, int row){
			quick_sort(&distances[n_reads*row], n_reads, rows_swp[row]);

			vector<size_t> temp;
			
			//Create the vector of positions according to the requirments in the query
			if(q.set==query_struct::IDENTITY){//insert all on the vector
				temp.resize(mask.cardinality());
				int q=0;
				for(size_t i=0;i!=n_reads;i++)if(rows_swp[row][i]<mask){temp[q]=rows_swp[row][i];q++;}
			}
			else if(q.set==query_struct::MIN){
				//Start from the bottom. Avoid the objects with p1 different distance values or stop when reads end up.
				//Keep going recording objects as long as they have p2-p1 different distances or reads ending up.
				
				size_t i=0;
				int slices=0;
				SCALAR old_min=0;
				
				//Look for the first value!
				for(;i<n_reads;){
					//When it is found save it in old_min and stop the cycle
					if(rows_swp[row][i]<mask){
						old_min=operator()(rows_swp[row][i],row);
						break;
					}
					//Else check the following!
					else i++;
				}
				
				//Skip q.parameter_1 slices at the beginning (if possible)
				for(;i<n_reads && slices<q.parameter_1;i++){
					//If the object is masked take the next one!
					if(!(rows_swp[row][i]<mask))continue;
					
					//If I'm changing value keep trace of the event
					if(operator()(row,rows_swp[row][i])>old_min){
						slices++;
						old_min=operator()(rows_swp[row][i],row);
					}
				}
				
				//Save the items inside the next q.parameter_2-q.parameter_1 slices.
				for(;i<n_reads && slices<q.parameter_2;i++){
					//If the object is masked take the next one!
					if(!(rows_swp[row][i]<mask))continue;
					
					if(operator()(rows_swp[row][i],row)>old_min){
						slices++;
						old_min=operator()(rows_swp[row][i],row);
						
						//What to do with my objecy?
						if(slices<q.parameter_2)temp.push_back(rows_swp[row][i]);
						else break;
					}
					else{
						temp.push_back(rows_swp[row][i]);
					}
				}
			
			}
			else if(q.set==query_struct::MAX){
				throw "Function not implemented";
			}
			else throw "Function not implemented";
			
			//Selection from the set according to selector
			
			if(q.limit==-1 or temp.size()==0)return temp;
			if((int)temp.size()<=q.limit)return temp;
			else if(q.selector==query_struct::ASC_ORDER or q.selector==query_struct::ANY){
				vector<size_t> t(q.limit);
				for(int i=0;i!=q.limit;i++)t[i]=temp[i];
				return t;
			}
			else if(q.selector==query_struct::DESC_ORDER){
				vector<size_t> t(q.limit);
				int pos=temp.size();
				for(int i=0;i!=q.limit;i++)t[i]=temp[pos-1-i];
				return t;
			}
			else if(q.selector==query_struct::RAND){
				vector<size_t> t(q.limit);
				
				unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
				default_random_engine generator(seed);
				
				uniform_int_distribution<int> distribution(0,temp.size()-1);
				for(int i=0;i!=q.limit;i++){t[i]=temp[distribution(generator)];}
				return t;
			}
			else if(q.selector==query_struct::RAND_P){
				vector<size_t> t(q.limit);
				
				vector<SCALAR> distmap(temp.size());
				SCALAR supersum=0;
				for(size_t i=0;i!=temp.size();i++){
					supersum+=distances[n_reads*row+temp[i]].lower_bound;
					distmap[i]=supersum;
				}
				
				unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
				default_random_engine generator(seed);
				uniform_int_distribution<SCALAR> distribution(0,supersum);
				
				for(int i=0;i!=q.limit;i++){
					int num=distribution(generator);
					int k=0;
					for(;distmap[k]<num;k++);
					t[i]=temp[k];
				}
				return t;
			}
			
			else if(q.selector==query_struct::RAND_P_INV){
				vector<size_t> t(q.limit);
				
				vector<SCALAR> distmap(temp.size());
				SCALAR supersum=0;
				for(size_t i=0;i!=temp.size();i++){
					supersum+=len_reads-distances[n_reads*row+temp[i]].lower_bound;
					distmap[i]=supersum;
				}
				
				unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
				default_random_engine generator(seed);
				uniform_int_distribution<SCALAR> distribution(0,supersum);
				
				for(int i=0;i!=q.limit;i++){
					int num=distribution(generator);
					int k=0;
					for(;distmap[k]<num;k++);
					t[i]=temp[k];

				}
				return t;
			}
			//ADD RANDOM BASED SELECTIONS!
			else {throw "Selector not implemented";}
			
		}
		
		enum fn_struct{
			NULL_FN,ITEMS,MEAN,MEDIAN,VARIANCE,SUM,MODE
		};
		
		double eval(fn_struct fn, int row){
			
			if(fn==NULL_FN)return 0;
			else if(fn==ITEMS){return mask.cardinality();}
			else if(fn==SUM){
				SCALAR partial=0;
				for(size_t i=0;i!=n_reads;i++)if(i<mask)partial+=operator()(row,i);
				return partial;
			}
			else if(fn==MEAN){
				SCALAR partial=0;int times=0;
				for(size_t i=0;i!=n_reads;i++)if(i<mask){partial+=operator()(row,i);times++;}
				return (double)partial/(double)times;
			}
			else if(fn==VARIANCE){
				SCALAR partial=0;SCALAR power=0;int times=0;
				for(size_t i=0;i!=n_reads;i++)if(i<mask){
					SCALAR t=operator()(row,i);
					partial+=t;
					power+=t*t;
					times++;
					}
				return (double)power/(double)times - (double)partial*(double)partial/(double)times/(double)times;
			}
			else throw "Function not implemented";
			
			return 0;
		}
		
		//Evaluate the object in a complete way
		inline SCALAR operator()(size_t x, size_t y){
			if(x>n_reads or y>n_reads){cerr<<"d"<<x<<" "<<y<<"\\\\\\n";throw "Out of Bounds {operator()}";}
			READ* ptr=(READ*)read_data;
			if(distances[y*n_reads+x].lower_bound!=distances[y*n_reads+x].upper_bound)distances[y*n_reads+x]=dist_type(ptr[x],ptr[y],0,len_reads);
			return distances[y*n_reads+x].lower_bound;
		}
		
		inline typename DISTANCE::naked distance(size_t x, size_t y){
			if(x>n_reads or y>n_reads)throw "Out of Bounds {distance(...)}";
			READ* ptr=(READ*)read_data;
			if(distances[y*n_reads+x].lower_bound!=distances[y*n_reads+x].upper_bound)distances[y*n_reads+x]=dist_type(ptr[x],ptr[y],0,len_reads);
			return distances[y*n_reads+x];
		}
		
		inline void eval_matrix(){
			READ* ptr=(READ*)read_data;
			for(size_t x=0;x<n_reads;x++)for(size_t y=0;y<n_reads;y++)
			distances[y*n_reads+x]=dist_type(ptr[x],ptr[y],0,len_reads);
		}
		
		inline READ& get_read(size_t pos){if(pos>n_reads)throw "Out of Bounds {get_read(...)}";READ* ptr=(READ*)read_data;return ptr[pos];}
		
		inline size_t reads() const {return n_reads;}
		inline size_t read_len() const{return len_reads;}
		
		set<size_t> mask;

	private:
	
		DISTANCE& dist_type;
		const distance_matrix* parent;
		
		swapper<size_t> *rows_swp;
		typename DISTANCE::naked *distances;
			
		size_t n_reads;
		size_t len_reads;
		
		void* read_data;
		
};
