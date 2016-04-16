/* This code is shared with no explicit or implicit guarantee about 
 * functionalities. It's part of a thesis. For more information:
 * http://pikkolamakkia.com/thesis
 * 
 * Author: Carlo Masaia
 * Version: 1.0-thesis
 * Last update: 2016.Feb.11
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
#include <cmath>

#include "bitstream.h"

using namespace std;


template<typename DISTANCE_MATRIX>
class strip_compression{
	private:
		const size_t MAX_STRIP=256;
		
		size_t* optimal_sequence;
		DISTANCE_MATRIX& matrix;
		
		ostream& shift;
		ostream& symbols;
		ostream& strips;
		
	public:
		strip_compression(	DISTANCE_MATRIX& m, size_t* opt,
							ostream& pshift, ostream& psymbols, ostream& pstrips)
							:matrix(m),shift(pshift),symbols(psymbols),strips(pstrips){optimal_sequence=opt;}
	
		//It DOES work!
		void compress(){
			
			bit_ofstream<ostream> shift(this->shift);
			bit_ofstream<ostream> symbols(this->symbols);
			bit_ofstream<ostream> strips(this->strips);			
				
			size_t *positions=new size_t[matrix.read_len()];
			size_t *offset=new size_t[matrix.reads()];
			
			for(size_t i=0;i!=matrix.read_len();i++)positions[i]=0;
			
			//Build up the absolute offset sequence
			size_t toffset=0;
			for(size_t i=0;i!=matrix.reads();i++){
				offset[i]=(toffset=(toffset+matrix.distance(optimal_sequence[i],optimal_sequence[(i+matrix.reads()-1)%matrix.reads()]).pos)%matrix.read_len());
			}
			
			
			unsigned long long int counted=0;
			
			for(size_t column=matrix.read_len()-1;
				counted<matrix.reads()*matrix.read_len();){
					
					size_t c_column=(column+1)%matrix.read_len();

					size_t strip=0;
					//Consider the first base still not yet evaluated at position c_column
					unsigned int base=matrix.get_read(optimal_sequence[positions[c_column]])[(c_column+matrix.read_len()-offset[positions[c_column]])%matrix.read_len()];

					//Increase the strip length! If the strip is too big, just split it!
					for(;positions[c_column]+strip<matrix.reads() and 
						 matrix.get_read(optimal_sequence[positions[c_column]+strip])[(c_column+matrix.read_len()-offset[positions[c_column]+strip])%matrix.read_len()]==base and
						 strip<MAX_STRIP;
						 strip++);
						 
					//Increase the position where to start next time!
					
					positions[c_column]+=strip;
					counted+=strip;
										
					//Save strip len and base in the streams!
					symbols.write(base,2);
					//strips<<strip<<";";
					strips.write(strip,8);	//log2(256)
					
					//Next column if condition in reached or !
					if(positions[c_column]>positions[column] or positions[c_column]>=matrix.reads())
						column=c_column;
			}
			
			delete[] positions;
			delete[] offset;
			
			int bits=log2(matrix.read_len());
			//Save the shifts table
			for(unsigned int i=1;i<matrix.reads();i++){
				auto data=matrix.distance(optimal_sequence[i],optimal_sequence[i-1]).pos;
				
				shift.write(data,bits);	//reads are usually small!
			}
		}

};


template<typename DISTANCE_MATRIX>
class intersect_compression{
	public:
		intersect_compression(	DISTANCE_MATRIX& m, size_t* opt,
								ostream& perrors, ostream& psymbols, ostream& ppointers)
								:matrix(m),errors(perrors),symbols(psymbols),pointers(ppointers){optimal_sequence=opt;}
	
		//I do not evaluate errors in this implementation, since the metric used for this example doesn't
		void compress(){
			unsigned long long int last_pointer=0;
			unsigned long long int pos=0;
			
			bit_ofstream<ostream> bsymbols(this->symbols);
			bit_ofstream<ostream> bpointers(this->pointers);


			int log_size=log2(matrix.reads()*matrix.read_len());
			for(size_t i=1;i<matrix.reads();i++){
				bpointers.write(pos-last_pointer,log_size);
				
				size_t distance=matrix(optimal_sequence[i],optimal_sequence[i-1]);
								
				for(size_t j=0;j<distance;j++){
					
					int k = matrix.get_read(optimal_sequence[i])[j];
					/*if(k==0)symbols<<'A';
					else if(k==1)symbols<<'C';
					else if(k==2)symbols<<'G';
					else symbols<<'T';*/
					bsymbols.write(k,2);
				}
				
				last_pointer=pos;
				pos+=distance;
			}
			
		}
		
	private:
		size_t* optimal_sequence;
		DISTANCE_MATRIX& matrix;
		
		
		
		ostream& errors;
		ostream& symbols;
		ostream& pointers;
};
