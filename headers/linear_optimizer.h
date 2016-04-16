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
#include <atomic>

using namespace std;

template<typename DISTANCE_MATRIX, typename DISTANCE>
class linear_optimizer{
	private:
	DISTANCE best;
	
	//Return a sequence of indices to be followed building the optimal or suboptimal string
	virtual DISTANCE optimize(atomic<bool>& kill, DISTANCE_MATRIX matrix, size_t hint, 
								size_t* seq, size_t* inc_seq,
								DISTANCE c_dist, DISTANCE& inc_dist,
								int depth=0){
									
		if(kill==0)throw "TimeOut";
		
		//cout<<"Instance "<<c_dist<<" "<<inc_dist<<" "<<depth<<"\n";
		//Hide the hint from the matrix since it's going to be visited
		matrix.mask-=hint;

		//Add the hint to the list/stack
		seq[depth]=hint;

		//Check if distance reached isn't too much! Otherwise just cut off the branch
		if(c_dist>=inc_dist){cout<<"        Worse solution: "<<c_dist<<" at "<<depth<<"\n";return c_dist;}
	
		//Check if I reached the maximum in depth
		if((int)matrix.reads()==depth+1){
			//Compare my distance with the min. If it's a new incubent just copy it!
			if(c_dist<inc_dist){memcpy(inc_seq,seq,sizeof(size_t)*matrix.reads());inc_dist=c_dist;}
			cout<<"        New good solution: "<<c_dist<<"\n";			
			//Anyway return back;
			return c_dist;
		}								
									
		//Else continue visiting the nodes
		
		//Ask the matrix for the query
		vector<size_t> indeces=matrix.query(QUERY,hint);
		
		//TODO:
		//SE È VUOTO SIGNIFICA CHE QUESTO RAMO È BLOCCATO! PREVEDI MECCANISMO DI RIMOZIONE!

		//Test for all the possible indeces		
		for_each(indeces.begin(), indeces.end(), [&](size_t index){
			//Save the distance of the path plus the one between hint and index
			optimize(kill,matrix,index,seq,inc_seq,c_dist+matrix(index,hint),inc_dist,depth+1);
			if(kill==0)throw "TimeOut";
		});
	
		return c_dist;
	}
	
	public:
	
	typename DISTANCE_MATRIX::query_struct QUERY;

	
	//Wrapper for optimize
	inline DISTANCE operator()(atomic<bool> &kill, DISTANCE_MATRIX& matrix, size_t* min_seq, size_t hint=0, DISTANCE incubent=0){
		size_t *seq=new size_t[matrix.reads()];
		
		if(incubent==0)best=matrix.reads()*matrix.read_len();
		else best=incubent;
		
		DISTANCE t;
		
		try{
		t=optimize(kill,matrix,hint,seq,min_seq,0,best,0);
		}
		catch(const char* c){cerr<<"    Exception: Thread is going to be closed for {"<<c<<"}\n";}		
		delete[] seq;
		
		return t;
	}
	
	inline DISTANCE result() const{return best;}
	
	linear_optimizer(){best=0;}
};
