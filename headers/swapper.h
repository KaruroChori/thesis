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

template<typename T=size_t>
struct swapper{
	swapper(){map=nullptr;}
	swapper(size_t len){
		this->len=len;
		map=new T[len];
		for(size_t i=0;i!=len;i++)map[i]=i;
	}
	
	~swapper(){
		delete[] map;
	}
	
	//Comment out the exception when not required any more!
	
	inline T operator[](size_t pos) const{
		if(pos>=len or pos<0) throw "OutOfBoundsException";
		return map[pos];
	}
	
	inline void swap(size_t a, size_t b){
		if(a>=len or a<0) throw "OutOfBoundsException";
		if(b>=len or b<0) throw "OutOfBoundsException";		
		T t=map[a];
		map[a]=map[b];
		map[b]=t;
	}
	
	inline void reset(){
		for(size_t i=0;i!=len;i++)map[i]=i;		
	}
	
	private:
		T *map;
		size_t len;
};
