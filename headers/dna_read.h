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
using std::ostream;

#include "generic_distance.h"
#include "countpop.h"

template<typename WORD>
struct fast_dna_read{
	
	enum BASE{A,C,G,T};
	
	fast_dna_read(size_t size){
		if(size<1)throw "Wrong Size";
		
		len=size;
		
		//Keep the space free as long as I do not need it!
		c_bases=nullptr;
		
		//Magic number used in order to avoid this computation in future
		kappa=	(len%(sizeof(WORD)*4)==0)?
				(len/(sizeof(WORD)*4)):(len/(sizeof(WORD)*4)+1);
				
		//I need to have enough space to save all the basic rotations
		bases=new bitfield[(kappa)*len]();
	}
	
	fast_dna_read(size_t size, const void* res, bool binary=false){
		if(size<1)throw "Wrong Size";
		len=size;
		c_bases=nullptr;
		kappa=	(len%(sizeof(WORD)*4)==0)?
				(len/(sizeof(WORD)*4)):(len/(sizeof(WORD)*4)+1);
		bases=new bitfield[(kappa)*len]();	
		
		if(binary==false){
			text_to_binary((const char*)res);
			rotations();
		}
		else{
			const WORD* res2=(const WORD*)res;
			for(int i=0;i<kappa;i++)bases[i].data=res2[i];
			rotations();	
		}
	}

	fast_dna_read(const fast_dna_read& copy){
		len=copy.len;
		kappa=copy.kappa;
		c_bases=nullptr;
		bases=new bitfield[(kappa)*len]();	

		for(size_t i=0;i<len*len;i++)
			bases[i].data=copy.bases[i].data;
	}
	
	~fast_dna_read(){delete[] c_bases;delete[] bases;len=0;}
	
	int read(const char* res){
		text_to_binary(res);
		rotations();
		return 0;
	}

	int read(const WORD* res){
		for(int i=0;i<kappa;i++)bases[i].data=res[i];
		rotations();
		return 0;
	}
	
	int read(const fast_dna_read& copy){
		if(copy.len!=len)return -1;
		for(size_t i=0;i<len*len;i++)
			bases[i].data=copy.bases[i].data;		
		return 0;
	}
	
	const char* write(){
		if(c_bases==nullptr){
			c_bases=new char[len+1];
		}
		
		binary_to_text();
		c_bases[len]=0;
		
		return c_bases;
	}
	
	static inline int kappas(size_t len){
		return 	(len%(sizeof(WORD)*4)==0)?
				(len/(sizeof(WORD)*4)):(len/(sizeof(WORD)*4)+1);
	}
	
	struct hamming_struct:generic_distance<int, const fast_dna_read&>{
		//virtual generic_distance<int, const fast_dna_read&>& operator()
		virtual hamming_struct& operator()
		(const fast_dna_read& a, const fast_dna_read& b, int inf=0, int sup=10000){
			if(a.len!=b.len)throw "Different size!";
			this->lower_bound=this->upper_bound=phamming(a.bases,b.bases,a.kappa);
			return *this;
		}
		virtual bool zero_implies_same() const{return true;}
		virtual bool same_implies_zero() const{return true;}
		virtual bool symmetric() const{return true;}
		virtual bool triangular_inequality() const{return true;}
		virtual bool non_negative() const{return true;}
	}; 
	static hamming_struct hamming;

	struct shift_hamming_struct:generic_distance<int, const fast_dna_read&>{
		
		//virtual generic_distance<int, const fast_dna_read&>& operator()
		virtual shift_hamming_struct& operator()
		(const fast_dna_read& wA, const fast_dna_read& wB, int inf=0, int sup=10000){
			if(wA.len!=wB.len){throw "Different size!";}

			this->upper_bound=(int)wA.len;
			this->pos=0;
			for(size_t i=0;i<wA.len;i++){
				int d=phamming(&wA.bases[i*wA.kappa],wB.bases,wA.kappa);
				if(d<this->upper_bound){
					//this->upper_bound=phamming(&wA.bases[i*wA.kappa],wB.bases,wA.kappa);
					this->upper_bound=d;
					this->pos=i;
				}
				if(d==0)break;
			}
			this->lower_bound=this->upper_bound;
			return *this;
		}
		virtual bool zero_implies_same() const{return true;}
		virtual bool same_implies_zero() const{return true;}
		virtual bool symmetric() const{return true;}
		virtual bool triangular_inequality() const{return true;}
		virtual bool non_negative() const{return true;}
		
		size_t pos;
		
		struct naked:public generic_distance<int, const fast_dna_read&>::naked{
			size_t pos;
			naked(){this->lower_bound=0;this->upper_bound=100000;pos=0;}
			naked(const shift_hamming_struct& r){this->lower_bound=r.lower_bound;this->upper_bound=r.upper_bound;pos=r.pos;}
			
			naked& operator=(const shift_hamming_struct& r){
				this->lower_bound=r.lower_bound;this->upper_bound=r.upper_bound;pos=r.pos;
				return *this;
			}
			
			friend std::ostream& operator<<(std::ostream& out, const naked& ref) {
				return ref.print(out);
			}
			
			std::ostream& print(std::ostream& out) const{
				if(this->lower_bound==this->upper_bound)out<<"("<<this->lower_bound<<", +"<<this->pos<<")" ;
				else out<<"(l:"<<this->lower_bound<<", h:"<<this->upper_bound<<")";
				return out;
			}
			
			inline friend bool operator>=(const naked&a, const naked&b){
				if(a.lower_bound>=b.upper_bound)return true;
				//else if(a.lower_bound>b.lower_bound)return true;
				else return false;
			}

			inline friend bool operator>(const naked&a, const naked&b){
				if(a.lower_bound>b.upper_bound)return true;
				//else if(a.lower_bound>b.lower_bound)return true;
				else return false;
			}

			inline friend bool operator<=(const naked&a, const naked&b){
				if(a.upper_bound<=b.lower_bound)return true;
				//else if(a.lower_bound>b.lower_bound)return true;
				else return false;
			}
			inline friend bool operator<(const naked&a, const naked&b){
				if((a.upper_bound)<(b.lower_bound))return true;
				//else if(a.lower_bound>b.lower_bound)return true;
				else return false;
			}
		};
						
		virtual void print(std::ostream& out) const{
			if(this->lower_bound==this->upper_bound)out<<this->lower_bound;
			else out<<"(l:"<<this->lower_bound<<", h:"<<this->upper_bound<<")";
		}

	}; 
	static shift_hamming_struct shift_hamming;

	//TODO:
	struct head_tail_struct:generic_distance<int, const fast_dna_read&>{
		
		virtual head_tail_struct& operator()
		(const fast_dna_read& wB, const fast_dna_read& wA, int inf=0, int sup=10000){
			if(wA.len!=wB.len){throw "Different size!";}

			this->upper_bound=std::min((int)wA.len,sup);
			int i=0;
			for(;i<this->upper_bound;i++){
				int d=phamming_masked(&wA.bases[i*wA.kappa],wB.bases,&(*current_mask).bases[i*wA.kappa],wA.kappa);
				if(d==0)break;
			}
			this->lower_bound=this->upper_bound=i;
			return *this;
		}
			
		virtual bool zero_implies_same() const{return true;}
		virtual bool same_implies_zero() const{return true;}
		virtual bool symmetric() const{return false;}
		virtual bool triangular_inequality() const{return true;}
		virtual bool non_negative() const{return true;}
		

		struct naked:public generic_distance<int, const fast_dna_read&>::naked{
			naked(){this->lower_bound=0;this->upper_bound=100000;}
			naked(const shift_hamming_struct& r){this->lower_bound=r.lower_bound;this->upper_bound=r.upper_bound;}
			
			naked& operator=(const head_tail_struct& r){
				this->lower_bound=r.lower_bound;this->upper_bound=r.upper_bound;
				return *this;
			}
			
			friend std::ostream& operator<<(std::ostream& out, const naked& ref) {
				return ref.print(out);
			}
			
			std::ostream& print(std::ostream& out) const{
				if(this->lower_bound==this->upper_bound)out<<this->lower_bound ;
				else out<<"(l:"<<this->lower_bound<<", h:"<<this->upper_bound<<")";
				return out;
			}
			
			inline friend bool operator>=(const naked&a, const naked&b){
				if(a.lower_bound>=b.upper_bound)return true;
				//else if(a.lower_bound>b.lower_bound)return true;
				else return false;
			}

			inline friend bool operator>(const naked&a, const naked&b){
				if(a.lower_bound>b.upper_bound)return true;
				//else if(a.lower_bound>b.lower_bound)return true;
				else return false;
			}

			inline friend bool operator<=(const naked&a, const naked&b){
				if(a.upper_bound<=b.lower_bound)return true;
				//else if(a.lower_bound>b.lower_bound)return true;
				else return false;
			}
			inline friend bool operator<(const naked&a, const naked&b){
				if((a.upper_bound)<(b.lower_bound))return true;
				//else if(a.lower_bound>b.lower_bound)return true;
				else return false;
			}
		};
						
		virtual void print(std::ostream& out) const{
			if(this->lower_bound==this->upper_bound)out<<this->lower_bound;
			else out<<"(l:"<<this->lower_bound<<", h:"<<this->upper_bound<<")";
		}
		
	};
	struct k_head_tail_struct:generic_distance<int, const fast_dna_read&>{
		virtual k_head_tail_struct& operator()
			(const fast_dna_read& a, const fast_dna_read& b, int inf=0, int sup=10000){return *this;}
			
		virtual bool zero_implies_same() const{return false;}
		virtual bool same_implies_zero() const{return true;}
		virtual bool symmetric() const{return false;}
		virtual bool triangular_inequality() const{return false;}
		virtual bool non_negative() const{return true;}
	};
	
	static head_tail_struct head_tail;
	static k_head_tail_struct k_head_tail;	
	
	
	void debug(ostream& out){
		for(size_t i=0;i<len;i++){
			//fast_dna_read t(len,(const WORD*)&bases[(kappa)*i])
			fast_dna_read t(kappa*sizeof(WORD)*4,(const WORD*)&bases[(kappa)*i],true);
			out<<i<<":\t"<<t.write()<<"\n";
		}	
	}

	inline BASE operator[](size_t pos) const{
		if(pos>=len or pos<0)throw "OutOfBounds";
		//else return bases[pos/(sizeof(WORD)*4/2)].get(pos%(sizeof(WORD)*4/2));
		else return bases[pos/(sizeof(WORD)*4)].get(pos%(sizeof(WORD)*4));
	}

	inline size_t length() const{return len;}

	private:
	
	//Unsafe structure on index i, internal use only!
	struct bitfield{
		WORD data;
		inline void set(BASE a, int i){
			data&= ~((WORD)3<<(i<<1));		//Erase a pair of bits, dec(3)=bin(11)
			data|=((WORD)a<<(i<<1));		//Set a pair of bits, a in set {00, 01, 10, 11}
		}
		inline BASE get(int i) const{
			//Mask two bits in a position and align them to be the first two bits.
			return (BASE)((data&((WORD)3<<(i<<1)))>>(i<<1));
		}
	};
	
	void text_to_binary(const char *wA){
		for(size_t i=0;i<len;i++){
			BASE B;
			switch(wA[i]){
				case 'A':
					B=BASE::A;
					break;
				case 'C':
					B=BASE::C;
					break;
				case 'G':
					B=BASE::G;
					break;
				case 'T':
					B=BASE::T;
					break;
				default:
					throw "Wrong Symbol";
					
			};
			bases[i/(sizeof(WORD)*4)].set(B,i%(sizeof(WORD)*4));
		}
	}
		
	void binary_to_text(){
		for(size_t i=0;i<len;i++){
			char B;
			switch(bases[i/(sizeof(WORD)*4)].get(i%(sizeof(WORD)*4))){
				case 0:
					B='A';
					break;
				case 1:
					B='C';
					break;
				case 2:
					B='G';
					break;
				case 3:
					B='T';
					break;
				default:
					throw "Wrong Symbol";
					
			};
			c_bases[i]=B;
		}
	}
	
	public:
	void create_mask(){
		for(int i=0;i!=kappa;i++)bases[i].data= ~(WORD)(0);
		
		bases[kappa-1].data&= ~(WORD)(0)>>(2*(WORD)(sizeof(WORD)*4-(len%(sizeof(WORD)*4))));	
			
		//for(size_t i=1;i<(sizeof(WORD)*4);i++){
		for(size_t i=1;i<len;i++){
			
			//Risolvo il caso base
			bases[(i)*kappa].data=bases[(i-1)*kappa].data<<(WORD)2;
			//bases[(i)*kappa].set(bases[(i-1)*kappa+kappa-1].get((len-1)%(sizeof(WORD)*4)),0);
			//[bases[(i)*kappa].set(bases[(i-1)*kappa+kappa-1].get((sizeof(WORD)*4)-1),0);
			
			for(int j=1;j<kappa;j++){
				bases[(i)*kappa+j].data=bases[(i-1)*kappa+j].data<<(WORD)2;
				bases[(i)*kappa+j].set(bases[(i-1)*kappa+j-1].get((sizeof(WORD)*4)-1),0);
			}
			
			//Mi assicuro che fuori dalla lunghezza della stringa ci sia solo 00 ovvero A.
			bases[(i)*(kappa)+kappa-1].data&= ~(WORD)(0)>>(2*(WORD)(sizeof(WORD)*4-(len%(sizeof(WORD)*4))));		
			
		}
	}
	
	private:
	
	void rotations(){
		
		bases[kappa-1].data&= ~(WORD)(0)>>(2*(WORD)(sizeof(WORD)*4-(len%(sizeof(WORD)*4))));	
			
		//for(size_t i=1;i<(sizeof(WORD)*4);i++){
		for(size_t i=1;i<len;i++){
			
			//Risolvo il caso base
			bases[(i)*kappa].data=bases[(i-1)*kappa].data<<(WORD)2;
			bases[(i)*kappa].set(bases[(i-1)*kappa+kappa-1].get((len-1)%(sizeof(WORD)*4)),0);
			//[bases[(i)*kappa].set(bases[(i-1)*kappa+kappa-1].get((sizeof(WORD)*4)-1),0);
			
			for(int j=1;j<kappa;j++){
				bases[(i)*kappa+j].data=bases[(i-1)*kappa+j].data<<(WORD)2;
				bases[(i)*kappa+j].set(bases[(i-1)*kappa+j-1].get((sizeof(WORD)*4)-1),0);
			}
			
			//Mi assicuro che fuori dalla lunghezza della stringa ci sia solo 00 ovvero A.
			bases[(i)*(kappa)+kappa-1].data&= ~(WORD)(0)>>(2*(WORD)(sizeof(WORD)*4-(len%(sizeof(WORD)*4))));		
			
		}
	}
	
	static int phamming(const bitfield* wA, const bitfield* wB, int kappa);
	static int phamming_masked(const bitfield* wA, const bitfield* wB, const bitfield* masks, int kappa);
	
	bitfield* 	bases;
	size_t 		len;
	char*		c_bases;
	int 		kappa;
	
	public:
	static fast_dna_read* current_mask;
}; 

template<typename WORD>
fast_dna_read<WORD> *fast_dna_read<WORD>::current_mask=nullptr;

template<typename WORD>
typename fast_dna_read<WORD>::hamming_struct fast_dna_read<WORD>::hamming;

template<typename WORD>
typename fast_dna_read<WORD>::shift_hamming_struct fast_dna_read<WORD>::shift_hamming;

template<typename WORD>
typename fast_dna_read<WORD>::head_tail_struct fast_dna_read<WORD>::head_tail;

template<typename WORD>
typename fast_dna_read<WORD>::k_head_tail_struct fast_dna_read<WORD>::k_head_tail;

template<>
int fast_dna_read<uint32_t>::phamming(const bitfield* wA, const bitfield* wB, int kappa){
	static uint32_t mask_1=0x55555555;
	static uint32_t mask_2=0xAAAAAAAA;
		
	int t=0;
	for(int i=0;i<kappa;i++){
		t+=__builtin_popcountll(
			((wA[i].data xor wB[i].data)&mask_1) |
			((wA[i].data xor wB[i].data)&mask_2)>>(uint32_t)1);
			
	}
	return t;	
}

template<>
int fast_dna_read<uint64_t>::phamming(const bitfield* wA, const bitfield* wB, int kappa){
	static uint64_t mask_1=0x5555555555555555;
	static uint64_t mask_2=0xAAAAAAAAAAAAAAAA;
		
	int t=0;
	for(int i=0;i<kappa;i++){
		t+=__builtin_popcountll(
			((wA[i].data xor wB[i].data)&mask_1) |
			((wA[i].data xor wB[i].data)&mask_2)>>(uint64_t)1);
			
	}
	return t;	
}


template<>
int fast_dna_read<uint64_t>::phamming_masked(const bitfield* wA, const bitfield* wB, const bitfield* mask, int kappa){
	static uint64_t mask_1=0x5555555555555555;
	static uint64_t mask_2=0xAAAAAAAAAAAAAAAA;
		
	int t=0;
	for(int i=0;i<kappa;i++){
		t+=__builtin_popcountll((
			((wA[i].data xor wB[i].data)&mask_1) |
			((wA[i].data xor wB[i].data)&mask_2)>>(uint64_t)1)&mask[i].data);
			
	}
	return t;	
}

template<>
int fast_dna_read<uint32_t>::phamming_masked(const bitfield* wA, const bitfield* wB, const bitfield* mask, int kappa){
	static uint32_t mask_1=0x55555555;
	static uint32_t mask_2=0xAAAAAAAA;
		
	int t=0;
	for(int i=0;i<kappa;i++){
		/*
		register uint32_t wwA=wA[i].data&mask[i].data;
		register uint32_t wwB=wB[i].data&mask[i].data;
		t+=__builtin_popcountll(
			((wwA xor wwB)&mask_1) |
			((wwA xor wwB)&mask_2)>>(uint32_t)1);
		*/
		t+=__builtin_popcountll((
			((wA[i].data xor wB[i].data)&mask_1) |
			((wA[i].data xor wB[i].data)&mask_2)>>(uint32_t)1)&mask[i].data);
			
	}
	return t;	
}
