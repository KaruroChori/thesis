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

#include <ostream>

//In order to speed it up it is non-safe in checking the same length
//Tutte le funzioni sono testate eccetto quelle di comparazione per le quali mancano ancora dei tests
template<typename T>
struct set{
	public:
		typedef uint64_t WORD;

		set(size_t items){
			len=items;
			kappa=	(items%(sizeof(WORD)*4)==0)?
					(items/(sizeof(WORD)*4)):(items/(sizeof(WORD)*4)+1);
			data=new WORD[kappa];
			
			memset(data,0,kappa*sizeof(WORD));
		}

		set(const set& ref){
			kappa=ref.kappa;
			len=ref.len;
			data=new WORD[kappa];
			
			memcpy(data,ref.data,kappa*sizeof(WORD));
		}
		
		~set(){delete[] data;}

		inline set& operator=(const set& b){
			if(len!=b.len)throw "Exceptional len!";
			memcpy(data,b.data,kappa*sizeof(WORD));
			return *this;
		}
		
		friend inline set& operator+=(set& a, T b){
			size_t wb=(size_t)b;
			if(wb>a.len)throw "Exceptional len!";
			
			a.data[wb/(sizeof(WORD)*8)]|=((WORD)1)<<(wb%(sizeof(WORD)*8));
			 
			return a;
		}
		
		friend inline set& operator-=(set& a, T b){
			size_t wb=(size_t)b;
			if(wb>a.len)return a;
			
			a.data[wb/(sizeof(WORD)*8)]&=~(((WORD)1)<<(wb%(sizeof(WORD)*8)));
			 
			return a;
		}
		
		friend inline set& operator+=(set& a, const set& b){
			size_t limit=std::min(a.kappa,b.kappa);
			for(size_t i=0;i<limit;i++)a.data[i]|=b.data[i];
			return a;
		}
		
		friend inline set& operator*=(set& a, const set& b){
			size_t limit=std::min(a.kappa,b.kappa);
			for(size_t i=0;i<limit;i++)a.data[i]&=b.data[i];
			return a;
		}
		
		friend inline set& operator-=(set& a, const set& b){return a*=-b;}

		friend inline set operator+(const set& a, const set& b){set t(a);t+=b;return t;}
		friend inline set operator*(const set& a, const set& b){set t(a);t*=b;return t;}
		friend inline set operator-(const set& a, const set& b){set t(a);t-=b;return t;}

		friend inline set operator+(const set& a, T b){set t(a);t+=b;return t;}
		friend inline set operator-(const set& a, T b){set t(a);t-=b;return t;}
		
		friend inline set operator-(const set& a){
			set t(a);
			for(size_t i=0;i<t.kappa;i++)t.data[i]=~a.data[i];
			t.data[t.kappa-1]&=~((WORD)0)>>(WORD)(2*(sizeof(WORD)*8-t.len%(sizeof(WORD)*8)));
			return t;
		}

		friend inline bool operator==(const set& a, const set&b){
			size_t limit=std::min(a.kappa,b.kappa);
			for(size_t i=0;i<limit;i++)if(a.data[i]!=b.data[i])return false;
			return true;
		}
		friend inline bool operator!=(const set& a, const set&b){return !(a==b);}
		
		friend inline bool operator<(const set& a, const set&b){
			return (a<=b)&&(a!=b);
		}
		friend inline bool operator<=(const set& a, const set&b){
			size_t limit=std::min(a.kappa,b.kappa);
			for(size_t i=0;i<limit;i++)if((a.data[i]|b.data[i])!=b.data[i])return false;
			return true;
		}
		
		friend inline bool operator>(const set& a, const set&b){
			return (a>=b)&&(a!=b);			
		}
		
		friend inline bool operator>=(const set& a, const set&b){
			size_t limit=std::min(a.kappa,b.kappa);
			for(size_t i=0;i<limit;i++)if((a.data[i]|b.data[i])!=a.data[i])return false;
			return true;	
		}
		
		friend inline bool operator>(const set& a, T b){
			size_t wb=(size_t)b;
			if(wb>a.len)return false;
			return ((a.data[wb/(sizeof(WORD)*8)]>>(wb%(sizeof(WORD)*8)))&(WORD)1)==1;
		}
		friend inline bool operator<(T b, const set& a){return a>b;}

		void dump(std::ostream& out){
			out<<"{";
			for(size_t i=0;i<len;i++)if(*this>i)out<<(T)i<<",";
			out<<"}";
		}
		
		friend std::ostream& operator<<(std::ostream& out, const set & d){
			d.dump(out);
			return out;
		}		
		
		inline int cardinality(){
			int ret=0;
			for(size_t i=0;i<len;i++){
				if((data[i/(sizeof(WORD)*8)]>>(i%(sizeof(WORD)*8))&1)==1)ret++;
			}
			return ret;
		}
		
		inline void erase(){memset(data,0,kappa*sizeof(WORD));}
				
	private:
		WORD* data;
		size_t len;
		size_t kappa;
};
