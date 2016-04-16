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

template<typename SCALAR, typename OBJECT>
struct generic_distance{
	SCALAR lower_bound;
	SCALAR upper_bound;
	
	struct naked{
		SCALAR lower_bound;
		SCALAR upper_bound;
		naked(){}
		naked(const generic_distance& r){lower_bound=r.lower_bound;upper_bound=r.upper_bound;}
	};
	
	virtual generic_distance& operator()(OBJECT a, OBJECT b, SCALAR inf, SCALAR sup)=0;
	virtual bool zero_implies_same() const =0;
	virtual bool same_implies_zero() const=0;
	virtual bool symmetric() const=0;
	virtual bool triangular_inequality() const=0;
	virtual bool non_negative() const=0;
		
	virtual void print(std::ostream& out) const{
		if(this->lower_bound==this->upper_bound)out<<this->lower_bound;
		else out<<"(l:"<<this->lower_bound<<", h:"<<this->upper_bound<<")";
	}


	friend std::ostream& operator<<(std::ostream& out, const generic_distance & d){
		d.print(out);
		return out;
	}
	
};
