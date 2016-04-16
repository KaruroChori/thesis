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

#include "swapper.h"
//Both quicksort and insertionsort are written in an inplace fashion.

template <typename OBJ>
void insertion_sort(OBJ *a, size_t n) {
	for(size_t i = 1; i < n; ++i) {
		OBJ tmp = a[i];
		size_t j = i;
		while(j > 0 && tmp < a[j - 1]) {
			a[j] = a[j - 1];
			--j;
		}
		a[j] = tmp;
	}
}

template <typename OBJ>
void insertion_sort(OBJ *a, size_t n, swapper<size_t>& swp, size_t step=0) {
	for(size_t i = 1; i < n; ++i) {
		OBJ tmp = a[swp[i+step]];
		size_t j = i;
		while(j > 0 && tmp < a[swp[j - 1+step]]) {
			swp.swap(j+step,j-1+step);
			--j;
		}
	}
}

template <typename OBJ>
void quick_sort (OBJ *a, size_t n) {
	//Threshold evaluated with direct analysis of clock timing
	
    int i, j;
    OBJ p, t;
    if (n < 2)
        return;
	//if(n<175){insertion_sort(a,n);return;}
    p = a[n / 2];
    for (i = 0, j = n - 1;; i++, j--) {
        while (a[i] < p)
            i++;
        while (p < a[j])
            j--;
        if (i >= j)
            break;
        t = a[i];
        a[i] = a[j];
        a[j] = t;
    }
    quick_sort(a, i);
    quick_sort(a + i, n - i);
}

template <typename OBJ>
void quick_sort (OBJ *a, size_t n, swapper<size_t>& swp, size_t step=0) {
	//Threshold evaluated with direct analysis of clock timing		
    int i, j;
    OBJ p;
    if(n<138){insertion_sort(a,n,swp,step);return;}
    p = a[swp[n/2+step]];
    for (i = 0, j = n - 1;; i++, j--) {
        while (a[swp[i+step]] < p)
            i++;
        while (p < a[swp[j+step]])
            j--;
        if (i >= j)
            break;
        swp.swap(i+step,j+step);
    }

    quick_sort(a, i, swp, step);
    quick_sort(a, n - i, swp, step+i);
}

//Testing code!
/*
 * 	{
	const int SIZE=200;
	
	int datam[30000],datam2[30000];
	uint64_t results[2][SIZE];
	for(int i=0;i!=2;i++)for(int j=0;j!=SIZE;j++)results[i][j]=0;
	
	for(int rpos=0;rpos!=SIZE;rpos++){
		for(int t=0;t!=1000;t++){
			for(int i=0;i!=rpos;i++)datam[i]=datam2[i]=rand()%200;
			uint64_t tt=rdtsc();
			insertion_sort(datam,rpos);
			results[0][rpos]+=rdtsc()-tt;
			
			tt=rdtsc();
			quick_sort(datam2,rpos);
			results[1][rpos]+=rdtsc()-tt;
		}
		
		cout<<rpos<<":\t"<<(double)results[0][rpos]/1000<<"\t"<<(double)results[1][rpos]/1000<<"\n";
	}
	
	}
	{
	const int SIZE=200;
	
	int datam[30000],datam2[30000];
	swapper<size_t> sw1(30000),sw2(30000);
	uint64_t results[2][SIZE];
	for(int i=0;i!=2;i++)for(int j=0;j!=SIZE;j++)results[i][j]=0;
	
	for(int rpos=0;rpos!=SIZE;rpos++){
		for(int t=0;t!=1000;t++){
			for(int i=0;i!=rpos;i++)datam[i]=datam2[i]=rand()%200;
			uint64_t tt=rdtsc();
			insertion_sort(datam,rpos,sw1);
			results[0][rpos]+=rdtsc()-tt;
			
			tt=rdtsc();
			quick_sort(datam2,rpos,sw2);
			results[1][rpos]+=rdtsc()-tt;
			
			sw2.reset();sw1.reset();
		}
		cout<<rpos<<":\t"<<(double)results[0][rpos]/1000<<"\t"<<(double)results[1][rpos]/1000<<"\n";
	}
	}
*/
