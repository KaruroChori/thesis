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
#include <fstream>

using namespace std;

/* I have to warn thee joung wayfarer. In this file thou findest how not
 * to implement stream wrappers in the proper way. The swoopstake usage
 * of these classes may be full o'dangers.
 * Be prepared. Stay away, stay alive.
*/

template<typename OSTREAM>
class bit_ofstream{
    public:
        bit_ofstream(OSTREAM& s):flux(s){buffer=0;pos=0;};
        ~bit_ofstream(){close();}
       
        friend bit_ofstream& operator<<(bit_ofstream& out, bool data){
            if(out.pos==8){
                out.flush();
            }
            out.buffer|=data<<out.pos;
            out.pos++;
            return out;
        }
       
        int write(uint64_t data, int bits){
            if(bits>64)    throw "Max 64 bits";
           
            for(;bits!=0;bits--){operator<<(*this,data&1);data>>=1;}
            return 0;
        }
       
        int flush(){flux<<buffer;buffer=0;pos=0;return 0;}
        int close(){flush();return 0;}
       
        OSTREAM& stream(){return flux;}
       
    private:
        OSTREAM& flux;
        uint8_t buffer;
        int pos;
};

template<typename ISTREAM>
class bit_ifstream{
    public:
        bit_ifstream(ISTREAM& s):flux(s){buffer=0;pos=0;};
        ~bit_ifstream(){close();}
       
        friend bit_ifstream& operator>>(bit_ifstream& in, bool& data){
            if(in.pos==8){
                in.sync();
            }
            data=(in.buffer&(1<<in.pos))>>in.pos;
            in.pos++;
            return in;
        }
       
        uint64_t read(int bits){
            if(bits>64)    throw "Max 64 bits";       
            uint64_t temp=0;
            bool value;
            for(int i=0;i!=bits;bits++){operator>>(*this,value);temp|=value<<i;}
            return temp;
        }
       
        int sync(){flux>>buffer;pos=0;return 0;}
        int close(){buffer=0;pos=0;return 0;}
       
        ISTREAM& stream(){return flux;}
       
    private:
        ISTREAM& flux;
        uint8_t buffer;
        int pos;
};
