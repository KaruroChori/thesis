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

template<typename WORD>
inline size_t countpop(WORD data);

template<>
inline size_t countpop<uint64_t>(uint64_t data){return __builtin_popcountll(data);}

template<>
inline size_t countpop<uint32_t>(uint32_t data){return __builtin_popcount(data);}
