#ifndef __HUFMAN_H__
#define __HUFMAN_H__

#include<iostream>
#include<vector>
#include<algorithm>
#include<sstream>
#include<fstream>
#include<iomanip>
#include<unordered_map>
#include<utility>

using namespace std;
// make sperate class hufman inside pl.h
class hufman
{
	public:
		string symbol;
		int val;
		bool is_leaf;
		hufman *left;
		hufman *right;
		hufman(bool is_l,string name)
		{
			symbol.assign(name);
			is_leaf=is_l;
			left=NULL;
			right=NULL;
			val=0;
		}		
};
#endif // HUFMAN