#pragma once

#include <tuple>


struct Edge{
	unsigned parent;
	unsigned child;
	Edge(unsigned p, unsigned c) : parent(p), child(c) {};
	friend bool operator<(const Edge & lhs, const Edge & rhs);
};



bool operator<(const Edge & lhs, const Edge & rhs){
	if (lhs.child == rhs.child){
		return lhs.parent < rhs.parent;
	}
	return lhs.child < rhs.child;
}