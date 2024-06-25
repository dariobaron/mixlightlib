#pragma once

#include <tuple>


struct Edge{
	unsigned parent;
	unsigned child;
	Edge(unsigned p, unsigned c) : parent(p), child(c) {};
	friend bool operator<(const Edge & lhs, const Edge & rhs);
};



bool operator<(const Edge & lhs, const Edge & rhs){
	return std::tie(lhs.parent, rhs.child) < std::tie(rhs.parent, rhs.child);
}