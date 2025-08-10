#pragma once

#include <tuple>
#include "node.hpp"


struct Edge{
	Node::NAME parent;
	Node::NAME child;
	Edge(Node::NAME p, Node::NAME c) : parent(p), child(c) {};
	friend bool operator<(const Edge & lhs, const Edge & rhs);
};



bool operator<(const Edge & lhs, const Edge & rhs){
	if (lhs.child == rhs.child){
		return lhs.parent < rhs.parent;
	}
	return lhs.child < rhs.child;
}