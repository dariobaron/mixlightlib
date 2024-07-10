#pragma once

#include <vector>
#include <set>
#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <random>
#include "node.hpp"
#include "edge.hpp"
#include "../rootfinder.hpp"
#include "../algorithms.hpp"
#include "../functions.hpp"


class Tree{
	const std::vector<Node> nodes_;
	const std::set<Node::ID> leaves_;
	std::vector<unsigned> depths_;
	std::vector<double> probabilities_;
	std::vector<unsigned> nL_subtree_;
	double B2_, B2norm_;
	unsigned long coph_; double cophnorm_;
public:
	Tree(const std::vector<Edge> & edges);
	Tree(const Edge * ptr, unsigned n);
	template<typename RndEng>
	static std::vector<Edge> generateYuleEdges(RndEng & rng, unsigned nL);
	template<typename RndEng>
	static std::vector<Edge> randomizeEdges(RndEng & rng, const Tree & source);
	std::vector<unsigned>& computeDepths();
	std::vector<double>& computeProbabilities();
	std::vector<unsigned>& computeNLeavesSubtree();
	double computeB2();
	double computeB2Norm();
	unsigned long computeCophenetic();
	double computeCopheneticNorm();
	std::vector<unsigned> computeNChildrenPerNode() const;
private:
	std::vector<Node> initNodes(const Edge * ptr, unsigned n);
	std::set<Node::ID> initLeaves(const std::vector<Node> & nodes);
};


Tree::Tree(const std::vector<Edge> & edges) : Tree(edges.data(), edges.size()) {};


Tree::Tree(const Edge * ptr, unsigned n) : nodes_(initNodes(ptr,n)), leaves_(initLeaves(nodes_)), B2_(-1), B2norm_(-1), coph_(-1), cophnorm_(-1) {}


template<typename RngEng>
std::vector<Edge> Tree::generateYuleEdges(RngEng & rng, unsigned nL){
	std::vector<Edge> edges = {{0,1}, {0,2}};
	std::vector<Node::ID> leaves = {1, 2};
	while (leaves.size() < nL){
		unsigned idx = std::uniform_int_distribution<>(0, leaves.size()-1)(rng);
		Node::ID parent = leaves[idx];
		Node::ID child1 = edges.size() + 1;
		Node::ID child2 = edges.size() + 2;
		edges.emplace_back(parent,child1);
		edges.emplace_back(parent,child2);
		eraseWithoutOrder(leaves, idx);
		leaves.push_back(child1);
		leaves.push_back(child2);
	}
	return edges;
}


template<typename RngEng>
std::vector<Edge> Tree::randomizeEdges(RngEng & rng, const Tree & source){
	unsigned Nnodes = source.nodes_.size();
	std::vector<unsigned> NsChildren(Nnodes);
	for (unsigned i = 0; i < Nnodes; ++i){
		NsChildren[i] = source.nodes_[i].nChildren();
	}
	std::shuffle(NsChildren.begin(), NsChildren.end(), rng);
	// check that can eventually be removed
	if (Nnodes - 1 != std::accumulate(NsChildren.begin(), NsChildren.end(), 0)){
		std::string error = "Error in the algorithm of randomize:";
		error += "Nnodes=" + std::to_string(Nnodes) + " vs NtotChildren=";
		error += std::to_string(std::accumulate(NsChildren.begin(),NsChildren.end(),0));
		throw std::runtime_error(error);
	}
	// end check
	std::vector<Edge> edges;
	Node::ID current_id = 0;
	std::vector<Node::ID> leaves = {current_id};
	++current_id;
	for (auto nchildren : NsChildren){
		unsigned leaf_idx = std::uniform_int_distribution<>(0, leaves.size()-1)(rng);
		Node::ID leaf_reproducing = leaves[leaf_idx];
		for (unsigned i = 0; i < nchildren; ++i){
			edges.emplace_back(leaf_reproducing, current_id);
			leaves.push_back(current_id);
			++current_id;
		}
		eraseWithoutOrder(leaves, leaf_idx);
	}
	std::sort(edges.begin(), edges.end());
	return edges;
}


std::vector<unsigned>& Tree::computeDepths(){
	if (depths_.size() != nodes_.size()){
		depths_ = std::vector<unsigned>(nodes_.size(), -1);
		depths_[0] = 0;
		for (const Node & node : nodes_){
			Node::ID parent = node.id();
			for (unsigned i = 0; i < node.nChildren(); ++i){
				Node::ID child = node.child(i)->id();
				depths_[child] = depths_[parent] + 1;
			}
		}
	}
	return depths_;
}


std::vector<double>& Tree::computeProbabilities(){
	if (probabilities_.size() != nodes_.size()){
		probabilities_ = std::vector<double>(nodes_.size(), -1);
		probabilities_[0] = 1;
		for (const Node & node : nodes_){
			Node::ID parent = node.id();
			for (unsigned i = 0; i < node.nChildren(); ++i){
				Node::ID child = node.child(i)->id();
				probabilities_[child] = probabilities_[parent] / node.nChildren();
			}
		}
	}
	return probabilities_;
}


std::vector<unsigned>& Tree::computeNLeavesSubtree(){
	if (nL_subtree_.size() != nodes_.size()){
		nL_subtree_ = std::vector<unsigned>(nodes_.size(), 0);
		for (auto & l : leaves_){
			nL_subtree_[l] = 1;
		}
		unsigned already_computed = leaves_.size();
		const Node * current = nodes_.data();
		while (already_computed == nodes_.size()){
			unsigned nChildren = current->nChildren();
			if (nChildren == 0){
				current = current->parent();
			}
			else{
				unsigned nL = 0;
				for (auto child : current->children()){
					unsigned nL_child = nL_subtree_[child->id()];
					if (nL_child == 0){
						current = child;
						nL = 0;
						break;
					}
					else{
						nL += nL_child;
					}
				}
				if (nL != 0){
					Node::ID idx = std::distance(std::as_const(nodes_).data(), current);
					nL_subtree_[idx] = nL;
					++already_computed;
				}
			}
		}
		for (auto & nL : nL_subtree_){
			if (nL == 0)	{	throw std::runtime_error("Algorithm failed");	}
		}
	}
	return nL_subtree_;
}


double Tree::computeB2(){
	if (B2_ == -1){
		computeProbabilities();
		const unsigned nL = leaves_.size();
		std::vector<double> plogp(nL);
		Node::ID idx = 0;
		for (auto & leaf : leaves_){
			double px = probabilities_[leaf];
			plogp[idx] = px * std::log2(px);
			++idx;
		}
		std::sort(plogp.begin(), plogp.end());
		B2_ = - std::accumulate(plogp.begin(), plogp.end(), 0.);
	}
	return B2_;
}


double Tree::computeB2Norm(){
	if (B2norm_ == -1){
		double B2 = computeB2();
		unsigned nL = leaves_.size();
		auto function = [B2,nL](double b){ return -B2 + 2*(b-1)/(b+0.0918) + std::log2(nL)/b; };
		RootFinder finder(function, std::log2(nL)/B2, std::log2(nL)/(B2-2));
		double b = finder.root();
		double log2e = 1 / std::log(2);
		B2norm_ = std::pow(2, (-b+1)/(log2e-1));
	}
	return B2norm_;
}


unsigned long Tree::computeCophenetic(){
	if (coph_ == -1){
		computeNLeavesSubtree();
		coph_ = 0;
		for (auto it = nL_subtree_.begin() + 1; it != nL_subtree_.end(); ++it){
			if ((*it) != 1){
				coph_ += binomCoeff(*it, 2);
			}
		}
	}
	return coph_;
}


double Tree::computeCopheneticNorm(){
	if (cophnorm_ == -1){
		double coph = computeCophenetic();
		cophnorm_ = coph / binomCoeff(leaves_.size(), 3);
	}
	return cophnorm_;
}


std::vector<unsigned> Tree::computeNChildrenPerNode() const{
	std::vector<unsigned> nchildren(nodes_.size());
	for (unsigned i = 0; i < nodes_.size(); ++i){
		nchildren[i] = nodes_[i].nChildren();
	}
	return nchildren;
}


std::vector<Node> Tree::initNodes(const Edge * ptr, unsigned n){
	std::vector<Node> nodes(n+1);
	if (ptr->parent != 0){
		throw std::runtime_error("Edgelist must be sorted and begin with the root node Id = 0");
	}
	nodes[0].id(0);
	for (unsigned i = 0; i < n; ++i){
		Edge edge = ptr[i];
		if ((edge.parent >= edge.child) || (edge.child > i+1)){
			throw std::runtime_error("Tree edgelist ill-formed!");
		}
		if (nodes[edge.parent].id() != edge.parent){
			throw std::runtime_error("Error in the algorithm!");
		}
		nodes[edge.child].id(edge.child);
		nodes[edge.parent].child(&nodes[edge.child]);
		nodes[edge.child].parent(&nodes[edge.parent]);
	}
	return nodes;
}


std::set<Node::ID> Tree::initLeaves(const std::vector<Node> & nodes){
	std::set<Node::ID> leaves;
	for (const Node & node : nodes){
		if (node.nChildren() == 0){
			leaves.insert(leaves.end(), node.id());
		}
	}
	return leaves;
}