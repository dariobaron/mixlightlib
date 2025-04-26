#pragma once

#include <vector>
#include <set>
#include <map>
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
	std::vector<unsigned> widths_;
	std::vector<double> probabilities_;
	std::vector<unsigned> nL_subtree_;
	double B2_, B2norm_;
	unsigned long coph_; double cophnorm_;
public:
	Tree(const std::vector<Edge> & edges);
	Tree(const Edge * ptr, unsigned n);
	Tree(const Tree & other);
	Tree(Tree && other);
	const std::vector<Node>& getNodes() const;
	std::vector<Edge> getEdgelist() const;
	template<typename RndEng>
	static std::vector<Edge> generateYuleEdges(RndEng & rng, unsigned nL);
	template<typename RndEng>
	static std::vector<Edge> randomizeEdges(RndEng & rng, const Tree & source);
	template<typename FlexEdge,typename NameType>
	static std::vector<Edge> renameEdgelist(std::vector<FlexEdge> edges, NameType root);
	std::vector<unsigned>& computeDepths();
	std::vector<unsigned>& computeWidths();
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

Tree::Tree(const Tree & other) : Tree(other.getEdgelist()) {};

Tree::Tree(Tree && other) : Tree(other.getEdgelist()) {};

const std::vector<Node>& Tree::getNodes() const{
	return nodes_;
}


std::vector<Edge> Tree::getEdgelist() const{
	std::vector<Edge> edgelist;
	for (auto & node : nodes_){
		for (auto child_ptr : node.children()){
			edgelist.emplace_back(node.id(), child_ptr->id());
		}
	}
	std::sort(edgelist.begin(), edgelist.end());
	return edgelist;
}


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


template<typename FlexEdge,typename NameType>
std::vector<Edge> Tree::renameEdgelist(std::vector<FlexEdge> edges, NameType root){
	auto copy = edges;
	bool root_notin_edges = true;
	for (const auto & edge : edges){
		if (edge.parent == root){
			root_notin_edges = false;
		}
		if (edge.child == root){
			throw std::invalid_argument("Root requested cannot be a child");
		}
	}
	if (root_notin_edges){
		throw std::invalid_argument("Root requested is not in the edgelist");
	}
	auto initial_size = edges.size();
	std::map<NameType,Node::ID> node_names;
	node_names[root] = 0;
	std::vector<Edge> renamed_edges;
	renamed_edges.reserve(edges.size());
	while (edges.size()){
		std::vector<unsigned> to_erase;
		to_erase.reserve(edges.size());
		for (unsigned i = 0; i < edges.size(); ++i){
			auto it_first = node_names.find(edges[i].parent);
			if (it_first != node_names.end()){
				auto it_second = node_names.find(edges[i].child);
				if (it_second == node_names.end()){
					node_names[edges[i].child] = node_names.size();
				}
				renamed_edges.emplace_back(node_names[edges[i].parent], node_names[edges[i].child]);
				to_erase.push_back(i);
			}
		}
		eraseWithoutOrder(edges, to_erase);
	}
	if (renamed_edges.size() != initial_size){
		std::cout << "renamed_edges.size()=" << renamed_edges.size() << " edges.size()=" << initial_size << std::endl;
		throw std::runtime_error("Error in the algorithm");
	}
	if (node_names.size() != initial_size + 1){
		std::set<Edge> edges_set;
		for (const auto & edge : copy){
			auto [it, inserted] = edges_set.insert(edge);
			if (!inserted){
				std::cout << "duplicate edge found=" << edge.parent << " -> " << edge.child << std::endl;
				auto it = edges_set.find(edge);
				std::cout << "duplicate edge found=" << it->parent << " -> " << it->child << std::endl;
			}
			if (node_names.find(edge.parent) == node_names.end()){
				std::cout << "parent not found=" << edge.parent << std::endl;
			}
			if (node_names.find(edge.child) == node_names.end()){
				std::cout << "child not found=" << edge.child << std::endl;
			}
		}
		std::cout << " edges_set.size()=" << edges_set.size() << std::endl;
		std::cout << "node_names.size()=" << node_names.size() << " initial_size=" << initial_size << std::endl;
		throw std::runtime_error("Error in the algorithm");
	}
	std::sort(renamed_edges.begin(), renamed_edges.end());
	return renamed_edges;
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


std::vector<unsigned>& Tree::computeWidths(){
	if (widths_.size() == 0){
		widths_.push_back(1);
		std::vector<const Node*> this_level = nodes_[0].children();
		std::vector<const Node*> next_level;
		while (this_level.size()){
			widths_.push_back(this_level.size());
			for (const Node * node_ptr : this_level){
				auto tmp = node_ptr->children();
				next_level.insert(next_level.end(), tmp.begin(), tmp.end());
			}
			this_level = next_level;
			next_level.clear();
		}
	}
	return widths_;
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
		unsigned already_computed = 0;
		for (auto l : leaves_){
			const Node * parent_ptr = nodes_[l].parent();
			while (parent_ptr){
				Node::ID parent_idx = parent_ptr->id();
				++nL_subtree_[parent_idx];
				parent_ptr = nodes_[parent_idx].parent();
			}
			++already_computed;
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
		if (nL <= 4){
			B2norm_ = -1;	
		}
		else if (B2 <= 2 * (1 - std::pow(2,1-nL)) + 1e-5){
			B2norm_ = 0;
		}
		else{
			double log2nL = std::log2(nL);
			auto function = [B2,log2nL](double b){ return -B2 + 2*(b-1)/(b+0.0918) + log2nL/b; };
			RootFinder finder(function, log2nL/B2, log2nL/(B2-2*(1-std::pow(2,1.-nL))));
			double b = finder.root();
			double log2e = 1 / std::log(2);
			double result = std::pow(2, (-b+1)/(log2e-1));
			double floored = std::floor(log2nL);
			double discretization_correction = log2nL / (floored + (nL - std::pow(2,floored)) / std::pow(2,floored));
			B2norm_ = result * std::pow(discretization_correction,2.5);
		}
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
		throw std::runtime_error("Edgelist must be begin with the root node Id = 0");
	}
	nodes[0].id(0);
	for (unsigned i = 0; i < n; ++i){
		Edge edge = ptr[i];
		if (edge.parent >= edge.child)	{	throw std::runtime_error("Tree edgelist ill-formed!");	}
		else if (edge.child > n)		{	throw std::runtime_error("Tree edgelist ill-formed!");	}
		nodes[edge.child].id(edge.child);
		nodes[edge.parent].child(&nodes[edge.child]);
		nodes[edge.child].parent(&nodes[edge.parent]);
	}
	for (unsigned i = 0; i < nodes.size(); ++i){
		if (nodes[i].id() != i)	{	throw std::runtime_error("Edgelist is incomplete, some nodes are missing");	}
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