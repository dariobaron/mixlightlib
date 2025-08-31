#pragma once

#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <stack>
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
	const std::map<Node::NAME,Node::ID> name_to_id_;
	const std::vector<Node> nodes_;
	const std::set<Node::ID> leaves_;
	std::vector<unsigned> depths_;
	std::vector<unsigned> widths_;
	std::vector<double> probabilities_;
	std::vector<unsigned> subtree_size_;
	std::vector<unsigned> nL_subtree_;
	double B2_, B2norm_;
	unsigned long coph_; double cophnorm_;
public:
	Tree(const std::vector<Edge> & edges);
	Tree(const Edge * ptr, unsigned n);
	Tree(const Tree & other);
	Tree(Tree && other);
	const std::vector<Node>& getNodes() const;
	std::vector<Node::NAME> getNodeNames() const;
	std::vector<Edge> getEdgelist() const;
	template<typename FlexEdge>
	static Node::NAME checkTreeIntegrity(const std::vector<FlexEdge> & edges);
	template<typename FlexEdge>
	static std::map<Node::NAME,Node::ID> getMapperNameToId(std::vector<FlexEdge> edges);
	template<typename RndEng>
	static std::vector<Edge> generateYuleEdges(RndEng & rng, unsigned nL);
	template<typename RndEng>
	static std::vector<Edge> randomizeEdges(RndEng & rng, const Tree & source);
	Tree prune(const std::vector<Node::ID> & nodes_to_keep);
	std::vector<unsigned>& computeDepths();
	std::vector<unsigned>& computeWidths();
	std::vector<double>& computeProbabilities();
	std::vector<unsigned>& computeSubtreeSizes();
	std::vector<unsigned>& computeNLeavesSubtree();
	double computeB2();
	double computeB2Norm();
	unsigned long computeCophenetic();
	double computeCopheneticNorm();
	std::vector<unsigned> computeNChildrenPerNode() const;
	unsigned computeDistanceFromAncestor(Node::NAME ancestor, Node::NAME subject) const;
private:
	std::vector<Node> initNodes(const Edge * ptr, const unsigned n);
	std::set<Node::ID> initLeaves(const std::vector<Node> & nodes);
};


Tree::Tree(const std::vector<Edge> & edges) : Tree(edges.data(), edges.size()) {};


Tree::Tree(const Edge * ptr, unsigned n) : name_to_id_(getMapperNameToId(std::vector<Edge>(ptr, ptr+n))), nodes_(initNodes(ptr,n)), leaves_(initLeaves(nodes_)), B2_(-1), B2norm_(-1), coph_(-1), cophnorm_(-1) {}

Tree::Tree(const Tree & other) : Tree(other.getEdgelist()) {};

Tree::Tree(Tree && other) : Tree(other.getEdgelist()) {};

const std::vector<Node>& Tree::getNodes() const{
	return nodes_;
}


std::vector<Node::NAME> Tree::getNodeNames() const{
	std::vector<Node::NAME> names(nodes_.size());
	for (unsigned i = 0; i < nodes_.size(); ++i){
		names[i] = nodes_[i].name();
	}
	return names;
}


std::vector<Edge> Tree::getEdgelist() const{
	std::vector<Edge> edgelist;
	for (auto & node : nodes_){
		for (auto child_ptr : node.children()){
			edgelist.emplace_back(node.name(), child_ptr->name());
		}
	}
	std::sort(edgelist.begin(), edgelist.end());
	return edgelist;
}


template<typename FlexEdge>
Node::NAME Tree::checkTreeIntegrity(const std::vector<FlexEdge> & edges){
	std::unordered_set<Node::NAME> root_candidates;
	std::unordered_set<Node::NAME> non_roots;
	std::unordered_set<Node::NAME> convergent_nodes;
	for (auto & edge : edges){
		auto [it, inserted] = non_roots.insert(edge.child);
		if (!inserted){
			convergent_nodes.insert(edge.child);
		}
		root_candidates.erase(edge.child);
		if (!non_roots.contains(edge.parent)){
			root_candidates.insert(edge.parent);
		}
	}
	if (root_candidates.size() > 1){
		std::string error_message = "Tree edgelist must have a single root node, but candidates are:";
		for (const auto & candidate : root_candidates){
			error_message += " " + std::to_string(candidate);
		}
		throw std::invalid_argument(error_message);
	}
	if (root_candidates.size() == 0){
		throw std::invalid_argument("Tree edgelist cannot be have loops");
	}
	if (convergent_nodes.size() > 0){
		std::string error_message = "Tree edgelist cannot have convergent branches, but found at nodes:";
		for (const auto & node : convergent_nodes){
			error_message += " " + std::to_string(node);
		}
		throw std::invalid_argument(error_message);
	}
	return *(root_candidates.begin());
}


template<typename FlexEdge>
std::map<Node::NAME,Node::ID> Tree::getMapperNameToId(std::vector<FlexEdge> edges){
	Node::NAME root = checkTreeIntegrity(edges);
	std::map<Node::NAME,Node::ID> name_to_id;
	name_to_id[root] = 0;
	while (edges.size()){
		std::vector<unsigned> to_erase;
		to_erase.reserve(edges.size());
		for (unsigned i = 0; i < edges.size(); ++i){
			auto it_first = name_to_id.find(edges[i].parent);
			if (it_first != name_to_id.end()){
				auto it_second = name_to_id.find(edges[i].child);
				if (it_second == name_to_id.end()){
					name_to_id[edges[i].child] = name_to_id.size();
				}
				to_erase.push_back(i);
			}
		}
		eraseWithoutOrder(edges, to_erase);
	}
	return name_to_id;
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


Tree Tree::prune(const std::vector<Node::NAME> & nodes_to_keep){
	std::set<Node::NAME> nodes_to_keep_unique(nodes_to_keep.begin(), nodes_to_keep.end());
	//<< check that all the nodes to keep belong to this tree
	for (auto name : nodes_to_keep_unique){
		name_to_id_.at(name);
	}
	//>>
	std::vector<Edge> new_edges;
	for (Node::NAME name : nodes_to_keep_unique){
		const Node & node = nodes_[name_to_id_.at(name)];
		const Node * parent_ptr = node.parent();
		if (!parent_ptr){ // if this node is the root of the tree
			continue;
		}
		auto iter = nodes_to_keep_unique.find(parent_ptr->name());
		while (iter == nodes_to_keep_unique.end()){
			parent_ptr = parent_ptr->parent();
			if (!parent_ptr){ // if this node is the root of the pruned tree
				break;
			}
			iter = nodes_to_keep_unique.find(parent_ptr->name());
		}
		if (!parent_ptr){
			continue;
		}
		new_edges.emplace_back(parent_ptr->name(), name);
	}
	std::sort(new_edges.begin(), new_edges.end());
	return Tree(new_edges);
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


std::vector<unsigned>& Tree::computeSubtreeSizes(){
	if (subtree_size_.size() != nodes_.size()){
		std::map<Node::ID,unsigned> subtrees;
		for (Node::ID leaf : leaves_){
			subtrees[leaf] = 0;
		}
		std::stack<Node::ID> to_process;
		to_process.push(0);
		while (to_process.size() > 0){
			Node::ID current_id = to_process.top();
			const Node& current = nodes_[current_id];
			bool processable = true;
			for (const Node* child : current.children()){
				if (!subtrees.contains(child->id())){
					processable = false;
					to_process.push(child->id());
				}
			}
			if (processable){
				to_process.pop();
				unsigned accumulator = current.nChildren();
				for (const Node* child : current.children()){
					accumulator += subtrees.at(child->id());
				}
				subtrees[current_id] = accumulator;
			}
		}
		subtree_size_.clear();
		subtree_size_.reserve(nodes_.size());
		for (auto & [id, size] : subtrees){
			subtree_size_.push_back(size);
		}
	}
	return subtree_size_;
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


unsigned Tree::computeDistanceFromAncestor(Node::NAME ancestor, Node::NAME subject) const{
	unsigned distance = 1;
	const Node * parent_ptr = nodes_[name_to_id_.at(subject)].parent();
	while (parent_ptr->name() != ancestor){
		++distance;
		parent_ptr = parent_ptr->parent();
		if (!parent_ptr){
			throw std::invalid_argument("The specified ancestor is not an ancestor of the subject");
		}
	}
	return distance;
}


std::vector<Node> Tree::initNodes(const Edge * ptr, const unsigned n){
	std::vector<Node> nodes(n+1);
	for (unsigned i = 0; i < n; ++i){
		const Edge edge = ptr[i];
		Node::ID parent_id = name_to_id_.at(edge.parent);
		Node::ID child_id = name_to_id_.at(edge.child);
		Node & parent = nodes[parent_id];
		Node & child = nodes[child_id];
		parent.id(parent_id);
		parent.name(edge.parent);
		parent.child(&child);
		child.id(child_id);
		child.name(edge.child);
		child.parent(&parent);
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