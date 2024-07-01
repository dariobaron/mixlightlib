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


class Tree{
	std::vector<Node> nodes_;
	std::set<Node::ID> leaves_;
	std::vector<unsigned> depths_;
	std::vector<double> probabilities_;
	double B2_, B2norm_;
public:
	Tree(std::vector<Edge> & edges);
	Tree(const Edge * ptr, unsigned n);
	template<typename RndEng>
	static std::vector<Edge> generateYuleEdges(RndEng & rng, unsigned nL);
	std::vector<unsigned>& computeDepths();
	std::vector<double>& computeProbabilities();
	double computeB2();
	double computeB2Norm();
};


Tree::Tree(std::vector<Edge> & edges) : Tree(edges.data(), edges.size()) {};


Tree::Tree(const Edge * ptr, unsigned n) : nodes_(n+1), B2_(-1), B2norm_(-1) {
	if (ptr->parent != 0){
		throw std::runtime_error("Edgelist must be sorted and begin with the root node Id = 0");
	}
	nodes_[0].id(0);
	for (unsigned i = 0; i < n; ++i){
		Edge edge = ptr[i];
		if ((edge.parent >= edge.child) || (edge.child > i+1)){
			throw std::runtime_error("Tree edgelist ill-formed!");
		}
		if (nodes_[edge.parent].id() != edge.parent){
			throw std::runtime_error("Error in the algorithm!");
		}
		nodes_[edge.child].id(edge.child);
		nodes_[edge.parent].child(&nodes_[edge.child]);
		nodes_[edge.child].parent(&nodes_[edge.parent]);
	}
	for (Node & node : nodes_){
		if (node.nChildren() == 0){
			leaves_.insert(leaves_.end(), node.id());
		}
	}
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


std::vector<unsigned>& Tree::computeDepths(){
	if (depths_.size() != nodes_.size()){
		depths_ = std::vector<unsigned>(nodes_.size(), -1);
		depths_[0] = 0;
		for (Node & node : nodes_){
			unsigned parent = node.id();
			for (unsigned i = 0; i < node.nChildren(); ++i){
				unsigned child = node.child(i)->id();
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
		for (Node & node : nodes_){
			unsigned parent = node.id();
			for (unsigned i = 0; i < node.nChildren(); ++i){
				unsigned child = node.child(i)->id();
				probabilities_[child] = probabilities_[parent] / node.nChildren();
			}
		}
	}
	return probabilities_;
}


double Tree::computeB2(){
	if (B2_ == -1){
		computeProbabilities();
		const unsigned nL = leaves_.size();
		std::vector<double> plogp(nL);
		unsigned idx = 0;
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