#pragma once

#include <vector>


class Node{
public:
	using ID = unsigned;
private:
	ID id_;
	Node * parent_;
	std::vector<Node*> children_;
public:
	Node();
	Node(ID id);
	ID id() const;
	void id(ID i);
	const Node * parent() const;
	void parent(Node * p);
	unsigned nChildren() const;
	std::vector<const Node*> children() const;
	void children(std::vector<Node*> & c);
	const Node * child(unsigned i) const;
	void child(Node * c);
};



Node::Node() : id_(-1), parent_(nullptr) {}


Node::Node(ID id) : id_(id), parent_(nullptr) {}


Node::ID Node::id() const{
	return id_;
}

void Node::id(ID i){
	id_ = i;
}


const Node * Node::parent() const{
	return parent_;
}

void Node::parent(Node * p){
	parent_ = p;
}


unsigned Node::nChildren() const{
	return children_.size();
}


std::vector<const Node*> Node::children() const{
	std::vector<const Node*> ret_val(children_.size());
	for (unsigned i = 0; i < children_.size(); ++i){
		ret_val[i] = children_[i];
	}
	return ret_val;
}

void Node::children(std::vector<Node*> & c){
	children_ = c;
}


const Node * Node::child(unsigned i) const{
	return children_[i];
}

void Node::child(Node * c){
	children_.push_back(c);
}