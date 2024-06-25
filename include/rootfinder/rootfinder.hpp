#pragma once

#include <concepts>
#include <stdexcept>
#include <cmath>
#include <string>
#include "brent.hpp"
#include "stopper.hpp"


template<typename T>
concept Callable = std::invocable<T,double>;

enum class SolverMethods{
	Brent
};


template<Callable Function>
class RootFinder{
	Function function_;
	double a_;
	double b_;
	SolverMethods solver_type_;
	Stopper * tol_;
	double root_;
	bool computed_;
	bool default_stopper_;
public:
	RootFinder(Function function, double a, double b);
	RootFinder(Function function, double a, double b, double rel_tol, double abs_tol, unsigned max_iter);
	RootFinder(Function function, double a, double b, Stopper * tol);
	~RootFinder();
	void setSolver(std::string method);
	double root();
	unsigned iterations() const;
};



template<Callable Function>
RootFinder<Function>::RootFinder(Function function, double a, double b)
				: function_(function), a_(a), b_(b), solver_type_(SolverMethods::Brent), computed_(false), default_stopper_(true) {
	if (a_ == b_){
		throw std::runtime_error("The interval endpoints cannot be equal!");
	}
	if (function_(a_) * function_(b_) > 0){
		throw std::runtime_error("The function must evaluate with different signs in the two endpoints!");
	}
	if (a_ > b_){
		std::swap(a_, b_);
	}
	tol_ = new StopperDefault();
}


template<Callable Function>
RootFinder<Function>::RootFinder(Function function, double a, double b, double rel_tol, double abs_tol, unsigned max_iter) : RootFinder(function, a, b) {
	delete tol_;
	tol_ = StopperDefault(rel_tol, abs_tol, max_iter);
}


template<Callable Function>
RootFinder<Function>::RootFinder(Function function, double a, double b, Stopper * tol) : RootFinder(function, a, b) {
	tol_ = tol;
	default_stopper_ = false;
}


template<Callable Function>
RootFinder<Function>::~RootFinder(){
	if (default_stopper_)	{	delete tol_;	}
}


template<Callable Function>
double RootFinder<Function>::root(){
	if (!computed_){
		computed_ = true;
		switch (solver_type_){
		case SolverMethods::Brent:
			root_ = brent(function_, a_, b_, tol_);
			break;
		}
	}
	return root_;
}


template<Callable Function>
unsigned RootFinder<Function>::iterations() const{
	return tol_->iterations();
}
