#pragma once

#include <cmath>


class Stopper{
public:
	virtual bool operator()(double xold, double xnew) = 0;
	virtual unsigned iterations() const = 0;
};


class StopperDefault : public Stopper{
	double rel_tol_;
	double abs_tol_;
	unsigned max_iter_;
	unsigned iter_;
public:
	StopperDefault() : rel_tol_(1e-10), abs_tol_(1e-10), max_iter_(1000), iter_(0) {};
	StopperDefault(double rel, double abs, unsigned maxit) : rel_tol_(rel), abs_tol_(abs), max_iter_(maxit), iter_(0) {};
	bool operator()(double xold, double xnew){
		++iter_;
		double delta = std::abs(xnew - xold);
		bool rel_reached = delta / std::abs(xold) < rel_tol_;
		bool abs_reached = delta < abs_tol_;
		bool max_iter_reached = max_iter_ < iter_;
		return rel_reached || abs_reached || max_iter_reached;
	};
	unsigned iterations() const{
		return iter_;
	};
};