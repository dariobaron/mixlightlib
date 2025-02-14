#pragma once

#include <list>
#include <numeric>
#include <stdexcept>
#include <algorithm>
#include "gausskronrod.hpp"
#include "mapper.hpp"


template <typename Func, typename T=double>
class Integrator{
public:
	using Solver = GaussKronrod<Func,T>;
	using FuncWrapper = Mapper<Func,T>;
private:
	const Func & f_;
	double a_, b_;
	double tolerance_;
	double result_;
	double error_;
	unsigned maxiter_;
	unsigned nSubintervals_;
	bool success_;
public:
	Integrator(const Func & f, double a, double b, double tolerance, unsigned maxiter=1000);
	double a() const;
	double b() const;
	double tolerance() const;
	double result() const;
	double error() const;
	unsigned nSubintervals() const;
	bool success() const;
private:
	void integrate();
};


template <typename Func, typename T>
Integrator<Func,T>::Integrator(const Func & f, double a, double b, double tolerance, unsigned maxiter)
												: f_(f), a_(std::min(a,b)), b_(std::max(a,b)), tolerance_(tolerance),
												result_(0), error_(0), maxiter_(maxiter), nSubintervals_(0), success_(false)
{
	integrate();
}


template <typename Func, typename T>
double Integrator<Func,T>::b() const{
	return b_;
}


template <typename Func, typename T>
double Integrator<Func,T>::tolerance() const{
	return tolerance_;
}


template <typename Func, typename T>
double Integrator<Func,T>::result() const{
	return result_;
}


template <typename Func, typename T>
double Integrator<Func,T>::error() const{
	return error_;
}


template <typename Func, typename T>
unsigned Integrator<Func,T>::nSubintervals() const{
	return nSubintervals_;
}


template <typename Func, typename T>
bool Integrator<Func,T>::success() const{
	return success_;
}


template <typename Func, typename T>
void Integrator<Func,T>::integrate(){
	std::list<Solver> subintegrals;
	{
		FuncWrapper mapper(f_, a_, b_);
		subintegrals.emplace_back(mapper, tolerance_);
	}
	auto current = subintegrals.begin();
	while (current != subintegrals.end()){
		current->integrate();
		if (current->success()){
			++current;
		}
		else{
			auto to_remove = current;
			double a = current->a();
			double b = current->b();
			FuncWrapper mapper(f_, a, b);
			double a_t = mapper.transformVariable(a);
			double b_t = mapper.transformVariable(b);
			double m_t = 0.5 * (a_t + b_t);
			double m = mapper.inverseTransformVariable(m_t);
			double new_tolerance = 0.5 * current->tolerance();
			current = subintegrals.insert(current, {Solver(FuncWrapper(f_, a, m), new_tolerance), Solver(FuncWrapper(f_, m, b), new_tolerance)});
			subintegrals.erase(to_remove);
		}
		if (subintegrals.size() > maxiter_){
			while (current != subintegrals.end()){
				current->integrate();
				++current;
			}
			break;
		}
	}
	subintegrals.sort([](const Solver & left, const Solver & right){return left.result() < right.result();});
	result_ = std::accumulate(subintegrals.begin(), subintegrals.end(), 0.0, [](double sum, const Solver & integral){return sum + integral.result();});
	subintegrals.sort([](const Solver & left, const Solver & right){return left.error() < right.error();});
	error_ = std::accumulate(subintegrals.begin(), subintegrals.end(), 0.0, [](double sum, const Solver & integral){return sum + integral.error();});
	nSubintervals_ = subintegrals.size();
	success_ = std::all_of(subintegrals.begin(), subintegrals.end(), [](const Solver & integral){return integral.success();});
}


template <typename Func, typename T>
double Integrator<Func,T>::a() const{
	return a_;
}