#pragma once

#include <limits>
#include <stdexcept>
#include <cmath>

////////////////////////
// This class allows to use the GaussKronrod class
// that integrates the normalized variable t in (-1,1)
// in a general integral of the variable x over (a,b)
// as in $\int_{ginv(a)}^{ginv(b)} f(g(t)) g'(t) dt = \int_a^b f(x) dx$.
// The map x = g(t) is called the transformation and t = g^-1(x) is its inverse.
//
// Depeneding on the finiteness of a, b, the following transformations are used:
// - a finite, b finite:
//			x = t * (b-a)/2 + (a+b)/2
//			t = (2x-b-a) / (b-a)
//			dx/dt = (b-a)/2
// - a infinite, b finite:
//			x = b - scale * (1-t) / (1+t)
//			t = - (x-b+scale) / (x-b-scale)
//			dx/dt = 2*scale / (1+t)^2
// - a finite, b infinite:
//			x = a + scale * (1+t) / (1-t)
//			t = (x-a-scale) / (x-a+scale)
//			dx/dt = 2*scale / (1-t)^2
// - a infinite, b infinite:
//			x = loc + 3/4 * scale * t / (1-t^2)
//			t = (-3/4*scale + \sqrt{9/16*scale^2 + 4*(x-loc)^2}) / (2*(x-loc))
//			dx/dt = 3/4 * scale * (1+t^2) / (1-t^2)^2
//
// The parameter `scale` is only considered if the interval is not bounded,
// the parameter `loc` is only considered if the interval is unbounded.
// In case of half-bound interval, the parameter `scale` represents the distance
// from the finite bound at which the interval is splitted at the next iteration. It
// must be chosen such that the splitting happens within a relevant region of
// evaluation of the function.
// In case of unbounded interval, the parameter `loc` represents the location
// at which the interval is splitted at the next iteration. It must be chosen
// at the centre of the relevant region of evaluation of the function.
// In case of unbounded interval, the parameter `scale` represents the amplitude
// of the relevant region of evaluation of the function.
////////////////////////

template <typename Func, typename T=double>
class Mapper{
public:
	enum class Bounds{
		bounded,
		leftbounded,
		rightbounded,
		unbounded
	};
private:
	const Func & f_;
	double a_, b_;
	double scale_, loc_;
	Bounds bounds_;
public:
	Mapper(const Func & f, double a, double b, double scale=1, double loc=0);
	double a() const;
	double b() const;
	double operator()(double t) const;
	double inverseTransformVariable(double x) const;
	double transformVariable(double t) const;
	double jacobian(double t) const;
};


template <typename Func, typename T>
Mapper<Func,T>::Mapper(const Func & f, double a, double b, double scale, double loc) : f_(f), a_(a), b_(b), scale_(scale), loc_(loc) {
	bool is_a_inf = a_ == -std::numeric_limits<double>::infinity();
	bool is_b_inf = b_ == std::numeric_limits<double>::infinity();
	bounds_ = static_cast<Bounds>(2 * is_a_inf + is_b_inf);
}


template <typename Func, typename T>
double Mapper<Func,T>::a() const{
	return a_;
}


template <typename Func, typename T>
double Mapper<Func,T>::b() const{
	return b_;
}


template <typename Func, typename T>
double Mapper<Func,T>::operator()(double t) const{
	return jacobian(t) * f_(transformVariable(t));
}


template <typename Func, typename T>
double Mapper<Func,T>::inverseTransformVariable(double x) const{
	if (x == a_){
		return -1;
	} else if (x == b_){
		return 1;
	}
	double t;
	switch (bounds_){
	case Bounds::bounded:
		t = (2*x - b_ - a_) / (b_ - a_);
		break;
	case Bounds::leftbounded:
		t = (x - a_ - scale_) / (x - a_ + scale_);
		break;
	case Bounds::rightbounded:
		t = - (x - b_ + scale_) / (x - b_ - scale_);
		break;
	case Bounds::unbounded:
		t = (-3./4.*scale_ + std::sqrt(std::pow(3./4.*scale_, 2) + 4*std::pow(x-loc_, 2))) / (2 * (x-loc_));
		break;
	default:
		throw std::runtime_error("Mapper::inverseTransformVariable: unknown bounds");
		break;
	}
	return t;
}


template <typename Func, typename T>
double Mapper<Func,T>::transformVariable(double t) const{
	double x;
	switch (bounds_){
	case Bounds::bounded:
		x = t * 0.5 * (b_ - a_) + 0.5 * (a_ + b_);
		break;
	case Bounds::leftbounded:
		x = a_ + scale_ * (1+t) / (1-t);
		break;
	case Bounds::rightbounded:
		x = b_ - scale_ * (1-t) / (1+t);
		break;
	case Bounds::unbounded:
		x = loc_ + 3./4. * scale_ * t / (1-t*t);
		break;
	default:
		throw std::runtime_error("Mapper::transformVariable: unknown bounds");
		break;
	}
	return x;
}

	
template <typename Func, typename T>
double Mapper<Func,T>::jacobian(double t) const{
	double j;
	switch (bounds_){
	case Bounds::bounded:
		j = 0.5 * (b_ - a_);
		break;
	case Bounds::leftbounded:
		j = 2 * scale_ / std::pow(1-t, 2);
		break;
	case Bounds::rightbounded:
		j = 2 * scale_ / std::pow(1+t, 2);
		break;
	case Bounds::unbounded:
		j = 3./4. * scale_ * (1+t*t) / std::pow(1-t*t, 2);
		break;
	default:
		throw std::runtime_error("Mapper::jacobian: unknown bounds");
		break;
	}
	return j;
}