#pragma once

#include <limits>
#include <stdexcept>


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
	Bounds bounds_;
public:
	Mapper(const Func & f, double a, double b);
	double a() const;
	double b() const;
	double operator()(double t) const;
	double jacobian(double t) const;
	double transformVariable(double x) const;
	double inverseTransformVariable(double t) const;
};


template <typename Func, typename T>
Mapper<Func,T>::Mapper(const Func & f, double a, double b) : f_(f), a_(a), b_(b){
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
	return jacobian(t) * f_(inverseTransformVariable(t));
}


template <typename Func, typename T>
double Mapper<Func,T>::jacobian(double t) const{
	double j;
	switch (bounds_){
	case Bounds::bounded:
		j = 0.5 * (b_ - a_);
		break;
	case Bounds::leftbounded:
		throw std::runtime_error("Mapper::jacobian: leftbounded not implemented");
		break;
	case Bounds::rightbounded:
		throw std::runtime_error("Mapper::jacobian: rightbounded not implemented");
		break;
	case Bounds::unbounded:
		throw std::runtime_error("Mapper::jacobian: unbounded not implemented");
		break;
	}
	return j;
}


template <typename Func, typename T>
double Mapper<Func,T>::transformVariable(double x) const{
	double t;
	switch (bounds_){
	case Bounds::bounded:
		t = (2*x - a_ - b_) / (b_ - a_);
		break;
	case Bounds::leftbounded:
		throw std::runtime_error("Mapper::jacobian: leftbounded not implemented");
		break;
	case Bounds::rightbounded:
		throw std::runtime_error("Mapper::jacobian: rightbounded not implemented");
		break;
	case Bounds::unbounded:
		throw std::runtime_error("Mapper::jacobian: unbounded not implemented");
		break;
	}
	return t;
}


template <typename Func, typename T>
double Mapper<Func,T>::inverseTransformVariable(double t) const{
	double x;
	switch (bounds_){
	case Bounds::bounded:
		x = t * 0.5 * (b_ - a_) + 0.5 * (a_ + b_);
		break;
	case Bounds::leftbounded:
		throw std::runtime_error("Mapper::jacobian: leftbounded not implemented");
		break;
	case Bounds::rightbounded:
		throw std::runtime_error("Mapper::jacobian: rightbounded not implemented");
		break;
	case Bounds::unbounded:
		throw std::runtime_error("Mapper::jacobian: unbounded not implemented");
		break;
	}
	return x;
}