#pragma once

#include <cmath>
#include "stopper.hpp"


template<typename Func>
double brent(Func f, double a, double b, Stopper * stopper){
	double fa = f(a);
	double fb = f(b);
	if (std::abs(fa) < std::abs(fb))	{	std::swap(a, b);	std::swap(fa, fb);	}
	double c = a;
	double fc = f(c);
	bool flag = true;
	double s, d, fs, delta;
	while(!(*stopper)(a, b)){
		delta = std::abs(a - b) / 2;
		if ((fa == fc) || (fb == fc))	{	s = b - fb * (b - a) / (fb - fa);	}
		else							{	s = (a*fb*fc) / ((fa-fb)*(fa-fc)) + (b*fa*fc) / ((fb-fa)*(fb-fc)) + (c*fa*fb) / ((fc-fa)*(fc-fb));	}
		bool cond1 = (s < (3*a+b)/4) || (s > b);
		bool cond2 = (flag) && (std::abs(s-b) >= std::abs(b-c)/2);
		bool cond3 = (!flag) &&  (std::abs(s-b) >= std::abs(c-d)/2);
		bool cond4 = (flag) && (std::abs(b-c) < std::abs(delta));
		bool cond5 = (!flag) && (std::abs(c-d) < std::abs(delta));
		if (cond1 || cond2 || cond3 || cond4 || cond5)	{	s = (a + b) / 2;	flag = true;	}
		else											{	flag = false;	}
		fs = f(s);
		d = c;
		c = b;
		if (fa * fs < 0)	{	b = s;	fb = f(b);	}
		else				{	a = s;	fa = f(a);	}
		if (std::abs(fa) < std::abs(fb))	{	std::swap(a, b);	std::swap(fa, fb);	}
		if ((fb == 0) || (fs == 0))	{	break;	}
	}
	if (std::abs(fb) < std::abs(fs))	{	return b;	}
	else								{	return s;	}
}