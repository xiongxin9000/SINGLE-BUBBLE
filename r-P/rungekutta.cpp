#include "RungeKutta.h"
#include <iostream>

RungeKutta::RungeKutta() : pf(0), pg(0),h(0){ }

RungeKutta::RungeKutta(func_ptr f, func_ptr g,double step) : pf(f), pg(g),h(step) {}
// RungeKutta::~RungeKutta() 
// {
// 	std::cout<<"deconstructor called"<<std::endl;
// }

RungeKutta& RungeKutta::setFunc(func_ptr f,func_ptr g)
{
	pf = f;
	pg = g;
	return *this;
}

RungeKutta& RungeKutta::setStep(double step)
{
	h = step;
	return *this;
}

double RungeKutta::operator()(double x, double y,double &z) const
{
	double k1 = h*pf(x,y,z);
	double l1 = h*pg(x,y,z);
	double k2 = h*pf(x+h/2.0,y+k1/2.0,z+l1/2.0);
	double l2 = h*pg(x+h/2.0,y+k1/2.0,z+l1/2.0);
	double k3 = h*pf(x+h/2.0,y+k2/2.0,z+l2/2.0);
	double l3 = h*pg(x+h/2.0,y+k2/2.0,z+l2/2.0);
	double k4 = h*pf(x+h,y+k3,z+l3);
	double l4 = h*pg(x + h, y + k3, z + l3);
	z=z+(1.0/6.0)*(l1 + 2*l2 + 2*l3 + l4);
	return y+(1.0/6.0)*(k1+2*k2+2*k3+k4);
}
