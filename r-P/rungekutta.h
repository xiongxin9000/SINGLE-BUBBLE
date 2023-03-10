typedef double (* func_ptr)( double x, double y,double z);


class RungeKutta {
	func_ptr pf;						// f(x,y,z)
	func_ptr pg;						// g(x,y,z)
	double h;							// step size
public:
	RungeKutta();							// default constructor  
	RungeKutta(func_ptr pf, func_ptr pg,double step); 	// alternate constructor
	// ~RungeKutta(); 	// deconstructor
	RungeKutta& setFunc(func_ptr pf,func_ptr pg);		// set the pointer
	RungeKutta& setStep(double step);		// set the step size
	double operator()(double x, double y,double &z) const;  // function call operator
};
