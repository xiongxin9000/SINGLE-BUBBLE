#include <iomanip>
#include <cmath>
#include "vector.h"
#include "rungekutta.h"
#include <fstream>
//reference: KirillGordievich/The-4th-order-Runge-Kutta-method-for-a-2nd-order-ODE
double f1(double t, double R,double z);
double f2(double t, double R,double z);

//double soln1(double x);

using namespace std;
char filename[]="setupfile.txt";
int myindex=0;
double S;//suface tension 
double R0;//initial bubble radius S/2R=0.000212  p_inf-p_v=0.00212 p_inf=0.00212+0.000692425=0.002812425
double P_inf;//boundary pressure
//P_inf[100]={0.00640246, 0.00640246, 0.00640246, 0.00640246, ...};//boundary pressure, time dependent
double rho_l;//liquid density
double nu_l;//liquid kinematic viscosity
//double P0=4.75623E-06;//initial liquid pressure
double P_v;//bubble pressure
double r_inf;//boundary distance 
// int SIZE = 500;
// double P_150[SIZE];
void Readinitialfile(double* S,double* R0,double *P_inf,double* rho_l,double* nu_l,double* P_v,double* r_inf
,double* z,double* time_start,double* time_end,double* step_size)
{ 
    FILE *file=NULL;
    file=fopen(filename,"r");
    if(!file)
    {
        std::cout<<"cannot find file"<<std::endl;
    }
    fscanf(file,"%lf",S);
    fscanf(file,"%lf",R0);
    fscanf(file,"%lf",P_inf);
    fscanf(file,"%lf",rho_l);
    fscanf(file,"%lf",nu_l);
    fscanf(file,"%lf",P_v);
    fscanf(file,"%lf",r_inf);
    fscanf(file,"%lf",z);
    fscanf(file,"%lf",time_start);
    fscanf(file,"%lf",time_end);
    fscanf(file,"%lf",step_size);


    
    fclose(file);
}
// void readData() {


//     string inFileName = "pressure.txt";
//     ifstream inFile;
//     inFile.open(inFileName.c_str());

//     if (inFile.is_open())
//     {
//         for (int i = 0; i < SIZE; i++)
//         {
//             inFile >> P_150[i];
//             //std::cout << P_150[i] << " "<<"\n";
//         }

//         inFile.close(); // CLose input file
//     }
//     else { //Error message
//         std::cerr << "Can't find input file " << inFileName << endl;
//     }
// }
std::string to_string_with_precision(const double a_value, const int n = 2)
{
    std::ostringstream out;
    out.precision(n);
    out << std::fixed << a_value;
    return out.str();
}
int main()
{
	double R,z;
	double t_start, t_end, step;

	// readData();
	Readinitialfile(&S,&R0,&P_inf,&rho_l,&nu_l,&P_v,&r_inf,&z,&t_start,&t_end,&step);
	try {
		// std::cout << "input start value for R\n";//1
		// std::cin >> R;

		// std::cout << "input start value for z\n";//2
		// std::cin >> z;

		// std::cout << "input the range for t\n";//0,1
		// std::cin >> t_start >> t_end;

		// std::cout << "input step size\n";
		// std::cin >> step;
		std::cout<<"S="<<S<<endl;
		std::cout<<"R0="<<R0<<endl;
		std::cout<<"P_inf="<<P_inf<<endl;
		std::cout<<"rho_l="<<rho_l<<endl;
		std::cout<<"nu_l="<<nu_l<<endl;
		std::cout<<"P_v="<<P_v<<endl;
		std::cout<<"r_inf="<<r_inf<<endl;
		std::cout<<"z="<<z<<endl;
		std::cout<<"t_start="<<t_start<<endl;
		std::cout<<"t_end="<<t_end<<endl;
		std::cout<<"step="<<step<<endl;
		// calculate the number of steps
		int n = (int)((t_end - t_start) / step);
		R=R0;
		// create vectors to store the y values and the x values 
		Vector v(n + 1);
		Vector t_step(n + 1);
		Vector vz(n+1);

		// initialise the start values
		t_step[0] = t_start;
		v[0] = R;
		vz[0] =z;	
		// perform the iteration
		for (int i = 1; i <= n; i++)
		{
			myindex=i-1;
			RungeKutta r(f1,f2, step);//pf pg
			// construct the runge kutta object
			R = r(t_step[i - 1], R,z);//operation (x,y,z)
			v[i] = R;
			vz[i]= z;
			t_step[i] = t_step[i - 1] + step;
		}

		// ouput the results
		// std::cout << "Solution to infinite P-R, " << t_start << " <= x  <= " << t_end << ", step size " << step << "\n";
		// std::cout << "**************************************************\n";
		// std::cout << "t " << std::setw(30) << "R(x_n)" <<"\n";
		double rounded_rho = std::round(rho_l * 100) / 100;
		// std::cout<<rounded_rho<<std::endl;
		std::string outputfile="R-P_"+std::to_string((int)(r_inf*2))+"_"+to_string_with_precision(rounded_rho)+".dat";
		std::ofstream out(outputfile);
		for (int i = 0; i <= n; i++) {
			out << std::scientific << std::setprecision(7) << std::setw(15) << t_step[i];
			out << std::scientific << std::setprecision(7) << std::setw(15) << v[i];
			out << std::endl;
		}
		out.close();
	} catch (std::exception& e) {
		// Catching other errors
		std::cerr << "std::exception caught" << std::endl;
		std::cerr << "Type: " << typeid(e).name() << std::endl;
		std::cerr << "What: " << e.what() << std::endl;
	}
	return 0;
}

// the functions to test
double f1(double t, double R,double z) { return z; }//pf
//double f2(double t, double R,double z) { return (P_B-P_inf)/(rho_l*R)-3*z*z/(2*R)-4*nu_l*z/(R*R)-2*S/(rho_l*R*R); }
//1.Numerical analysis of Rayleigh–Plesset equation for cavitating water jets
//2.Dynamics of an acoustically driven cavitation bubble cluster in the vicinity of a solid surface
// double f2(double t, double R,double z) { return (P_v-P_inf+(P0+2*S/R0-P_v)*pow((R0/R),3))/(rho_l*R)-3*z*z/(2*R)-4*nu_l*z/(R*R)-2*S/(rho_l*R*R); }
// double f2(double t, double R,double z) { return (-0.182876087)/(rho_l*R*(1-R/50))-z*z*(1.5-2*R/50+pow(R,4)/(2*pow(50,4)))/((1-R/50)*R); }


//A numerical study of the early-stage dynamics of a bubble cluster
// double f2(double t, double R,double z) 
// { return ((P_v-P_inf)/rho_l-4*nu_l*z/R-2*S/(rho_l*R)+z*z/2)/((1/R-1/r_inf)*R*R)-2/R*z*z;}

//Dynamics of an acoustically driven cavitation bubble cluster in the vicinity of a solid surface
// double f2(double t, double R,double z) { return (P_v-2*S/R-4*nu_l*z/R-P_inf)/(rho_l*R)*r_inf/(r_inf-R)
// -(1.5-2*R/r_inf+pow(R,4)/(2*pow(r_inf,4)))*r_inf/(R*(r_inf-R))*z*z; }

//simulation of multiple cavitation bubbles interaction with single-component multiphase lattice boltzmann method
double f2(double t, double R,double z) { return ((P_v-P_inf-S/R-2*nu_l*z/R)/rho_l+(1-pow((R/r_inf),2))/2*z*z-log(r_inf/R)*z*z)/(log(r_inf/R)*R); }
// double f2(double t, double R,double z) {
// 	return ((P_150[myindex]-P_inf-S/R-2*nu_l*z/R)/rho_l+(1-pow((R/r_inf),2))/2*z*z-log(r_inf/R)*z*z)/(log(r_inf/R)*R);
// }


// the analytic solutions
//double soln1(double x) { return 6.0/5.0*cos(2.0*x) - 1.0/5.0*cos(3.0*x) + sin(2.0*x); } 