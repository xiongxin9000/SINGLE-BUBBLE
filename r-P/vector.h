#ifndef VECTOR_H //Include guard
#define VECTOR_H

#include <iostream> //Generic IO operations
#include <fstream>  //File IO operations
#include <stdexcept> //provides exceptions

class Vector {
private:
    int num; // Number of elements
    double* pdata; // Pointer to the data
    void Init(int Num); // private function since user should not call it
public:
	// CONSTRUCTORS
    Vector(); // default constructor
    //alternate constructor, it is explicit since implicit
    //type conversion int -> vector doesn't make sense
    explicit Vector(int Num); //construct vector of given size
    Vector(const Vector& v); // copy constructor

	// DESTRUCTOR
    virtual ~Vector(); // destructor is virtual since it is
    //a good thing to have virtual destr. in base class

	// SIZE
    int size() const; //get number of elements in the vector

	// OVERLOADED OPERATORS
    Vector& operator= (const Vector& v); //overloaded assignment operator
    double& operator[] (int i); // overloaded array access operator for writing
    double operator[] (int i) const; // overloaded array access operator for reading

	// COMPARISON
	bool equiv(const Vector& v, double tol) const;

	// THREE NORMS
	double one_norm() const;
	double two_norm() const;
	double uniform_norm() const;

	// KEYBOARD/SCREEN INPUT AND OUTPUT
    friend std::istream& operator>>(std::istream& is, Vector& v);// keyboard input
    friend std::ostream& operator<<(std::ostream& os, const Vector& v);// screen output


	// FILE INPUT AND OUTPUT
    //the file output operator is compatible with file input operator,
    //ie. everything written can be read later.
    friend std::ifstream& operator>>(std::ifstream& ifs, Vector& v);// file input
    friend std::ofstream& operator<<(std::ofstream& ofs, const Vector& v);// file output
};

#endif