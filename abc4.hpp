/* ABC2 symplectic integrator time stepper class header file.

   Claudia Castro-Castro
   1D-DNLS @ SMU
   Fall 2015  */

#ifndef ABC4_DEFINED__
#define ABC4_DEFINED__

// Inclusions
#include <math.h>
#include "matrix.hpp"
#include "rhs.hpp"

// Declare abstract base class for action RHSs, to define what the 
// abc2_stepper expects; derived classes must at least 
// implement the Evaluate() routine


// abc2 time stepper class
class abc4_stepper {

 private:

  // private reusable local data
  Matrix *dqdtA;             // storage for q vector
  Matrix *dpdtA;             // storage for p vector
  RHSFunction *rhsLA;        // pointer to action RHS function

  Matrix *dqdtB;                  // storage for q vector
  Matrix *dpdtB;                  // storage for p vector
  RHSFunctionCouplingB *rhsLB;     // pointer to  action RHS function

  Matrix *dqdtC;                  // storage for q vector
  Matrix *dpdtC;                  // storage for p vector
  RHSFunctionCouplingC *rhsLC;     // pointer to  action RHS function

 public:

  // constructor (sets RHS function pointer, copies y for local data)
  abc4_stepper(RHSFunction &rhsLA_, RHSFunctionCouplingB &rhsLB_, RHSFunctionCouplingC &rhsLC_, Matrix &q0, Matrix &p0) {
    rhsLA = &rhsLA_;
    rhsLB = &rhsLB_;
    rhsLC = &rhsLC_;

    dqdtA = new Matrix(q0);
    dpdtA = new Matrix(p0);

    dqdtB = new Matrix(q0);
    dpdtB = new Matrix(p0);

    dqdtC = new Matrix(q0);
    dpdtC = new Matrix(p0);
  };

  // destructor (frees local data)
  ~abc4_stepper() { delete dqdtA; delete dpdtA; delete dqdtB; delete dpdtB; delete dqdtC; delete dpdtC;};


  // Evolve routine (evolves the solution via abc2)
  Matrix Evolve(Matrix zspan, Matrix &eps, Matrix &alpha, double delta, Matrix &x0, Matrix &x1, double beta, Matrix &q0, Matrix &p0);

};

#endif
