/* ODE RHS function abstract base class definitions.

   Claudia Castro Castro
   @ SMU
   Fall 2015 */

#ifndef ACTION_RHS_DEFINED__
#define ACTION_RHS_DEFINED__

// Inclusions
#include "matrix.hpp"

// Declare abstract base classes for ODE RHS and its Jacobian, to 
// define what the backward Euler solver expects from each.

//   Action of LB function abstract base class; derived classes 
//   must at least implement the Evaluate() routine


class RHSFunction {
 public: 
  //virtual int Evaluate(double alpha, Matrix &q0, Matrix &p0, Mat &LA) = 0;
  virtual int Evaluate(double delta, Matrix &alpha, Matrix &q0, Matrix &p0, Matrix &dqdtA, Matrix &dpdtA)=0;
};

//   Action of LB RHS function abstract base class; derived 
//   classes must at least implement the Evaluate() routine
class RHSFunctionCouplingB {
 public: 
  virtual int Evaluate(double delta, Matrix &q0, Matrix &p0, Matrix &dqdtB, Matrix &dpdtB)=0;
};

//   Action of LB RHS function abstract base class; derived 
//   classes must at least implement the Evaluate() routine
class RHSFunctionCouplingC {
 public: 
  virtual int Evaluate(double delta, Matrix &q0, Matrix &p0, Matrix &dqdtC, Matrix &dpdtC)=0;
};



#endif
