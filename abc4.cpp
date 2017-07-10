/* abc4_stepper class implementation file.

 abc_stepper
   This function calculates one time step using a second order symplectic
   integrator algorithm abc4, Forth order Symplectic Method of 3 stages and 13 steps

           ABC4(delta) = ABC2(x1 delta) ABC2(x0 delta) ABC2(x1 delta),

  where

  ABC2(delta) =
       exp(delta/2LA)*exp(delta/2LB)*exp(deltaLC)*exp(delta/2LB)*exp(delta/2LA),

  where the Hamiltonian was expressed as the sum of three integrable parts A,B, and C 

           H(q,p) = A(q,p) + B(q,p) + C(q,p).

   OUTPUT:
       q: evolution of generalized coordenates after one time step.
       p: evolution of momenta after one time step.
   INPUT:
       q0: initial condition of generalized coordinate
       p0: initial condition of momenta
       eps: on site energies
       beta : nonlineatity strength
       delta:  time step
       dqdtA : propagation of initial condition, q, given by the operator exp(delta LA)
       dpdtA : propagation of initial condition, p, given by the operator exp(delta LA)
       dpdtB : propagation of initial condition, p, given by the operator exp(delta LB)
       dqdtB : propagation of initial condition, q, given by the operator exp(delta LB)
       dqdtC : propagation of initial condition, q, given by the operator exp(delta LC)
       dpdtC : propagation of initial condition, p, given by the operator exp(delta LC)


       The return value is a row vector containing all internal 
       times at which the solution was computed,
                   [z0, z1, ..., zN]

   Claudia Castro-Castro
   @ SMU
   Fall 2015  */

#include "matrix.hpp"
#include "abc4.hpp"



Matrix abc4_stepper::Evolve(Matrix zspan, Matrix &eps, Matrix &alpha, double delta, Matrix &x0, Matrix &x1, double beta, Matrix &q0, Matrix &p0) {

  // check for legal inputs 
  if (delta <= 0.0) {
    std::cerr << "abc4_stepper: Illegal step size\n";
    return zspan;
  }
  if (zspan(1) <= zspan(0)) {
    std::cerr << "abc4_stepper: Illegal zspan\n";
    return zspan;	  
  }
  
  // figure out how many time steps
  long int N = (zspan(1)-zspan(0))/delta;
  if (zspan(1) > zspan(0)+N*delta)  N++;
    
  // create ouput Mat
  Matrix times(1, N+1);
  times(0) = zspan(0);

  // iterate over time steps
  for (int i=0; i<N; i++) {

    // last step only: update delta to stop directly at final time
    if (i == N-1) 
      delta = zspan(1)-times(i);
    

    // *******************************************************************


    //delta0 = (x0(0)*delta)/2.0;
    //delta1 = (x1(0)*delta)/2.0;
    //delta01 = ((x0(0) + x1(1))*delta)/2.0;

    // compute action of exp(LA) over current sol - 1st step
    // Compute alpha 
    for (long int i=0; i<p0.Size(); i++)
      alpha.data[0][i] = eps.data[0][i] + (0.5)*(beta)*(p0.data[0][i]*p0.data[0][i] + q0.data[0][i]*q0.data[0][i]);
 
 
    if (rhsLA->Evaluate((x1.data[0][0]*delta)/2.0, alpha, q0, p0, *dqdtA, *dpdtA)  != 0) {
      std::cerr << "abc4_stepper: Error in operator action A RHS function\n";
      return times;
    }
    
    // update solution  
    // compute action of exp(LB) over current sol - 2nd step
    if (rhsLB->Evaluate((x1.data[0][0]*delta)/2.0, *dqdtA, *dpdtA , *dqdtB, *dpdtB)  != 0) {
      std::cerr << "abc4_stepper: Error in operator action B RHS function\n";
      return times;
    }


    // compute action of exp(LC) over current sol - 3rd step
    if (rhsLC->Evaluate(x1.data[0][0]*delta, *dqdtB, *dpdtB , *dqdtC, *dpdtC)  != 0) {
      std::cerr << "abc4_stepper: Error  in operator action C RHS function\n";
      return times;
    }


    // compute action of exp(LB) over current sol - 4th step
    if (rhsLB->Evaluate((x1.data[0][0]*delta)/2.0, *dqdtC, *dpdtC , *dqdtB, *dpdtB)  != 0) {
      std::cerr << "abc4_stepper: Error in operator action B RHS function\n";
      return times;
    }

    // compute action of exp(LA) over current sol - 5th step 
    // Modify alpha
    for (long int i=0; i<q0.Size(); i++)
      alpha.data[0][i] = eps.data[0][i] + (0.5*beta)*(   (*dpdtB).data[0][i]*(*dpdtB).data[0][i] 
                                                       + (*dqdtB).data[0][i]*(*dqdtB).data[0][i]  );



    if (rhsLA->Evaluate( ((x0.data[0][0] + x1.data[0][0])*delta)/2.0, alpha, *dqdtB, *dpdtB, *dqdtA, *dpdtA)  != 0) {
      std::cerr << "abc4_stepper: Error in operator action A RHS function\n";
      return times;
    }

    // compute action of exp(LB) over current sol - 6th step
    if (rhsLB->Evaluate((x0.data[0][0]*delta)/2.0, *dqdtA, *dpdtA , *dqdtB, *dpdtB)  != 0) {
      std::cerr << "abc4_stepper: Error in operator action B RHS function\n";
      return times;
    }

    // compute action of exp(LC) over current sol  - 7th step
    if (rhsLC->Evaluate(x0.data[0][0]*delta, *dqdtB, *dpdtB , *dqdtC, *dpdtC)  != 0) {
      std::cerr << "abc4_stepper: Error  in operator action C RHS function\n";
      return times;
    }

    // compute action of exp(LB) over current sol - 8th step
    if (rhsLB->Evaluate((x0.data[0][0]*delta)/2.0, *dqdtC, *dpdtC , *dqdtB, *dpdtB)  != 0) {
      std::cerr << "abc4_stepper: Error in operator action B RHS function\n";
      return times;
    }

    // compute action of exp(LA) over current sol - 9th step
    // Modify alpha
    for (long int i=0; i<q0.Size(); i++)
      alpha.data[0][i] = eps.data[0][i] + (0.5*beta)*(   (*dpdtB).data[0][i]*(*dpdtB).data[0][i] 
                                                       + (*dqdtB).data[0][i]*(*dqdtB).data[0][i]  );

 
    if (rhsLA->Evaluate( ((x0.data[0][0] + x1.data[0][0])*delta)/2.0, alpha, *dqdtB, *dpdtB, *dqdtA, *dpdtA)  != 0) {
      std::cerr << "abc4_stepper: Error in operator action A RHS function\n";
      return times;
    }

    // compute action of exp(LB) over current sol - 10th step
    if (rhsLB->Evaluate((x1.data[0][0]*delta)/2.0, *dqdtA, *dpdtA , *dqdtB, *dpdtB)  != 0) {
      std::cerr << "abc4_stepper: Error in operator action B RHS function\n";
      return times;
    }

    // compute action of exp(LC) over current sol - 11th step
    if (rhsLC->Evaluate(x1.data[0][0]*delta, *dqdtB, *dpdtB , *dqdtC, *dpdtC)  != 0) {
      std::cerr << "abc4_stepper: Error  in operator action C RHS function\n";
      return times;
    }

    // compute action of exp(LB) over current sol - 12th step
    if (rhsLB->Evaluate((x1.data[0][0]*delta)/2.0, *dqdtC, *dpdtC , *dqdtB, *dpdtB)  != 0) {
      std::cerr << "abc4_stepper: Error in operator action B RHS function\n";
      return times;
    }


    // compute action of exp(LA) over current sol - 13th step
    // Modify alpha
    for (long int i=0; i<q0.Size(); i++)
      alpha.data[0][i] = eps.data[0][i] + (0.5*beta)*(   (*dpdtB).data[0][i]*(*dpdtB).data[0][i] 
                                                       + (*dqdtB).data[0][i]*(*dqdtB).data[0][i]  );
    if (rhsLA->Evaluate( (x1.data[0][0]*delta)/2.0, alpha, *dqdtB, *dpdtB, *dqdtA, *dpdtA)  != 0) {
      std::cerr << "abc4_stepper: Error in operator action A RHS function\n";
      return times;
    }

    // perform update
    q0 = (*dqdtA);
    p0 = (*dpdtA);

    // update current time, store in output array
    times(i+1) = times(i) + delta;
  }

  return times;
}
