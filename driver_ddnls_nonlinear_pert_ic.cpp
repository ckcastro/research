/* Main routine to test ABC4 symplectic integrator method -- O(h^4) accurate 
   Symplectic Integrator solver method on the 1D-valued ODE problem 

 Claudia Castro-Castro
 October 1, 2016  @SMU
 Second order symplectic integrator ABC2 for disordered Discrete nonlinear
 Schrodinger (DDNLS) with periodic boundary conditions (ring-like array)

  i d/dz \psi_n = \psi_{n+1} + \psi_{n-1} + |\psi_n|^2 \psi_n

  ABC2(h) =
       exp(h/2LA)*exp(h/2 LB)*exp(hLC)*exp(h/2 LB)*exp(h/2 LA),

  where the Hamiltonian was expressed as the sum of three integrable parts A,B, and C 

           H(q,p) = A(q,p) + B(q,p) + C(q,p).

  and h is the time step.

 Reference: 
       Skokos, Ch, E. Gerlach, J. D. Bodyfelt, G. Papamikos, and S. Eggl. 
       "High order three part split symplectic integrators: Efficient techniques for the long  
       time simulation of the disordered discrete nonlinear Schrödinger equation." 
       Physics Letters A 378, no. 26 (2014): 1809-1815.

*/

#include "matrix.hpp"
#include "rhs.hpp"
#include "abc4.hpp"
#include <string>     // string, to_string
#include <iostream>   // std
#include <sstream>
#include <fstream>
#include <iomanip>
#include <cstdio>
#include <stdio.h>
#include <cstring>
#include <ctime>
#include <math.h>       // sqrt  
#include <cmath>        // pow  
#include <complex>      // complex
#include <random>
//#include <omp.h>        // open mp 

/* >> echo $OMP_NUM_THREADS
   >> env // enviroment variables
   >> export OMP_NUM_THREADS = 4  */

using namespace std;

// Define classes to compute action of operators 

// Action of operator Exp(L_A)
class MyRHS_LA: public RHSFunction {
public:
  //double delta;                              // stores some local data
  int Evaluate(double delta, Matrix &alpha, Matrix &q0, Matrix &p0, Matrix &dqdtA, Matrix &dpdtA) {    

  // evaluates the RHS function
    #pragma omp parallel for
    for (long int k=0; k < q0.Size(); k++) {
      dqdtA.data[0][k] = q0.data[0][k]*cos(alpha.data[0][k]*delta) + p0.data[0][k]*sin(alpha.data[0][k]*delta);
      dpdtA.data[0][k] = p0.data[0][k]*cos(alpha.data[0][k]*delta) - q0.data[0][k]*sin(alpha.data[0][k]*delta);
    }
   
    return 0; 
  }// end evaluate
};

// Action of operator Exp(L_B)
class MyRHS_LB: public RHSFunctionCouplingB {
public:
  //double delta;                              // stores some local data
  int Evaluate(double delta, Matrix &q0, Matrix &p0, Matrix &dqdtB, Matrix &dpdtB) {    
   // evaluates the RHS function

    dpdtB = p0;

    // coupling end point: n = 1
    dqdtB.data[0][0] = q0.data[0][0] + (p0.data[0][p0.Size()-1] + p0.data[0][1])*delta;

    // coupling interior points
    for(long int k = 1; k < p0.Size()-1; k++)
      dqdtB.data[0][k] =  q0.data[0][k] + (p0.data[0][k-1] + p0.data[0][k+1])*delta ;
    
    // coupling end point: n = N
    dqdtB.data[0][p0.Size()-1] = q0.data[0][p0.Size()-1] + (p0.data[0][p0.Size()-2] + p0.data[0][0])*delta;


    return 0; 
  }//end evaluate
};


// Action of operator Exp(L_C)
class MyRHS_LC: public RHSFunctionCouplingC {
public:
  //double delta;                              // stores some local data
  int Evaluate(double delta, Matrix &q0, Matrix &p0, Matrix &dqdtC, Matrix &dpdtC) {    // evaluates the RHS function

    dqdtC = q0;

    // coupling end points: n = 1
    dpdtC.data[0][0] = p0.data[0][0] - (q0.data[0][p0.Size()-1] + q0.data[0][1])*delta;

    // coupling interior points
    for(long int k = 1; k < p0.Size()-1; k++)
      dpdtC.data[0][k] =  p0.data[0][k] - (q0.data[0][k-1] + q0.data[0][k+1])*delta ;
 
    // coupling end points: n = N
    dpdtC.data[0][p0.Size()-1]= p0.data[0][p0.Size()-1] - (q0.data[0][p0.Size()-2] + q0.data[0][0])*delta;
    
    return 0; 
  } // end evaluate
};



// main routine
int main() {

  /////////////// Three part split symplectic integrator - abc 2 : 1d DDNLS/////////////////
  cout << " ====================================================================";
  cout << "\n\tThree part split symplectic integrator - abc 4 : 1d DDNLS \n";
  cout << " ====================================================================\n";

  //printf(" Current thread number: %d\n",omp_get_thread_num());
  clock_t t;

  // Print current date and time
  time_t result = time(0);
  cout << " " << asctime(std::localtime(&result)) << endl;

  // ====================================================================";
  // set problem information
  // ====================================================================";
  Matrix delta("0.05");     //step size

  double z0 = 0.0;         // initial time 	
  long double zfinal =2e3;   //final time  
  long double captures0 = 1e3;
  long double zsteps = zfinal/delta(0) ;    // integration streps
  Matrix zspan = Linspace(z0, zfinal, 1, zsteps);    // linear span constructor from matrix class

  double zcur;
  double dzout = delta(0);

  // Set size of optical fiber array
  // Node ODD number of fibers so we can have only one middle fiber
  long int numfib = 4096;
  long double Nc = ( numfib + 1 )/ 2.0;  // central fiber

  // Nonlinearity strength
  Matrix betas("1.0"); 

  // Disorder strength
  Matrix W("1.0");

  //No randomness
  //Matrix eps(numfib,1);

  // Uniformly distributed on-site energies between [ - W/2, W/2 ]
  // eps = -W/2+ ( W/2 + W/2 )*Rand(1, size)
  Matrix eps  = Random(numfib,1);
  
  eps.Mul( W(0) );
  eps.Add(- W(0)/2.0 ) ;
 
  Matrix alpha(numfib,1); 

  // Initializing all arrays
  double newsize0 = 0;
  newsize0 = zsteps/captures0 + 0.5;
  // add + 0.5 before casting to int, because the compiler will always truncate.
  //int newsize = 1;
  //newsize = newsize0+1;
  int newsize = (int)(newsize0); 
  newsize = newsize+1;
  int captures = (int)(captures0); 

  Matrix Serr(newsize,1);      // Norm
  Matrix Eerr(newsize,1);      // Energy
  Matrix zcapture(newsize,1);  // time
  Matrix m2_t(1,1);            // second moment
  Matrix m2(newsize,1);        // second moment
  Matrix par(newsize,1);       // participation
  Matrix zeros1(numfib,1);     // zeros
  Matrix zeros2(newsize,1);     // zeros
  Matrix coupling(numfib-1,1); 
  Matrix SumCoupling(1,1);
  Matrix SumE0(1,1); 
  Matrix E0(1,1);
  Matrix tmpE0(numfib,1);

  // ====================================================================";
  //	Set up initial conditions \n";
  //===================================================================="
 
  Matrix q0(numfib,1); 
  Matrix p0(numfib,1); 

  //Matrix Z0data("0.044194174");  //  sqrt(2/N), numfib = 1024 
  //Matrix Z0data("0.127000127");
  //Matrix Z0 = Z0data.T();

  //cout << " Z0 = " << Z0(0) <<endl;

   
  // Uniformly distributed on-site initial conditions between [ a, b ]
  // Z0 = a + ( b - a )*Rand(1, size)
  // [ 0.3 - eps, 0.3 + eps ]

  Matrix Z0  = Random(numfib,1);
  Z0.Mul(0.02);
  Z0.Add(0.3-0.01) ; // 0.3/sqrt(2)
  
  // UNIFORM ALONG ALL THE ARRAY( REAL PART - q0)
  for (long int il=0; il<numfib; il++)
  	//q0(il,0) = il ;
  	//q0(il,0) = sqrt(2.0)*Z0(il) ;
        q0(il,0) = Z0(il) ;

  // store initial profiling amplitudes
  q0.Write("q0.dat"); // See diary_abc4.txt 
  p0.Write("p0.dat");


  
  // ====================================================================";
  // Screen output 
  // ====================================================================";


  q0.Trans();
  cout << " q0 = " << q0 << endl;//<< " p0 = " << p0 << endl;
  q0.Trans();


  cout << " initial time = "<< z0 << "\tfinal time = " << zfinal << endl;
  cout << " step size = " << delta(0) << "\tsteps = "<< zsteps <<endl;
  cout << " number of fibers = " << numfib << "   middle fiber = " << Nc << endl;
  eps.Trans();
  cout << " eps first site = " << eps.data[0][0] << endl;
  eps.Trans();
  cout << " disorder strength: W = "<< W(0) << endl;
  cout << " random on-site energies chosen uniformly distributed from interval [ -W/2, W/2 ] "<< endl;
  cout << " \n captures double: " << captures0 << endl;
  cout << " captures int : " << captures << endl;
  cout << " newsize = " << newsize << endl;
  cout << " Serr size " << Serr.Rows() << "Rows x " << Serr.Cols() <<"Cols " <<endl;

  // Write data to files 
  betas.Write("betas_file.dat"); 
  delta.Write("delta_file.dat");
  eps.Write("eps.dat"); // Write to a file to test Matlab code.
  ofstream zffile;
  zffile.open ("zfinal_file.dat"); zffile << zfinal ; zffile.close();
  ofstream numfibfile;
  numfibfile.open ("numfib_file.dat"); numfibfile << numfib ; numfibfile.close();
 
  // ====================================================================";
  //	Scheme coefficients for ABC4[Y] \n";
  //===================================================================="
  Matrix x0(1,1);
  Matrix x1(1,1);

  x0 = - pow(2.0, 1.0/3.0)/(2.0 - pow(2.0,1.0/3.0));
  x1 =  1.0/(2.0 - pow(2.0,1.0/3.0));

  // ====================================================================";
  //	Conserved quantities \n";
  //===================================================================="

  // define Initial Solution Norm as 1/2 (p0^2 + q0^2)
  Matrix A0(q0); 
  A0.Mul(q0);       // q0^2 
  Matrix tmpp0(p0); // tmp = p0
  tmpp0.Mul(tmpp0); // tmp = p0^2
  A0.Add(tmpp0);    // q0^2 +p0^2
  A0.Mul(1.0/2.0);  //(q0^2 +p0^2)/2

  //cout << "A0 = "<< A0 << endl;


  // Initial norm
  Matrix S0(1,1);

  S0 =  Sum(A0);


  // Store spatio-temporal info in array
  Matrix Normdata(numfib, newsize); 
  
  // initial  norm
  for (long int il=0; il<numfib-1; il++)
    Normdata(il,0) = A0(il,0) ;

  cout <<" Norm data " << Normdata.Rows() << "Rows x " << Normdata.Cols() <<"Cols " <<endl;
  cout << " A0 size:  "<< A0.Size() << endl;


  // =====================================================================
  // Second moment initializations
  // =====================================================================
  Matrix ll = Linspace(1, numfib, 1, numfib);
  ll.Size();
  ll.Trans(); 

  Matrix lbar(numfib,1);
  //#pragma omp parallel for
  for(long int il = 0; il < ll.Size(); il++ ){
    lbar.data[0][il] = (ll.data[0][il]- Nc)*(ll.data[0][il] - Nc);
  }
  //lbar.Trans(); cout << lbar << endl; lbar.Trans();

  // =====================================================================
  // create LA, LB, and LC operator RHS objects
  // =====================================================================
  MyRHS_LA rhsLA;
  MyRHS_LB rhsLB;
  MyRHS_LC rhsLC;
  
  // create abc2 stepper object
  abc4_stepper ABC4(rhsLA, rhsLB, rhsLC, q0, p0);

  // temporary variables
  Matrix A(numfib,1); // Norm (p^2 +q^2)/2
  Matrix tmpE(A);
  Matrix Z(1,1), tspan(1,2);
  Matrix S(1,1); // Norm error
  Matrix E(1,1); // Energy error




    long int istep;
    long int icapture;

    // loop over beta values
    for (int ib = 0; ib < betas.Size(); ib++) {
      cout << "\n    nonlinearity = " << betas(ib);

      //clock function to meassure cpu time for every beta
      t = clock();

      // set up the initial time
      zcur = z0;

      // reset diagnostics to zero
      Eerr = zeros2;
      Serr = zeros2;
      m2  = zeros2;

      /* define Total Initial Energy (Hamiltonian)
      H = \sum_{n}\epsilon_{n}|\psi_{n}|^{2} + \frac{\gamma}{2}|\psi_{n}|^{4}...
                        + c(\psi_{n+1}\psi_{n}^{*} + \psi_{n+1}^{*}\psi_{n})
      */ 
      tmpE0 = A0; 

      // coupling 
      for (long int il=0; il<numfib-1; il++)
        coupling.data[0][il] = p0.data[0][il+1]*p0.data[0][il] + q0.data[0][il+1]*q0.data[0][il];

      tmpE0.Mul(betas(ib)/2.0); 
      tmpE0.Add(eps); 
      tmpE0.Mul(A0);
      SumE0 = Sum(tmpE0);
      SumCoupling =  Sum(coupling);
      E0 = SumE0 + SumCoupling ;
      cout << "\n\t Initial Norm : S0 = sum((q0^2 +p0^2)/2) = " << S0 << endl;
      cout << "\t Initial Energy : E0 = " << E0 << endl;

      // reset counters
      istep = 0;
      icapture = 1;
    
      // Call solver and output error
      while (zcur < 0.9999999999999*zfinal) {
          
        // set the time interval for this solve
        tspan(0) = zcur;
        tspan(1) = zcur + dzout;
        if (tspan(1) > zfinal)  tspan(1) = zfinal;


        // call the solver, update current time
        Matrix zvals = ABC4.Evolve(tspan, eps, alpha, delta(0), x0, x1, betas(ib), q0, p0);

        zcur = zvals(zvals.Size()-1);   // last entry in tvals
        

        //cout << "istep congruent? captures : " << istep % captures<<endl; 
        if ( istep % captures == 0.0 )
        {

          // define Solution Norm at current time
          //#pragma omp parallel for
          for (long int il=0; il<numfib-1; il++)
            A.data[0][il] = p0.data[0][il]*p0.data[0][il] + q0.data[0][il]*q0.data[0][il];

          A.Mul(1.0/2.0);

          
          cout << " istep = " << istep <<"  " ;
          cout << " icapture = " << icapture <<"  " ;

          if (istep < zsteps-1){ 
            for (int il = 0; il < numfib-1; il++)
              Normdata(il,icapture-1) = A(il,0) ;    
          }

          S =  Sum(A);

          //cout << " S = " << S(0);

          // define relative Norm error at current time
          S.Add(-S0(0));    //cout << " Sr -> Sum(A)-S0" << Sr;
          S.Mul(1.0/S0(0)); //cout << " Sr -> (Sum(A)-S0)/S0" << Sr;
          S.Abs();          //cout << " Sr ->" << Sr << " i = " << icapture << endl;
          Serr(icapture-1,0) = S(0);
 
          cout <<  "   Serr(" << icapture-1 << ") = " << Serr(icapture-1,0);      
          //cout << " Serr vector" << Serr ;

          // define Total Energy (Hamiltonian) at current time
          for (long int il=0; il<numfib-1; il++)
            coupling.data[0][il] = p0.data[0][il+1]*p0.data[0][il] + q0.data[0][il+1]*q0.data[0][il];

          SumCoupling =  Sum(coupling);  
          tmpE = A; 
          tmpE.Mul(betas(ib)/2.0); // cout << endl <<  " beta/2(q^2 + p^2)/2= " ; tmpE.Write();    
          tmpE.Add(eps); //cout << endl <<  " eps + beta/2(q^2 + p^2)/2= " ; tmpE.Write();
          tmpE.Mul(A);  //cout << endl <<  " (q^2 + p^2)/2(eps + beta/2(q^2 + p^2)/2)= "; tmpE.Write();         
          E = Sum(tmpE) + SumCoupling ; //cout << " Energy = Sum(tmpE) + SumCoupling =  " << E << endl;
          cout << "  " << " E = " << E(0) ;

          // Relative energy error
          E.Add( - E0(0) );  //cout << " Er =   E - E0" << E;
          E.Mul(1.0/E0(0));  //cout << " Er =  (E - E0)/E0" << E;
          E.Abs();           //cout << " Er = |(E - E0)/E0|" << E << " i = " << icapture << endl;
          Eerr.data[0][icapture-1] = E.data[0][0];

          cout << "  " << " Eerr(" << icapture-1 <<") = " << Eerr(icapture-1,0)<<endl;

          zcapture.data[0][icapture-1] = zvals(1);
 
          icapture++;

          //cout << " icapture = " << icapture <<"\t" ;

        } // end if captures

        istep++;

      } //end while 

      
      t = clock() - t;
      printf ("\n\t Elapsed time %f seconds.\n",((float)t)/CLOCKS_PER_SEC);


      char Serrfilename[80]; 
      char Eerrfilename[80]; 
      char finalNormfilename[80]; 
      char numstrbetas[21]; // enough to hold all numbers up to 64-bits
      sprintf(numstrbetas, "%g", betas(ib));
      char numstrzf[21]; // enough to hold all numbers up to 64-bits
      sprintf(numstrzf, "%Lf", zfinal);

      strcpy (Serrfilename,"Serr_abc4_beta_"); strcat (Serrfilename, numstrbetas); strcat (Serrfilename, ".dat");
      Serr.Write(Serrfilename);

      strcpy (Eerrfilename,"Eerr_abc4_beta_"); strcat (Eerrfilename, numstrbetas); strcat (Eerrfilename, ".dat");
      Eerr.Write(Eerrfilename);

      // Norm distribution at final time
      strcpy (finalNormfilename,"Norm_abc4_beta_"); strcat (finalNormfilename, numstrbetas); 
      strcat (finalNormfilename, "_zf"); 
      //strcat (finalNormfilename, numstrzf); 
      strcat (finalNormfilename, ".dat");
      A.Write(finalNormfilename);

      // print date and time at the end
      time_t result2 = time(0);
      cout << endl << asctime(std::localtime(&result2)) << endl;

      // print spatio-temporal norm info
      Normdata.Write("Normdata.dat");
 
    } // end for over betas

    // Write data to files 
    zcapture.Write("timefile.dat"); // timefile

  return 0;
} // end main
