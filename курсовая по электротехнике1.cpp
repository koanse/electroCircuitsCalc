#include <complex>
#include <iostream>
#include <conio.h>
#include <stdio.h>


int main( )
{
   using namespace std;
   double pi = 3.14159265359;
   double E1m = 125, E2m = 100, E1 = E1m/sqrt(2.0), E2 = E2m/sqrt(2.0);
   double Jkm = 12, Jk = Jkm/sqrt(2.0);
   double psi[3] = {225.0/(2*pi), 10.0/(2*pi), -70.0/(2*pi)};
   double f = 50;
   double R[4] = {7, 7, 5, 6};
   double L[3] = {0.040, 0.075, 0.032};
   double C[3] = {0.000095, 0.000550, 0.000450};

   double xC[3], xL[3];
 
   complex <double> one ( 1.0 , 0.0 );
   complex <double> j ( 0.0 , 1.0 );
   complex <double> Z[6];
   complex <double> coFi[4];
   complex <double> Fi5, Fi6[2] ,Fi7[2] , Fi8[2], Fi9[2], Fi10[2], Fi11[2];
   complex <double> Y11, Y44, Y13, Y14, Y43;
   complex <double> Iuz1, Iuz4;
   complex <double> coE1, coE2, coJk;
   complex <double> delta, deltaFi1, deltaFi4;
   complex <double> coI[6];
   complex <double> Sist ( 0.0 , 0.0 ), Spr ( 0.0, 0.0 );
   
   for(int i = 0; i < 3; i++)
   {
	   xC[i] = 1/(2*pi*f*C[i]);
	   xL[i] = 2*pi*f*L[i];
	   printf("xC%d = %.5lf\n", i + 1, xC[i]);
	   printf("xL%d = %.5lf\n", i + 1, xL[i]);
   }

   coE1 = complex <double> (polar(sqrt(E1), psi[0]));
   coE2 = complex <double> (polar(sqrt(E2), psi[1]));
   coJk = complex <double> (polar(sqrt(Jk), psi[2]));
   
   coFi[1] = complex <double> (0.0, 0.0);
   coFi[2] = coE2;
   
   Z[0] = complex <double> (999999999.0, 0.0);
   Z[1] = complex <double> (R[2], xL[2] - xC[2]);
   Z[2] = complex <double> (0.0, 0.0);
   Z[3] = complex <double> (R[1], xL[0]);
   Z[4] = complex <double> (R[3], -xC[1]);
   Z[5] = complex <double> (0.0, xL[1] - xC[0]);

   Y11 = one/Z[5] + one/Z[1];
   Y44 = one/Z[3] + one/Z[5] + one/Z[4];
   Y13 = complex <double> (0.0, 0.0);
   Y14 = -one/Z[5];
   Y43 = -one/Z[3];

   cout << "Y11 = " << Y11 << endl;
   cout << "Y44 = " << Y44 << endl;
   cout << "Y13 = " << Y13 << endl;
   cout << "Y14 = " << Y14 << endl;
   cout << "Y43 = " << Y43 << endl;

   Iuz1 = coE1/Z[5] + coJk;
   Iuz4 = -one*coE1/Z[5];

   cout << "Iuz1 = " << Iuz1 << endl;
   cout << "Iuz4 = " << Iuz4 << endl;

   delta = Y11*Y44 - Y14*Y14;
   deltaFi1 = Iuz1*Y44 - (Iuz4-coFi[2]*Y43)*Y14;
   deltaFi4 = Y11*(Iuz4-coFi[2]*Y43) - Y14*Iuz1;

   coFi[0] = deltaFi1/delta;
   coFi[3] = deltaFi4/delta;

   coI[0] = coJk;
   coI[1] = (coFi[0] - coFi[1])/Z[1];
   coI[2] = (coFi[1] - coFi[2])/Z[2];
   coI[3] = (coFi[3] - coFi[2])/Z[3];
   coI[4] = (coFi[1] - coFi[3])/Z[4];
   coI[5] = (coFi[3] - coFi[0] + coE1)/Z[5];

   coI[2] = coI[0] - coI[3];
   
   cout << "Fi1 = " << coFi[0] << endl;
   cout << "Fi2 = " << coFi[1] << endl;
   cout << "Fi3 = " << coFi[2] << endl;
   cout << "Fi4 = " << coFi[3] << endl;

   for(i = 0; i < 6; i++)
	   cout << "I" << i+1 << " = " << coI[i] << endl;
   
// ѕотенциалы двум€ способами
   Fi5 = coFi[0] + coI[0]*R[0];
   
   Fi6[0] = coFi[0] - coI[1]*j*xL[2];
   Fi7[0] = Fi6[0] - coI[1]*R[2];

   Fi7[1] = coFi[1] - coI[1]*j*xC[2];
   Fi6[1] = Fi7[1] + coI[1]*R[2];

   Fi8[0] = coFi[2] + coI[3]*R[1];
   Fi8[1] = coFi[3] - coI[3]*j*xL[0];

   Fi9[0] = coFi[1] - coI[4]*R[3];
   Fi9[1] = coFi[3] - coI[4]*j*xC[1];

   Fi10[0] = coFi[3] + coI[5]*j*xC[0];
   Fi11[0] = Fi10[0] - coI[5]*j*xL[1];

   Fi11[1] = coFi[0] - coE1;
   Fi10[1] = Fi11[1] + coI[5]*j*xL[1];

 
// округление комплексных токов до 3-х знаков
 //  for(i = 0; i < 6; i++)
//	   coI[i] = complex <double> ((int(real(coI[i])*1000.0))/1000.0,
//										(int(imag(coI[i])*1000.0))/1000.0);


   Sist = coE1*conj(coI[5]) + coE2*conj(coI[2]) + conj(coJk)*(Fi5 - coFi[2]);
   
   for(i = 1; i < 6; i++)
	   Spr += abs(coI[i])*abs(coI[i])*Z[i];
   Spr += abs(coJk) * abs(coJk) * R[0];

      
   cout << "Sist = " << Sist << endl;
   cout << "Spr  = " << Spr << endl;
  
   cout << "Fi5  = " << Fi5 << endl;
   cout << "Fi6  = " <<  Fi6[0] << " and " << Fi6[1] << endl;
   cout << "Fi7  = " <<  Fi7[0] << " and " << Fi7[1] << endl;
   cout << "Fi8  = " <<  Fi8[0] << " and " << Fi8[1] << endl;
   cout << "Fi9  = " <<  Fi9[0] << " and " << Fi9[1] << endl;
   cout << "Fi10 = " <<  Fi10[0] << " and " << Fi10[1] << endl;
   cout << "Fi11 = " <<  Fi11[0] << " and " << Fi11[1] << endl;

   cout << "E1 = " <<  coE1 << endl;
   cout << "E2 = " <<  coE2 << endl;
   cout << "Jk = " <<  coJk << endl;

   for(i = 0; i < 6; i++)
   {
       cout << "I" << i + 1 << " = " << abs(coI[i])*sqrt(2.0) << endl;
	   cout << "argI" << i + 1 << " = " << arg(coI[i])/pi*180.0 << endl;
   }

   cout << "SE1 = " <<  coE1*conj(coI[5]) << endl;
   cout << "SE2 = " <<  coE2*conj(coI[2]) << endl;
   cout << "SJk = " <<  conj(coJk)*(Fi5 - coFi[2]) << endl;

   getch();
}
