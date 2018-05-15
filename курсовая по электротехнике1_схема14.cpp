#include <complex>
#include <iostream>
#include <conio.h>


int main( )
{
   using namespace std;
   double pi = 3.14159265359;
   double E1m = 100, E2m = 135, E1 = E1m/sqrt(2.0), E2 = E2m/sqrt(2.0);
   double Jkm = 18, Jk = Jkm/sqrt(2.0);
   double psi[3] = {-90.0/(2*pi), 160.0/(2*pi), -30.0/(2*pi)};
   double f = 50;
   double R[4] = {7, 3, 8, 10};
   double L[3] = {0.030, 0.150, 0.060};
   double C[3] = {0.000600, 0.000100, 0.000300};

   double xC[3], xL[3];
 
   complex <double> one ( 1.0 , 0.0 );
   complex <double> min_one ( -1.0 , 0.0 );
   complex <double> j ( 0.0 , 1.0 );
   complex <double> Z[6];
   complex <double> coFi[4];
   complex <double> Y11, Y12, Y13, Y32, Y33;
   complex <double> Iuz1, Iuz3;
   complex <double> coE1, coE2, coJk;
   complex <double> delta, deltaFi1, deltaFi3;
   complex <double> coI[6];
   complex <double> Sist ( 0.0 , 0.0 ), Spr ( 0.0, 0.0 );
   
   for(int i = 0; i < 3; i++)
   {
	   xC[i] = 1/(2*pi*f*C[i]);
	   xL[i] = 2*pi*f*L[i];
	   cout << "xC" << i+1 << " = " << xC[i] << endl;
	   cout << "xL" << i+1 << " = " << xL[i] << endl;
   }

   coE1 = complex <double> (polar(sqrt(E1), psi[0]));
   coE2 = complex <double> (polar(sqrt(E2), psi[1]));
   coJk = complex <double> (polar(sqrt(Jk), psi[2]));
   
   coFi[3] = complex <double> (0.0, 0.0);
   coFi[1] = min_one * coE2;
   
   Z[0] = complex <double> (R[0], xL[0]);
   Z[1] = complex <double> (9999999909999999.0, 9999990999999999.0);
   Z[2] = complex <double> (R[2], -xC[2]);
   Z[3] = complex <double> (R[1], -xC[0]);
   Z[4] = complex <double> (0.0, xL[1] - xC[1]);
   Z[5] = complex <double> (0.0, 0.0);

   Y11 = one/Z[0] + one/Z[4] + one/Z[1];
   Y33 = one/Z[0] + one/Z[3] + one/Z[2];
   Y12 = min_one/Z[1];
   Y13 = min_one/Z[0];
   Y32 = min_one/Z[2];

   cout << "Y11 = " << Y11 << endl;
   cout << "Y33 = " << Y33 << endl;
   cout << "Y12 = " << Y12 << endl;
   cout << "Y13 = " << Y13 << endl;
   cout << "Y32 = " << Y32 << endl;

   Iuz1 = coJk;
   Iuz3 = min_one*coE1/Z[2];

   cout << "Iuz1 = " << Iuz1 << endl;
   cout << "Iuz3 = " << Iuz3 << endl;

   delta = Y11*Y33 - Y13*Y13;
   deltaFi1 = (Iuz1 - coFi[1]*Y12)*Y33 - (Iuz3 - coFi[1]*Y32)*Y13;
   deltaFi3 = Y11*(Iuz3 - coFi[1]*Y32) - Y13*(Iuz1 - coFi[1]*Y12);

   coFi[0] = deltaFi1/delta;
   coFi[2] = deltaFi3/delta;

   coI[0] = (coFi[2] - coFi[0])/Z[0];
   coI[1] = -coJk;
   coI[2] = (coFi[1] - coFi[2] - coE1)/Z[2];
   coI[3] = (coFi[2] - coFi[3])/Z[3];
   coI[4] = (coFi[0] - coFi[3])/Z[4];
   //coI[5] = (coFi[1] - coFi[3])/Z[5];

   coI[5] = coI[1] - coI[2];
   
   cout << "Fi1 = " << coFi[0] << endl;
   cout << "Fi2 = " << coFi[1] << endl;
   cout << "Fi3 = " << coFi[2] << endl;
   cout << "Fi4 = " << coFi[3] << endl;

   for(i = 0; i < 6; i++)
      cout << "I" << i+1 << " = " << coI[i] << endl;
     
   Sist = -coE1*conj(coI[2]) + coE2*conj(coI[5]) +
	   conj(coJk)*(coFi[0] + coJk*j*xL[2] - (coFi[1]-coJk*R[3]));
   
   for(i = 0; i < 6; i++)
   {
	   if(i == 1) i++;
	   Spr += abs(coI[i])*abs(coI[i])*Z[i];
   }

   Spr += abs(coJk) * abs(coJk) * (R[3] + j*xL[2]);
   
   cout << "Sist = " << Sist << endl;
   cout << "Spr  = " << Spr << endl;
   
   cout << "I1 - I2 - I5  = " << coI[0] - coI[1] - coI[4] << endl;
   cout << "I3 - I1 - I4  = " << coI[2] - coI[0] - coI[3] << endl;
   
   getch();
}
