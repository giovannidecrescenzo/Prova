#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include "../Random/random.h"
#include "../Functions/functions.h"
#include "../Functions/analysis.h"
#include <vector>


using namespace std;
 
int main (int argc, char *argv[]){

   Random rnd;
   Functions fun;
   Analysis sys;


   ofstream myfile1;
   ofstream myfile2;
   ofstream myfile3;
   ofstream myfile4;
   ofstream myfile5;

   myfile1.open("random_test.txt");
   myfile2.open("sigma.txt");
   myfile3.open("chi.txt");
   myfile4.open("CLT.txt");
   myfile5.open("buffon.txt");

   int seed[4];
   int p1, p2;
   ifstream Primes("../Random/Primes");
   if (Primes.is_open()){
      Primes >> p1 >> p2 ;
   } else cerr << "PROBLEM: Unable to open Primes" << endl;
   Primes.close();

   ifstream input("../Random/seed.in");
   string property;
   if (input.is_open()){
      while ( !input.eof() ){
         input >> property;
         if( property == "RANDOMSEED" ){
            input >> seed[0] >> seed[1] >> seed[2] >> seed[3];
            rnd.SetRandom(seed,p1,p2);
         }
      }
      input.close();
   } else cerr << "PROBLEM: Unable to open seed.in" << endl;



   //Here are defined the values for the exercise
   vector<double> unif;
   vector<double> expo;
   vector<double> lor;
   int M = 1000000; //total number of throws
   int N = 100; //number of blocks
   int K = int(M/N); //number of throws in each block
   double lam = 1.;
   double gamma = 1.;
   double mu = 0.;

   
   for(int i=0; i<N; i++){
      for (int j=0; j<K; j++){
         unif.push_back(rnd.Rannyu());
         lor.push_back(rnd.Lorentz(gamma,mu));
         expo.push_back(rnd.Exp(lam));
      }
   }

   sys = fun.graph(fun.av(unif,N), fun.av2(unif,N));
   fun.print_on_file_errbar_expected(sys.sumprog, sys.errprog, myfile1, 0.5);

   myfile1.close();

   vector<double> integral;

   for (int i = 0; i<N; i++){
      for (int j = 0; j<K; j++){
         integral.push_back(pow((unif[j+i*K]-0.5),2));
      }
   }

   sys = fun.graph(fun.av(integral,N), fun.av2(integral,N));
   fun.print_on_file_errbar_expected(sys.sumprog, sys.errprog, myfile2, 1./12);

   myfile2.close();
   
   // Print on file


   int F = 100; //number of subintervals in which [0,1] is divided
   int N_ = 10000; //how many random numbers I use for each chi square evaluation
   int C = int(M/N_); //number of blocks of N_ numbers
   double exp = N_/F; //expected number of points in each bin
   double count[F]; //counting vector

   for(int i = 0; i<F; i++){
      count[i] = 0;
   }

   for (int k = 0; k<C; k++){
      double sum = 0;
      for(int i = k*N_; i<(k+1)*N_; i++){
         for(int j = 0; j<F; j++){
            if(unif[i]>j*(1./F) and unif[i]<(j+1)*(1./F)){
               count[j] += 1;
            }
         }
      }
      for(int h = 0; h<F; h++){
         sum = sum + (count[h]-exp)*(count[h]-exp);
         count[h] = 0;
      }
      myfile3 << sum/exp << "\t" << 100 << endl;
   }

   myfile3.close();

   int L[4] = {1,2,10,100};

   for (int i = 0; i<4; i++) {
      for (int h = 0; h<10000*L[i]; h+=L[i]){
         double suma = 0.;
         double sumb = 0.;
         double sumc = 0.;
         for (int j = h; j<h+L[i]; j++){
            suma += unif[j];
            sumb += expo[j];
            sumc += lor[j];
         }
         suma = suma/L[i];
         sumb = sumb/L[i];
         sumc = sumc/L[i];

         myfile4 << suma << "\t" << sumb << "\t" << sumc << endl;

      }
   }

   myfile4.close();


   double d = 2.;
   double l = 0.5;

   vector<double> pi;

    for (int i = 0; i<100; i++){
      double count = 0;
      for (int j = 0; j<10000; j++){
         double x = rnd.Rannyu();
         x = x*d;
         double a;
         double b;
         if (x >= d-l) {

            do{
               a = rnd.Rannyu();
               b = rnd.Rannyu();
            }
            while(a*a+b*b > 1.);
            double theta = atan(b/a);
            if (theta < acos((d-x)/l)){
               count +=1;
            }
         }
         if (x <= l){
            do{
               a = rnd.Rannyu();
               b = rnd.Rannyu();
            }
            while(a*a+b*b > 1.);
            double theta = atan(b/a);
            if (theta < acos(x/l)){
               count +=1;
            }
         }
         if (j>100) {
            pi.push_back(4*l*(j+1)/(count*d));
         }
      }
   }


   sys = fun.graph(fun.av(pi,N), fun.av2(pi,N));
   fun.print_on_file_errbar_expected(sys.sumprog, sys.errprog, myfile5, M_PI);
   myfile5.close();


   rnd.SaveSeed();
   
   return 0;
}
