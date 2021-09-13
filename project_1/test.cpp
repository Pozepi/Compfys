#include "functions.hpp"
int main()
{   
   // local variable declaration:
   int problem;
   std::cout << "Choose problem: \n";
   std::cout << "[1] Problem 2\n";
   std::cout << "[2] Problem 7\n";         // problem 2
   std::cout << "[3] Problem 8\n";
   std::cout << "[4] Problem 8c\n";
   std::cout << "[5] Problem 10\n";

   std::cin >> problem;
   switch(problem)
   {
      case 1:
      {
         // problem 2
         double start = 0;
         double end   = 1;
         int length;
         std::cout << "Size of array: \n";
         std::cin >> length;

         double x[length]; // define x and u
         double u[length]; 

         linspace(x, start, end, length); // fill x with linspace values
         Poisson_analytical(u, x, length);

         std::string filename = "values1";
         write_to_file(filename, u, x, length);
         
         break;
      }

      case 2: 
      {
         // problem 7

         int N[4] = {10, 100, 1000, 10000};
         for (int i=0; i<4; i++)
         {
            int n = N[i];

            double a_value = -1;
            double b_value =  2;
            double c_value = -1;

            // std::cin >> a_value;
            // std::cin >> b_value;
            // std::cin >> c_value;

            double x[n]; // empty x array
            double u[n]; // empty Poisson array

            double a_array[n-1]; // diagonal a
            double b_array[n];   // diagonal b
            double c_array[n];   // ...

            double v[n];   // solution

            linspace(x, 0, 1, n); // x array

            fill_array(a_array, n, a_value);
            fill_array(b_array, n, b_value);
            fill_array(c_array, n, c_value);

            solving_matrix(v, x, a_array, b_array, c_array, n);

            
            // std::cout << "-----V print-----" << '\n';
            // print_array(v_filled, n);
            // std::string filename = "values2";

            std::string filename = "u_values_N_";
            std::string s = std::to_string(i);
            write_to_file(filename+s, v, x, n);
         }

         break;
      }

      case 3:
      {
         // problem 8.a., 8.b.

         double start = 0;
         double end   = 1;
         
         int lengths[4] = {10, 100, 1000, 10000};
         for (int i=0;i<4;i++)
         {
            int length = lengths[i];

            double x[length]; // dummy x-array to be filled and abused

            linspace(x, start, end, length);

            /// ANALYTICAL ///
            double u[length];

            Poisson_analytical(u, x, length);

            /// NUMERICAL ///
            double v[length];

            double a_value = -1;
            double b_value =  2;
            double c_value = -1;

            double a_array[length-1]; // diagonal a
            double b_array[length];   // diagonal b
            double c_array[length];   // ...

            fill_array(a_array, length-1, a_value);
            fill_array(b_array, length, b_value);
            fill_array(c_array, length, c_value);

            solving_matrix(v, x, a_array, b_array, c_array, length);

            /// ERROR ///

            // u is true, v is descrete
            double error[length];
            
            Delta(error, u, v, length);
            

            /// RELATIVE ERROR ///
            double relerror[length];
            epsilon(relerror, u, v, length);

            /// chop chop, it seems the arrays are too big. Lets cut them, shall we ;)
            double x_[length-2];
            double relerror_[length-2];
            double error_[length-2];
            
            int j = 0;
            for (int k=1;k<length-1;k++)
            {
               x_[j] = x[k];
               relerror_[j] = relerror[k];
               error_[j] = error[k];
               j++;
            }

            std::string filename0 = "error_values_N_";
            std::string s = std::to_string(i);

            write_to_file(filename0+s, error_, x_, length-2);
            std::string filename1 = "rel_error_values_N_";
            write_to_file(filename1+s, relerror_, x_, length-2);
         }
         
         break;
      }
      
      case 4:
      {
         // Problem 8.c.
         double start = 0;
         double end   = 1;

         double N[5] = {10, 100, 1000, 10000, 100000};
         double maxeps[5];
         
         
         for (int i=0;i<5;i++)
         {
            int length = N[i];
            
            double x[length]; // dummy x array to be filled and abused
            
            linspace(x, start, end, length);
            
            /// ANALYTICAL ///
            double u[length];
            Poisson_analytical(u, x, length);

            /// NUMERICAL ///
            double v[length];

            solving_special(v, x, length);
            // print_array(v, length);

            /// RELATIVE ERROR ///
            
            double relerror[length];
            
            epsilon(relerror, u, v, length);

            // Cut down array (remove boundaries)
            double x_[length-2];
            double relerror_[length-2];
            
            int j = 0;
            for (int k=1;k<length-1;k++)
            {
               x_[j] = x[k];
               relerror_[j] = relerror[k];
               j++;
            }
            // Find max
            double max = 0;
            for (int l = 0; l < length-2; l++)
            {
               if (relerror_[l] > max)
               {
                  max = relerror_[l];
               }
            }
            maxeps[i] = max; 
            //maxeps[i] = *std::max_element(relerror, relerror+length);
            std::cout << maxeps[i] << '\n';
            
         }
         std::string filename = "max_eps";

         //double Nin[] = {10, pow(10,2), pow(10,3), pow(10,4), pow(10,5), pow(10,6), pow(10,7)};
         write_to_file(filename, maxeps, N, 5);

         break;
      }

      case 5:
      {
         // Problem 10
         double start = 0;
         double end   = 1;
         double N[5] = {10, 100, 1000, 10000, 100000};

         // General Algorithm

         double time_in[50];
         double N_in[50];
         
         int k = 0;
         for (int i=0;i<10;i++)
         {
            for (int j=0;j<5;j++)
            {
               // Start measuring time
               clock_t t1 = clock();

               // CODE : 
               int n = N[j];

               double a_value = -1;
               double b_value =  2;
               double c_value = -1;

               double x[n]; // empty x array
               double u[n]; // empty Poisson array

               double a_array[n-1]; // diagonal a
               double b_array[n];   // diagonal b
               double c_array[n];   // ...

               double v[n];   // solution

               linspace(x, 0, 1, n); // x array

               fill_array(a_array, n, a_value);
               fill_array(b_array, n, b_value);
               fill_array(c_array, n, c_value);

               solving_matrix(v, x, a_array, b_array, c_array, n);

               // Stop measuring time
               clock_t t2 = clock();

               // Calculate the elapsed time.
               double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

               std::cout << "\nRuntime (s) : " << duration_seconds;

               time_in[k] = duration_seconds;
               N_in[k]    = n; 
               k++;
            }
         }

         write_to_file("time_general", time_in, N_in, 50);

         // Special Algorithm


         double time_in2[50];
         double N_in2[50];
         
         int k2 = 0;
         for (int i=0;i<10;i++)
         {
            for (int j=0;j<5;j++)
            {
               // Start measuring time
               clock_t t1 = clock();

               // CODE : 
               int length = N[j];
               
               double x[length]; // dummy x array to be filled and abused
               
               linspace(x, start, end, length);

               /// NUMERICAL ///
               double v[length];

               solving_special(v, x, length);

               // Stop measuring time
               clock_t t2 = clock();

               // Calculate the elapsed time.
               double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

               std::cout << "\nRuntime (s) : " << duration_seconds;

               time_in2[k2] = duration_seconds;
               N_in2[k2]    = length; 
               k2++;
            }
         }

         write_to_file("time_special", time_in2, N_in2, 50);

      }
      default:
         std::cout << "Default statement";
         break;
   }  

   return 0;
}