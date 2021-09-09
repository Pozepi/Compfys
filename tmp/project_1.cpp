#include "functions.hpp"

int main()
{   
    // Start measuring time
    clock_t t1 = clock();

   // local variable declaration:
   int problem;
   std::cout << "Chose problem [2/7/8/9]:";
   std::cin >> problem;

   switch(problem) 
   {
      case 2:
      {
         // problem 2
         double start = 0;
         double end   = 1;
         int length;
         std::cin >> length;
         double lin[length]; // master array
         double u_empty[length];

         double x = *linspace(start, end, lin, length);
         double u = *Poisson_analytical(x, u_empty, length);
         std::string filename = "values1";
         write_to_file(filename, u, x, length);
         
         break;
      }

      case 7: 
      {
         // problem 7

         int n;
         double a_value = -1;
         double b_value =  2;
         double c_value = -1;

         std::cout << "length of arrays: ";
         std::cin >> n;
         // std::cin >> a_value;
         // std::cin >> b_value;
         // std::cin >> c_value;

         double lin[n]; // master array

         double u_zeros[n]; // Poisson
         double a_zeros[n-1]; // diagonal a
         double b_zeros[n];   // diagonal b
         double c_zeros[n];   // ...
         double v_zeros[n];   // solution

         double lins = *linspace(0, 1, lin, n); // x array

         double a_array = *fill_array(a_zeros, n, a_value);
         double b_array = *fill_array(b_zeros, n, b_value);
         double c_array = *fill_array(c_zeros, n, c_value);

         double v_filled = *solving_matrix(lins, a_array, b_array, c_array, v_zeros, n);

         
         // std::cout << "-----V print-----" << '\n';
         // print_array(v_filled, n);
         std::string filename = "values2";
         write_to_file(filename, v_filled, lins, n);

         break;
      }

      case 8:
      {
         // problem 8.a., 8.b.

         double start = 0;
         double end   = 1;

         int length;
         std::cout << "length of arrays: ";
         std::cin >> length; // input N length

         double lin[length]; // dummy array to be filled and abused
         double x = *linspace(start, end, lin, length);

         /// ANALYTICAL ///


         double u = *Poisson_analytical(x, lin, length);
         //x = linspace(start, end, lin, length);
         print_array(x, length);

         /// NUMERICAL ///

         double a_value = -1;
         double b_value =  2;
         double c_value = -1;

         double a_zeros[length-1]; // diagonal a
         double b_zeros[length];   // diagonal b
         double c_zeros[length];   // ...

         double a_array = *fill_array(a_zeros, length-1, a_value);
         double b_array = *fill_array(b_zeros, length, b_value);
         double c_array = *fill_array(c_zeros, length, c_value);

         double v = *solving_matrix(x, a_array, b_array, c_array, lin, length);

         /// ERROR ///

         // u is true, v is descrete
         double error = *Delta(u, v, lin, length);
         
         std::string filename0 = "error";
         write_to_file(filename0, error, x, length);

         /// RELATIVE ERROR ///

         double relerror = *epsilon(u, v, lin, length);

         std::string filename1 = "relerror";
         write_to_file(filename1, relerror, x, length);
         
         break;
      }
      /*
      case 9:
      {
         // Problem 8.c.
         int Nn = 7;
         double N[7] = {10, pow(10, 2), pow(10, 3), pow(10, 4), pow(10, 5), pow(10, 6), pow(10, 7)};
         double maxeps[Nn];
         
         for (int i=0;i<Nn;i++)
         {
            double start = 0;
            double end   = 1;
            int length = N[i];
            

            double lin[length]; // dummy array to be filled and abused
            
            double x = *linspace(start, end, lin, length);

            /// ANALYTICAL ///

            double u = *Poisson_analytical(x, lin, length);

            /// NUMERICAL ///

            double a_value = -1;
            double b_value =  2;
            double c_value = -1;

            double a_zeros[length-1]; // diagonal a
            double b_zeros[length];   // diagonal b
            double c_zeros[length];   // ...

            double a_array = *fill_array(a_zeros, length, a_value);
            double b_array = *fill_array(b_zeros, length, b_value);
            double c_array = *fill_array(c_zeros, length, c_value);

            double v = *solving_matrix(x, a_array, b_array, c_array, lin, length);

            /// RELATIVE ERROR ///

            double relerror = *epsilon(u, v, lin, length);
            double max = 0;
            
            for (int j=0;j<length;j++)
            {
               std::cout << relerror[j] << '\n';
               if (max < relerror[j])
               {
                  std::cout << "before";
                  max = relerror[j];
                  std::cout << "after";
               }
            }
            maxeps[i] = max;
         }

         std::string filename = "maxepsofN";
         write_to_file(filename, maxeps, N, Nn);
         
         break;
         
      }*/
   }  


    // Stop measuring time
    clock_t t2 = clock();

    // Calculate the elapsed time.
    double duration_seconds = ((double) (t2 - t1)) / CLOCKS_PER_SEC;

    return 0;
}