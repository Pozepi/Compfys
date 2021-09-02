#include "functions.hpp"

int main()
{   
    int n;
    double a_value;
    double b_value;
    double c_value;

    std::cin >> n;
    std::cin >> a_value;
    std::cin >> b_value;
    std::cin >> c_value;

    double lin[n];
    double u_zeros[n];
    double a_zeros[n-1];
    double b_zeros[n];
    double c_zeros[n];
    double v_zeros[n];

    double * lins = linspace(0, 1, lin, n);

    double * a_array = fill_array(a_zeros, n, a_value);
    double * b_array = fill_array(b_zeros, n, b_value);
    double * c_array = fill_array(c_zeros, n, c_value);

    double * v_filled = solving_matrix(lins, a_array, b_array, c_array, v_zeros, n);

    
    std::cout << "-----V print-----" << '\n';
    print_array(v_filled, n);

    write_to_file(v_filled, lins, n);

    return 0;
}