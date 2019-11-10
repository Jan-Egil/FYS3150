#include <iostream>
#include <fstream>
#include <armadillo>
#include <omp.h>
#include <time.h>
#include <random>

using namespace std;


//Foreward declaration
tuple<double,double,arma::mat> Init_State_Ordered(int L);
tuple<double,double,arma::mat> Init_State_Random(int L);
tuple<double,double> Calc_Init_Vals(arma::mat A, int L);
//End foreward declaration


int main()
{
    int L = 10; //Size of lattice
    int N = 100000; //Number of cycles
    return 0;
}

tuple<double,double,arma::mat> Init_State_Ordered(int L)
{
    arma::mat A(L,L,arma::fill::ones);
    tuple<double,double> Init_Vals = Calc_Init_Vals(A,L);
    double E_init = get<0>(Init_Vals); double M_init = get<1>(Init_Vals);
    return make_tuple(E_init,M_init,A);
}

tuple<double,double,arma::mat> Init_State_Random(int L)
{
    arma::mat A(L,L); double ranval;
    for (int i = 0;i < L;i++)
    {
        for (int j = 0;j < L;j++)
        {
            int k;
        }
    }
}

tuple<double,double> Calc_Init_Vals(arma::mat A, int L){
    double E_init = 0; double M_init = 0;
    for (int i = 0;i < L;i++)
    {
        for (int j = 0;j < L;j++)
        {
            M_init += A(i,j);
            if (i == 0)
            {
                if (j == 0)
                {
                    E_init += A(i,j)*(A(i,j+1)+A(i+1,j)+A(i,L)+A(L,j));
                }
                else if (j == L)
                {
                    E_init += A(i,j)*(A(i+1,j)+A(L,j)+A(i,j-1)+A(i,0));
                }
                else
                {
                    E_init += A(i,j)*(A(i,j+1)+A(i,j-1)+A(i+1,j)+A(L,j));
                }
            }
            else if (i == L)
            {
                if (j == 0)
                {
                    E_init += A(i,j)*(A(i,j+1)+A(i,L)+A(i-1,j)+A(0,j));
                }
                else if (j == L)
                {
                    E_init += A(i,j)*(A(i,j-1)+A(i,0)+A(i-1,j)+A(0,j));
                }
                else
                {
                    E_init += A(i,j)*(A(i,j+1)+A(i,j-1)+A(i-1,j)+A(0,j));
                }
            }
            else
            {
                if (j == 0)
                {
                    E_init += A(i,j)*(A(i,j+1)+A(i,L)+A(i+1,j)+A(i-1,j));
                }
                else if (j == L)
                {
                    E_init += A(i,j)*(A(i,j-1)+A(i,0)+A(i-1,j)+A(i+1,j));
                }
                else
                {
                    E_init += A(i,j)*(A(i,j+1)+A(i,j-1)+A(i+1,j)+A(i-1,j));
                }
            }
        }
    }
    return make_tuple(E_init,M_init);
}
