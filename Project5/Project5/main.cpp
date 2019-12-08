#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include <tuple>
#include <random>
#include <armadillo>

using namespace std;

//START Foreward Decleration

void ExplicitIntegrator1D(int SpaceN, double DelX, double TotTime);
void ImplicitIntegrator1D(int SpaceN, int TimeN, double DelT, double DelX, double alpha);
void ImplicitCrankNicolson1D(int SpaceN, double DelX, double TotTime);
arma::rowvec TriDiagSolver(int N, double b, double e, arma::rowvec g);
void ExplicitIntegrator2D(int Nx, int Ny, double h, double TotTime);
void ExplicitIntegratorPGP2D(int SpaceN, double LenX, double LenY, double TotTime);

//END Foreward Declaration

int main()
{
    //ExplicitIntegrator1D(100,1,1);
    //ImplicitIntegrator1D(10,10000,0.1,0.005,0.5);
    ExplicitIntegrator2D(10,10,0.01,1);
    return 0;
}

void ExplicitIntegrator1D(int SpaceN, double LenX, double TotTime)
{
    ofstream outfile; outfile.open("Explicit1D.txt");

    double DelX = LenX/SpaceN;
    double DelT = (DelX*DelX)/2;
    double alpha = DelT/(DelX*DelX);
    int TimeN = TotTime/DelT;

    arma::rowvec unow = arma::zeros<arma::rowvec>(SpaceN+1);
    arma::rowvec uprev = arma::zeros<arma::rowvec>(SpaceN+1);
    uprev(0) = unow(0) = 0; uprev(SpaceN) = unow(SpaceN) = 1;

    outfile << 0 << " " << unow;

    for (int j = 1;j <= TimeN ; j++)
    {
        for (int i = 1;i < SpaceN;i++)
        {
            unow(i) = alpha*uprev(i-1)+(1-2*alpha)*uprev(i)+alpha*uprev(i+1);
        }
        outfile << j*DelT << " " << unow;
        uprev = unow;
    }
    outfile.close();
}

void ImplicitIntegrator1D(int SpaceN, int TimeN, double DelT, double DelX, double alpha)
{
    double b = 2*alpha+1; double e = -alpha;
    arma::rowvec unow = arma::zeros<arma::rowvec>(SpaceN+1);
    arma::rowvec uprev = arma::zeros<arma::rowvec>(SpaceN+1); uprev(SpaceN) = 1;

    ofstream outfile; outfile.open("Implicit.txt");
    outfile << 0 << " " << uprev;
    for (int i = 0 ; i < TimeN ; i++)
    {
        unow = TriDiagSolver(SpaceN,b,e,uprev);
        outfile << DelT*(i+1) << " " << unow;
        uprev = unow; uprev(0) = 0; uprev(SpaceN) = 1;
        cout << i << endl;
    }
    outfile.close();
}

void ImplicitCrankNicolson1D(int SpaceN, double DelX, double TotTime)
{
    int a = 0;
}

arma::rowvec TriDiagSolver(int N, double b, double e, arma::rowvec g)
{

    /*Defining the vectors to be used*/
    arma::rowvec btilde = arma::zeros<arma::rowvec>(N+1); // Modified diagonal element
    arma::rowvec gtilde = arma::zeros<arma::rowvec>(N+1); //Modified function value
    arma::rowvec u = arma::zeros<arma::rowvec>(N+2); u(N+1) = 1;

    btilde(0) = b; gtilde(0) = g(0);
    for (int i = 1; i <= N ; i++) //The foreward substitution
    {
        btilde(i) = b - ((e*e)/(double)btilde(i-1));
        gtilde(i) = g(i) - ((gtilde(i-1)*e)/(double)btilde(i-1));
    }


    for (int i = N;i > 0;i--){ //The backwards substitution
        u(i) = ((gtilde(i-1) - e*u(i+1)))/(double)(btilde(i-1));
    }
    return u;
}

void ExplicitIntegrator2D(int Nx, int Ny, double h, double TotTime)
{
    ofstream outfile; outfile.open("Explicit2D.txt");

    double DelT = (h*h)/4;
    double alpha = DelT/(h*h);
    int TimeN = TotTime/DelT;
    cout << TimeN << endl;

    arma::mat unow(Nx+1,Ny+1,arma::fill::zeros);
    arma::mat uprev(Nx+1,Ny+1,arma::fill::zeros);
    uprev.row(0) = arma::ones<arma::rowvec>(Ny+1); unow.row(0) = arma::ones<arma::rowvec>(Ny+1);
    uprev.row(Nx) = arma::ones<arma::rowvec>(Ny+1); unow.row(Nx) = arma::ones<arma::rowvec>(Ny+1);

    outfile << uprev << endl << endl;

    for (int k = 1;k <= TimeN ; k++)
    {
        for (int i = 1;i < Nx;i++)
        {
            for (int j = 1;j < Ny;j++)
            {
                unow(i,j) = uprev(i,j) + alpha*(uprev(i+1,j) + uprev(i-1,j) + uprev(i,j+1) + uprev(i,j-1) - 4*uprev(i,j));
            }
        }
        //outfile << k*DelT << " " << unow << endl;
        outfile << unow << endl << endl;
        uprev = unow;
        cout << k << endl;
    }
    outfile.close();
}
void ExplicitIntegratorPGP2D(int SpaceN, double LenX, double LenY, double TotTime)
{
    int a = 0;
}

