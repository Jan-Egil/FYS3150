#include <iostream>
#include <cmath>
#include <fstream>
#include <time.h>
#include <armadillo>

using namespace std;

//START Foreward Decleration

void ExplicitIntegrator1D(int SpaceN, double LenX, double TotTime);
void ImplicitIntegrator1D(int SpaceN, double LenX, double TotTime);
void ImplicitCrankNicolson1D(int SpaceN, double LenX, double TotTime);
arma::rowvec TriDiagSolver(int N, double b, double e, arma::rowvec g);
void ExplicitIntegrator2D(int Nx, int Ny, double h, double TotTime);

double QradHalfLife(double x, double y, double t);
void ExplicitIntegratorPGP2D(int Nx, int Ny, double h, double TotTime);

//END Foreward Declaration

int main()
{
    clock_t start = clock(); //START clock

    /*
    Remove comment in program you want to run. Set your own parameters
    */

    //ExplicitIntegrator1D(10,1,1);
    //ImplicitIntegrator1D(10,1,1);
    //ImplicitCrankNicolson1D(10,1,1);
    //ExplicitIntegrator2D(10,10,0.01,1);
    //ExplicitIntegratorPGP2D(351,121,1,1);

    clock_t stop = clock(); //STOP clock
    double time = (double)(stop-start)/CLOCKS_PER_SEC; //Calculate time
    cout << "Time spent: " << time*1000 << " ms" << endl;

    return 0;
}

void ExplicitIntegrator1D(int SpaceN, double LenX, double TotTime)
{
    /*
    Takes:
    SpaceN - Number of (spatial) integration points
    LenX - Length of "simulation box".
    TotTime - Total time to simulate
    */

    ofstream outfile; outfile.open("Explicit1D.txt"); //Open output file

    /*
    Initiate variables and arrays, based on Explicit Algorithm limitations
    */
    double DelX = LenX/SpaceN;
    double DelT = (DelX*DelX)/2;
    double alpha = DelT/(DelX*DelX);
    int TimeN = TotTime/DelT;

    arma::rowvec unow = arma::zeros<arma::rowvec>(SpaceN+1);
    arma::rowvec uprev = arma::zeros<arma::rowvec>(SpaceN+1);
    uprev(0) = unow(0) = 0; uprev(SpaceN) = unow(SpaceN) = 1;

    arma::rowvec actual = arma::zeros<arma::rowvec>(SpaceN+1); //Analytic value array

    for (int k = 0;k <= SpaceN ;k++)//Set analytic values (for calculating error)
    {
        actual(k) = k/(double)(SpaceN/LenX);
    }
    double error;
    /*
    Start "simulation"
    */
    outfile << 0 << " " << 1 << " " << unow;

    for (int j = 1;j <= TimeN ; j++) //Time-iteration
    {
        for (int i = 1;i < SpaceN;i++) //Space-iteration
        {
            unow(i) = alpha*uprev(i-1)+(1-2*alpha)*uprev(i)+alpha*uprev(i+1); //Calculate new values
        }
        error = sqrt(arma::sum(arma::dot(unow-actual,unow-actual)))/sqrt(arma::sum(arma::dot(actual,actual)));

        outfile << j*DelT << " " << error << " " << unow; //Output to file.

        uprev = unow;
    }
    outfile.close();
}

void ImplicitIntegrator1D(int SpaceN, double LenX, double TotTime)
{
    /*
    Takes:
    SpaceN - Number of (spatial) integration points
    LenX - Length of "simulation box".
    TotTime - Total time to simulate
    */

    /*
    Initiate variables and arrays, based on Explicit Algorithm limitations
    */
    double DelX = LenX/(double)SpaceN;
    double DelT = (DelX*DelX)/(double)2;
    double alpha = DelT/(double)(DelX*DelX);
    int TimeN = TotTime/(double)DelT;
    double error;

    double b = (2*alpha)+1; double e = -alpha;
    arma::rowvec unow = arma::zeros<arma::rowvec>(SpaceN+1); unow(SpaceN) = 1;
    arma::rowvec uprev = arma::zeros<arma::rowvec>(SpaceN+1); uprev(SpaceN) = 1;

    ofstream outfile; outfile.open("Implicit.txt"); //Open output file

    arma::rowvec actual = arma::zeros<arma::rowvec>(SpaceN+1); //Analytic value array

    for (int k = 0;k <= SpaceN ;k++)//Set analytic values (for calculating error)
    {
        actual(k) = k/(double)(SpaceN/LenX);
    }

    outfile << 0 << " " << 1 << " " << unow << endl;

    /*
    Start simulation
    */
    for (int i = 0 ; i < TimeN ; i++)
    {
        unow = TriDiagSolver(SpaceN,b,e,uprev);

        error = sqrt(arma::sum(arma::dot(unow-actual,unow-actual)))/sqrt(arma::sum(arma::dot(actual,actual)));

        outfile << i*DelT << " " << error << " " << unow;

        uprev = unow;
    }
    outfile.close();
}

void ImplicitCrankNicolson1D(int SpaceN, double LenX, double TotTime)
{
    /*
    Takes:
    SpaceN - Number of (spatial) integration points
    LenX - Length of "simulation box".
    TotTime - Total time to simulate
    */

    /*
    Initiate variables and arrays, based on Explicit Algorithm limitations
    */
    double DelX = LenX/(double)SpaceN;
    double DelT = (DelX*DelX)/(double)2;
    double alpha = DelT/(double)(DelX*DelX);
    int TimeN = TotTime/(double)DelT;

    double b = (2*alpha)+2; double e = -alpha;
    arma::rowvec unow = arma::zeros<arma::rowvec>(SpaceN+1); unow(SpaceN) = 1;
    arma::rowvec uprev = arma::zeros<arma::rowvec>(SpaceN+1); uprev(SpaceN) = 1;
    arma::rowvec r = arma::zeros<arma::rowvec>(SpaceN+1); r(SpaceN) = 1;

    double error;

    arma::rowvec actual = arma::zeros<arma::rowvec>(SpaceN+1); //Analytic-value array

    ofstream outfile; outfile.open("CrankNicolson.txt"); //Output file

    outfile << 0 << " " << 1 << " " << unow; //Initial state

    for (int k = 0;k <= SpaceN ;k++)//Set analytic values (for calculating error)
    {
        actual(k) = k/(double)(SpaceN/LenX);
    }

    for(int i = 0; i <= TimeN ; i++)
    {
        for(int j = 1;j < SpaceN;j++)
        {
            r(j) = alpha*uprev(j-1)+(2-2*alpha)*uprev(j)+alpha*uprev(j+1);
        }
        unow = TriDiagSolver(SpaceN,b,e,r);

        error = sqrt(arma::sum(arma::dot(unow-actual,unow-actual)))/sqrt(arma::sum(arma::dot(actual,actual))); //Finite Difference Scheme

        outfile << i*DelT << " " << error << " " << unow;

        uprev = unow;
    }
    outfile.close();
}

arma::rowvec TriDiagSolver(int N, double b, double e, arma::rowvec g) //Tridiagonal solver for Ax = b-system.
{

    /*Defining the vectors to be used*/
    arma::rowvec btilde = arma::zeros<arma::rowvec>(N+1); // Modified diagonal element
    arma::rowvec gtilde = arma::zeros<arma::rowvec>(N+1); //Modified function value
    arma::rowvec u = arma::zeros<arma::rowvec>(N+1); u(0) = g(0); u(N) = g(N);

    btilde(0) = b; gtilde(0) = g(0);
    for (int i = 1; i <= N ; i++) //The foreward substitution
    {
        btilde(i) = b - ((e*e)/(double)btilde(i-1));
        gtilde(i) = g(i) - ((gtilde(i-1)*e)/(double)btilde(i-1));
    }


    for (int i = (N-1);i > 0;i--){ //The backwards substitution
        u(i) = ((gtilde(i) - e*u(i+1)))/(double)(btilde(i));
    }
    return u;
}

void ExplicitIntegrator2D(int Nx, int Ny, double h, double TotTime)
{
    /*
    Takes:
    Nx - Number of (spatial) integration points in the x-direction
    Ny - Number of (spatial) integration points in the y-direction
    h - Spatial step size
    TotTime - Total time to simulate
    */

    ofstream outfile; outfile.open("Explicit2D.txt"); //Outfile for results

    /*
    Initiate variables and arrays, based on Explicit Algorithm (2D) limitations
    */
    double DelT = (h*h)/4;
    double alpha = DelT/(h*h);
    int TimeN = TotTime/DelT;

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
        outfile << unow << endl << endl;
        uprev = unow;
        cout << k << endl;
    }
    outfile.close();
}

double QradHalfLife(double x, double y, double t) //Calculates the Heat Production Q (taking Half Life into account)
{
    double returnval;
    if (y < 20)
    {
        returnval = 1.4;
    }
    else if (y < 40)
    {
        returnval = 0.35;
    }
    else
    {
        if (x > 100 && x < 250)
        {
            returnval = 0.05+0.5*(0.4*pow(0.5,t/4.47)+0.4*pow(0.5,t/14.0)+0.2*pow(0.5,t/1.25));
        }
        else
        {
            returnval = 0.05;
        }
    }
    return returnval;
}

void ExplicitIntegratorPGP2D(int Nx, int Ny, double h, double TotTime)
{
    /*
    Takes:
    Nx - Number of (spatial) integration points in the x-direction
    Ny - Number of (spatial) integration points in the y-direction
    h - Spatial step size
    TotTime - Total time to simulate
    */

    ofstream outfile; outfile.open("PGPExplicit2D.txt"); //Outfile for results

    /*
    Initiate variables and arrays, based on Explicit Algorithm (2D w/ units) limitations
    */
    double DelT = 1e-5;
    double k = 2500; double rho = 3.5e12; double cp = 1000;
    int TimeN = TotTime/DelT;
    double alpha = (k*DelT)/(rho*cp*h*h);

    arma::mat Tnow(Nx,Ny,arma::fill::zeros);
    arma::mat Tprev(Nx,Ny,arma::fill::zeros);
    Tprev.row(0) = arma::ones<arma::rowvec>(Ny)*8; Tnow.row(0) = arma::ones<arma::rowvec>(Ny)*8;
    Tprev.row(Nx-1) = arma::ones<arma::rowvec>(Ny)*1300; Tnow.row(Nx-1) = arma::ones<arma::rowvec>(Ny)*1300;

    double x; double y; double t;
    /*
    Start Simulation
    */
    for (int k = 1; k <= TimeN ; k++)
    {
        t = k*DelT;
        for (int i = 0; i <= (Nx-1); i++)
        {
            x = i*h;
            for (int j = 1; j < (Ny-1); j++)
            {
                y = j*h;
                if (i == 0) //Different cases, to account for periodic boundary conditions
                {
                    Tnow(i,j) = Tprev(i,j) + DelT*QradHalfLife(x,y,t)/(rho*cp) + alpha*(Tprev(i+1,j)+Tprev(Nx-1,j)+Tprev(i,j+1)+Tprev(i,j-1)-4*Tprev(i,j));
                }
                else if (i == (Nx-1))
                {
                    Tnow(i,j) = Tprev(i,j) + DelT*QradHalfLife(x,y,t)/(rho*cp) + alpha*(Tprev(0,j)+Tprev(i-1,j)+Tprev(i,j+1)+Tprev(i,j-1)-4*Tprev(i,j));
                }
                else
                {
                    Tnow(i,j) = Tprev(i,j) + DelT*QradHalfLife(x,y,t)/(rho*cp) + alpha*(Tprev(i+1,j)+Tprev(i-1,j)+Tprev(i,j+1)+Tprev(i,j-1)-4*Tprev(i,j));
                }
            }
        }
        cout << k << endl;
        Tprev = Tnow;
        if (k%10==0) //Only output every 10 iteration. (To avoid too large files)
        {
            outfile << Tnow << endl;
        }

    }

    outfile.close();
}

