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
void flip_n_calc_ordinary(int L, int N, arma::mat Lat, double E, double M, double Temp);
void flip_n_calc_equilibrium(int L, int N, arma::mat Lat, double E, double M, double Temp, string filename);
void flip_n_calc_probdist(int L, int N, arma::mat Lat, double E, double Temp, string filename);
tuple<double,double,double,double> flip_n_calc_ordinary_output(int L, int N, arma::mat Lat, double E, double M, double Temp);
//End foreward declaration




int main()
{
    /*
    The program basically guides you through
    */
    string b;
    cout << "Do you want to do the long-ass simulation that takes forever?" << endl;
    cout << "or do you want to do some more basic simulations?" << endl;
    cout << endl << "'a' for big ass simulation. (Beware, this takes hours), 'b' for anything else" << endl;
    cin >> b;

    if((b == "b") || (b == "B")) //Run a smaller simulation (too be determined)
    {
        int L;
        cout << "How big do you want your (square) lattice to be?" << endl;
        cout << "L = {2,20,40,60,80,100} recommended." << endl;
        cin >> L;

        int N;
        cout << "How many MC sweeps to you want to do?" << endl;
        cin >> N;

        string a;
        cout << "Do you want to create an ordered lattice, or a random lattice?" << endl;
        cout << "Type 'a' for ordered. Type 'b' for random." << endl;
        cin >> a;
        double E; double M; arma::mat Lat;
        if ((a == "a") || (a == "A"))
        {
            tuple<double,double,arma::mat> InitVals = Init_State_Ordered(L);
            E = get<0>(InitVals); M = get<1>(InitVals); Lat = get<2>(InitVals);
        }
        else if ((a == "b") || (a == "B"))
        {
            tuple<double,double,arma::mat> InitVals = Init_State_Random(L);
            E = get<0>(InitVals); M = get<1>(InitVals); Lat = get<2>(InitVals);
        }
        else
        {
            cout << "You didn't type a nor b.. Ya noob.." << endl;
            cout << "Aborting program.." << endl << endl;
            return 0;
        }
        cout << "What temperature do you want to simulate for?" << endl;
        double Temp; cin >> Temp;

        cout << "For simulations where you write to file. What filename do you want to use?" << endl;
        cout << "(By default saved as .txt - input 'testfile' returns 'testfile.txt'" << endl;
        string filename; cin >> filename;

        string sim;
        cout << "Which simulation do you want to run?" << endl << endl;
        cout << "- Press 'a' for a simple simulation returning the final expectation values in terminal" << endl << endl;
        cout << "- Press 'b' for a simulation writing expectation values as a function of time to file" << endl << endl;
        cout << "- Press 'c' for a simulation counting the number of times a state has happened (for probability distribution" << endl << endl;
        cin >> sim;
        if ((sim == "a") || (sim == "A"))
        {
            cout << "Please be aware that this simulation doesn't write anything to file.." << endl;
            cout << "Starting simulation.." << endl;
            flip_n_calc_ordinary(L,N,Lat,E,M,Temp);
            cout << endl << endl << "Simulation finished!" << endl;
        }
        else if((sim == "b") || (sim == "B"))
        {
            cout << "Starting simulation.." << endl;
            flip_n_calc_equilibrium(L,N,Lat,E,M,Temp,filename);
            cout << endl << endl << "Simulation finished!" << endl;
        }
        else if((sim == "c") || (sim == "C"))
        {
            cout << "Starting simulation.." << endl;
            flip_n_calc_probdist(L,N,Lat,E,M,Temp,filename);
            cout << endl << endl << "Simulation finished!" << endl;
        }
        else
        {
            cout << "So close, yet so far.." << endl;
            cout << "You didn't press a, b or c!" << endl;
            cout << "Aborting program.." << endl << endl;
            return 0
        }

    }
    else if ((b == "a") || (b == "A")) //Run big simulation
    {
        int N = 10000000;
        for(int L = 40;L <= 100;L+=20){
            ofstream outfile;

            outfile.open(to_string(L)+"X"+to_string(L)+"_2.txt");
            outfile << "T - <E> - <|M|> - Cv - Chi" << endl;

            tuple<double,double,double,double> tuplevals;
            tuple<double,double,arma::mat> tupleinit;
            double E_start; double M_start; arma::mat Lat;
            double Eret; double Mret; double Cvret; double chiret; double T;

            #pragma omp parallel for default(shared) private(tuplevals,tupleinit,E_start,M_start,Lat,Eret,Mret,Cvret,chiret,T);
            for (T = 2.0;T <= 2.3;T += 0.05)
            {
                tupleinit = Init_State_Ordered(L);
                E_start = get<0>(tupleinit); M_start = get<1>(tupleinit); Lat = get<2>(tupleinit);
                cout << T << " " << L << " Start" << endl;
                tuplevals = flip_n_calc_ordinary_output(L,N,Lat,E_start,M_start,T);
                Eret = get<0>(tuplevals); Mret = get<1>(tuplevals);
                Cvret = get<2>(tuplevals); chiret = get<3>(tuplevals);
                outfile << T << " " << Eret << " " << Mret << " " << Cvret << " " << chiret << endl;
                cout << T << " " << L << " End" << endl;
            }
            outfile.close();
        }
    }
    else
    {
        cout << "You didn't press a or b.." << endl;
        cout << "Learn to follow instructions.." << endl << endl;
        return 0;
    }
    return 0;
}



tuple<double,double,arma::mat> Init_State_Ordered(int L) //Initiates an ordered lattice (all spins same direction)
{
    arma::mat A(L,L,arma::fill::ones); //Creates lattice
    tuple<double,double> Init_Vals = Calc_Init_Vals(A,L); //Calculate Initial values
    double E_init = get<0>(Init_Vals); double M_init = get<1>(Init_Vals); //Extracts and returns energy and magnetization
    return make_tuple(E_init,M_init,A);
}



tuple<double,double,arma::mat> Init_State_Random(int L) //Initiates a random lattice (all spins random directions)
{
    arma::mat A(L,L);

    double randomval;
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0); //Initialize a uniform PDF RNG
    for (int i = 0;i < L;i++) //for-loops giving random spins to all elements
    {
        for (int j = 0;j < L;j++)
        {
            randomval = RandomNumberGenerator(gen);
            if (randomval > 0.5)
            {
                A(i,j) = 1;
            }
            else
            {
                A(i,j) = -1;
            }
        }
    }
    tuple<double,double> Init_Vals = Calc_Init_Vals(A,L); //Calculate initial values of lattice
    double E_init = get<0>(Init_Vals); double M_init = get<1>(Init_Vals); //Extracts and returns energy and magnetization
    return make_tuple(E_init,M_init,A);
}



tuple<double,double> Calc_Init_Vals(arma::mat A, int L) //Calculates the initial energy of a full LXL lattice
{
    double E_init = 0; double M_init = 0; //Initiate energy/magnetization variables
    for (int i = 0;i < L;i++)
    {
        for (int j = 0;j < L;j++)
        {
            M_init += A(i,j);
            if (i == 0) //if-tests checking if boundary-conditions apply
            {
                if (j == 0)
                {
                    E_init -= A(i,j)*(A(i,L-1)+A(L-1,j));
                }
                else
                {
                    E_init -= A(i,j)*(A(i,j-1)+A(L-1,j));
                }
            }
            else
            {
                if (j == 0)
                {
                    E_init -= A(i,j)*(A(i-1,j)+A(i,L-1));
                }
                else
                {
                    E_init -= A(i,j)*(A(i-1,j)+A(i,j-1));
                }
            }
        }
    }
    return make_tuple(E_init,M_init);
}


void flip_n_calc_ordinary(int L, int N, arma::mat Lat, double E, double M, double Temp)
{
    /*
    Calculates <E>, <|M|>, heat capacity and susceptibility of a lattice LXL after N sweeps

    Input:
    L: Lattice Size
    N: Number of MC Sweeps
    Lat: The initial lattice
    E: Initial Energy of lattice
    M Initial Magnetization of lattice
    Temp: Temperature of simulation

    Output:
    Print to terminal with values
    */
    double E_tot = 0; double E2_tot = 0;
    double M_tot = M; double M2_tot = M*M;
    double abs_M_tot = abs(M_tot); double abs_M2_tot = abs(M_tot)*abs(M_tot); //Initiate values

    random_device rand; mt19937_64 gen(rand()); uniform_int_distribution<int> RandintGen(0,L-1); uniform_real_distribution<double> RandFloatGen(0,1); //create RNG

    int Row; int Column; double DeltaE; double DeltaM; double MetroReq; double Abs_M = 0; //initiate more values

    arma::vec DeltaEList(5);
    DeltaEList(0) = exp(8/Temp); DeltaEList(1) = exp(4/Temp);
    DeltaEList(2) = 1;
    DeltaEList(3) = exp(-4/Temp); DeltaEList(4) = exp(-8/Temp); //pre-calculates probabilities
    for (int i = 1;i <= N;i++)
    {
        for (int j = 0;j < L*L;j++)
        {
            Row = RandintGen(gen); Column = RandintGen(gen); //Picks random spin in lattice
            DeltaE = 0; DeltaM = 0; //Just in case resetting values

            /*
            Checks if boundary conditions apply. Calculates DeltaE and DeltaM
            */
            if (Row == 0)
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(L-1,Column)+Lat(Row+1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(L-1,Column)+Lat(Row+1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row+1,Column)+Lat(L-1,Column))*(-2);
                }
            }
            else if (Row == L-1)
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(0,Column)+Lat(Row-1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(0,Column)+Lat(Row-1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row-1,Column)+Lat(0,Column))*(-2);
                }
            }
            else
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(Row+1,Column)+Lat(Row-1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(Row+1,Column)+Lat(Row-1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row-1,Column)+Lat(Row+1,Column))*(-2);
                }
            }
            DeltaM = Lat(Row,Column)*(-2);

            /*
            Metropolis algorithm
            */
            if (DeltaE <= 0)
            {
                Lat(Row,Column) *= (-1); //Swap orientation
            }
            else
            {
                MetroReq = RandFloatGen(gen); //Metropolis requirement
                if (DeltaE == 4)
                {
                    if (MetroReq <= DeltaEList(3))
                    {
                        Lat(Row,Column) *= (-1); //Swap orientation
                    }
                    else
                    {
                        DeltaE = 0; DeltaM = 0;
                    }
                }
                else
                {
                    if (MetroReq <= DeltaEList(4))
                    {
                        Lat(Row,Column) *= (-1);
                    }
                    else
                    {
                        DeltaE = 0; DeltaM = 0;
                    }
                }
            }
            E += DeltaE; M += DeltaM; Abs_M = abs(M);
        }
        E_tot += E; E2_tot += E*E;
        M_tot += M; M2_tot += M*M;
        abs_M_tot += Abs_M; abs_M2_tot += Abs_M*Abs_M; //Add to total
    }

    /*
    Calculate values
    */

    double mean_E = E_tot/(N+1);
    double mean_M = M_tot/(N+1);
    double mean_E2 = E2_tot/(N+1);
    double mean_M2 = M2_tot/(N+1);
    double absmean_M = abs_M_tot/(N+1);
    double absmean_M2 = abs_M2_tot/(N+1);

    double cv = (mean_E2-(mean_E*mean_E))/(Temp*Temp);
    double chi = (mean_M2-(mean_M*mean_M))/Temp;

    /*
    Print to terminal
    */

    cout << "<E> = " << mean_E/(double)(L*L) << endl;
    cout << "<M> = " << mean_M/(double)(L*L) << endl;
    cout << "<|M|> = " << absmean_M/(double)(L*L) << endl;
    cout << "cv = " << cv/(double)(L*L) << endl << "chi = " << chi/(double)(L*L) << endl;
}

void flip_n_calc_equilibrium(int L, int N, arma::mat Lat, double E, double M, double Temp, string filename)
{
    /*
    Calculates <E>, <|M|> and <M> of a lattice LXL after N sweeps. Returns product to file.

    Input:
    L: Lattice Size
    N: Number of MC Sweeps
    Lat: The initial lattice
    E: Initial Energy of lattice
    M Initial Magnetization of lattice
    Temp: Temperature of simulation


    Output:
    Returns values to file for each sweep.
    */

    ofstream outfile;
    outfile.open(filename+".txt"); //initiate file

    double E_tot = E; double M_tot = M; double Abs_M_tot = abs(M); //initiate variables

    outfile << "<E> - <M> - <accepted> - <|M|>" << endl;

    random_device rand; mt19937_64 gen(rand()); uniform_int_distribution<int> RandintGen(0,L-1); uniform_real_distribution<double> RandFloatGen(0,1); //initiate RNG

    int Row; int Column; double DeltaE; double DeltaM; double MetroReq; double Abs_M; //Initiate more variabels

    arma::vec DeltaEList(5); //Pre-calculate probabilities
    DeltaEList(0) = exp(8/Temp); DeltaEList(1) = exp(4/Temp);
    DeltaEList(2) = 1;
    DeltaEList(3) = exp(-4/Temp); DeltaEList(4) = exp(-8/Temp);
    for (int i = 1;i <= N;i++) //Simulation time!
    {
        for (int j = 0;j < L*L;j++)
        {
            Row = RandintGen(gen); Column = RandintGen(gen); //Picks a random point
            DeltaE = 0; DeltaM = 0; //Just in case resets values

            /*
            Checks if boundary conditions apply. Sets DeltaM and DeltaE
            */

            if (Row == 0)
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(L-1,Column)+Lat(Row+1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(L-1,Column)+Lat(Row+1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row+1,Column)+Lat(L-1,Column))*(-2);
                }
            }
            else if (Row == L-1)
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(0,Column)+Lat(Row-1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(0,Column)+Lat(Row-1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row-1,Column)+Lat(0,Column))*(-2);
                }
            }
            else
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(Row+1,Column)+Lat(Row-1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(Row+1,Column)+Lat(Row-1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row-1,Column)+Lat(Row+1,Column))*(-2);
                }
            }
            DeltaM = Lat(Row,Column)*(-2);

            /*
            Metropolis Algorithm
            */
            if (DeltaE <= 0)
            {
                Lat(Row,Column) *= (-1); //Swap orientation
            }
            else
            {
                MetroReq = RandFloatGen(gen); //Maybe accept
                if (DeltaE == 4)
                {
                    if (MetroReq <= DeltaEList(3))
                    {
                        Lat(Row,Column) *= (-1);
                    }
                    else
                    {
                        DeltaE = 0; DeltaM = 0;
                    }
                }
                else
                {
                    if (MetroReq <= DeltaEList(4))
                    {
                        Lat(Row,Column) *= (-1);
                    }
                    else
                    {
                        DeltaE = 0; DeltaM = 0;
                    }
                }
            }
            E += DeltaE; M += DeltaM; Abs_M = abs(M);
        }
        E_tot += E;
        M_tot += M; Abs_M_tot += Abs_M; //Update values



        outfile << E_tot/(L*L*(i+1)) << " " << M_tot/(L*L*(i+1)) << " " << Abs_M_tot/(L*L*(i+1)) << endl; //Returns each sweep result to file

    }
    outfile.close();
}

void flip_n_calc_probdist(int L, int N, arma::mat Lat, double E, double Temp, string filename)
{
    /*
    Calculates <E^2> and E. Returns product to file.

    Input:
    L: Lattice Size
    N: Number of MC Sweeps
    Lat: The initial lattice
    E: Initial Energy of lattice
    M Initial Magnetization of lattice
    Temp: Temperature of simulation
    Filename: Name of file you want to return values to (without .txt at the end)


    Output:
    Returns values to file for each sweep.
    */

    ofstream outfile; //Start file
    outfile.open(filename+".txt");

    double E_tot = E; double E2_tot = E*E; //Initiate variables

    random_device rand; mt19937_64 gen(rand()); uniform_int_distribution<int> RandintGen(0,L-1); uniform_real_distribution<double> RandFloatGen(0,1); //initiate RNG

    int Row; int Column; double DeltaE; double MetroReq; double Abs_M; double E_mean; // Initiate more variables

    arma::vec DeltaEList(5);
    DeltaEList(0) = exp(8/Temp); DeltaEList(1) = exp(4/Temp);
    DeltaEList(2) = 1;
    DeltaEList(3) = exp(-4/Temp); DeltaEList(4) = exp(-8/Temp);
    for (int i = 1;i <= N;i++)
    {
        for (int j = 0;j < L*L;j++)
        {
            Row = RandintGen(gen); Column = RandintGen(gen); //Pick random spin in lattice
            DeltaE = 0; //Juuuust in case

            /*
            Check if boundary conditions are met. Update DeltaE
            */

            if (Row == 0)
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(L-1,Column)+Lat(Row+1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(L-1,Column)+Lat(Row+1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row+1,Column)+Lat(L-1,Column))*(-2);
                }
            }
            else if (Row == L-1)
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(0,Column)+Lat(Row-1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(0,Column)+Lat(Row-1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row-1,Column)+Lat(0,Column))*(-2);
                }
            }
            else
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(Row+1,Column)+Lat(Row-1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(Row+1,Column)+Lat(Row-1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row-1,Column)+Lat(Row+1,Column))*(-2);
                }
            }

            /*
            Metropolis Algorithm, go!
            */

            if (DeltaE <= 0)
            {
                Lat(Row,Column) *= (-1); //Swap orientation
            }
            else
            {
                MetroReq = RandFloatGen(gen); //Maybe accept
                if (DeltaE == 4)
                {
                    if (MetroReq <= DeltaEList(3))
                    {
                        Lat(Row,Column) *= (-1);
                    }
                    else
                    {
                        DeltaE = 0;
                    }
                }
                else
                {
                    if (MetroReq <= DeltaEList(4))
                    {
                        Lat(Row,Column) *= (-1);
                    }
                    else
                    {
                        DeltaE = 0;
                    }
                }
            }
            E += DeltaE; //New Energy of lattice
        }
        E_tot += E; E2_tot += E*E;
        if (i > N/(double)10) //Only start returning if steady state is reached
        {
            outfile << E << endl; //Okay boomer, now you can return.
        }

    }
    outfile << E2_tot/(N+1) << endl;
    outfile.close();
}

tuple<double,double,double,double> flip_n_calc_ordinary_output(int L, int N, arma::mat Lat, double E, double M, double Temp)
{
    /*
    Calculates <E>, <|M|>, heat capacity and susceptibility of a lattice (includes a "stabilization period")

    Input:
    L: Lattice Size
    N: Number of MC Sweeps
    Lat: The initial lattice
    E: Initial Energy of lattice
    M Initial Magnetization of lattice
    Temp: Temperature of simulation

    Output:
    tuple of values:
    [0]: Mean energy of lattice
    [1]: Mean (absolute) magnetization of lattice
    [2]: Heat capacity of lattice
    [3]: Susceptibility of lattice
    */

    double E_tot = E; double E2_tot = E*E;
    double M_tot = M; double M2_tot = M*M;
    double abs_M_tot = abs(M_tot); double abs_M2_tot = abs(M_tot)*abs(M_tot); //Initiate variables

    random_device rand; mt19937_64 gen(rand()); uniform_int_distribution<int> RandintGen(0,L-1); uniform_real_distribution<double> RandFloatGen(0,1); //Initiate RNG

    int Row; int Column; double DeltaE; double DeltaM; double MetroReq; double Abs_M = 0;

    arma::vec DeltaEList(5); //Pre-calculate probabilities (for Metropolis algo)
    DeltaEList(0) = exp(8/Temp); DeltaEList(1) = exp(4/Temp);
    DeltaEList(2) = 1;
    DeltaEList(3) = exp(-4/Temp); DeltaEList(4) = exp(-8/Temp);
    for (int i = 1;i <= N;i++)
    {
        for (int j = 0;j < L*L;j++)
        {
            Row = RandintGen(gen); Column = RandintGen(gen); //Find a random row and column for spin flip
            DeltaE = 0; DeltaM = 0; //Resets values (just in case)

            /*
            Set DeltaE and DeltaM.
            */

            if (Row == 0) //if-tests to check if we need to apply periodic boundary conditions
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(L-1,Column)+Lat(Row+1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(L-1,Column)+Lat(Row+1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row+1,Column)+Lat(L-1,Column))*(-2);
                }
            }
            else if (Row == L-1)
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(0,Column)+Lat(Row-1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(0,Column)+Lat(Row-1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row-1,Column)+Lat(0,Column))*(-2);
                }
            }
            else
            {
                if(Column == 0)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,L-1)+Lat(Row,Column+1)+Lat(Row+1,Column)+Lat(Row-1,Column))*(-2);
                }
                else if(Column == L-1)
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,0)+Lat(Row+1,Column)+Lat(Row-1,Column))*(-2);
                }
                else
                {
                    DeltaE = -Lat(Row,Column)*(Lat(Row,Column-1)+Lat(Row,Column+1)+Lat(Row-1,Column)+Lat(Row+1,Column))*(-2);
                }
            }
            DeltaM = Lat(Row,Column)*(-2);

            /*
            Start Metropolis Test / Algorithm
            */

            if (DeltaE <= 0)
            {
                Lat(Row,Column) *= (-1); //Swap orientation
            }
            else
            {
                MetroReq = RandFloatGen(gen); //Maybe accept
                if (DeltaE == 4)
                {
                    if (MetroReq <= DeltaEList(3)) //
                    {
                        Lat(Row,Column) *= (-1);
                    }
                    else
                    {
                        DeltaE = 0; DeltaM = 0;
                    }
                }
                else
                {
                    if (MetroReq <= DeltaEList(4))
                    {
                        Lat(Row,Column) *= (-1);
                    }
                    else
                    {
                        DeltaE = 0; DeltaM = 0;
                    }
                }
            }
            E += DeltaE; M += DeltaM; Abs_M = abs(M);
        }
        if(i > N/(double)10) //Check if we're still waiting for steady state
        {
            E_tot += E; E2_tot += E*E; //If yes, let's update values!
            M_tot += M; M2_tot += M*M;
            abs_M_tot += Abs_M; abs_M2_tot += Abs_M*Abs_M;
        }
    }

    /*
    Calculate means and return values
    */

    double mean_E = E_tot/(N+1-(N/(double)10));
    double mean_M = M_tot/(N+1-(N/(double)10));
    double mean_E2 = E2_tot/(N+1-(N/(double)10));
    double mean_M2 = M2_tot/(N+1-(N/(double)10));
    double absmean_M = abs_M_tot/(N+1-(N/(double)10));
    double absmean_M2 = abs_M2_tot/(N+1-(N/(double)10));

    double cv = (mean_E2-(mean_E*mean_E))/(Temp*Temp);
    double chi = (mean_M2-(mean_M*mean_M))/Temp;

    return make_tuple(mean_E/(double)(L*L),absmean_M/(double)(L*L),cv/(double)(L*L),chi/(double)(L*L));
}
