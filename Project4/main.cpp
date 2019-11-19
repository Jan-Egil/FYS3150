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
void flip_n_calc_probdist(int L, int N, arma::mat Lat, double E, double M, double Temp, string filename);
tuple<double,double,double,double> flip_n_calc_ordinary_output(int L, int N, arma::mat Lat, double E, double M, double Temp);
//End foreward declaration


int main()
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

    double Temp = 1;
    string filename = "histogramdata1";
    //flip_n_calc_equilibrium(L,N,Lat,E,M,Temp,filename);
    //flip_n_calc_ordinary(L,N,Lat,E,M,Temp);
    //flip_n_calc_probdist(L,N,Lat,E,M,Temp,filename);
    for(int L = 40;L <= 100;L+=20){
        ofstream outfile;
        outfile.open(to_string(L)+"X"+to_string(L)+".txt");
        outfile << "T - <E> - <|M|> - Cv - Chi" << endl;
        tuple<double,double,double,double> tuplevals;
        double Eret; double Mret; double Cvret; double chiret; double T;

        for (double T = 2.0;T <= 2.3;T += 0.05)
        {
            cout << T << endl;
            tuplevals = flip_n_calc_ordinary_output(L,N,Lat,E,M,T);
            Eret = get<0>(tuplevals); Mret = get<1>(tuplevals);
            Cvret = get<2>(tuplevals); chiret = get<3>(tuplevals);
            outfile << T << " " << Eret << " " << Mret << " " << Cvret << " " << chiret << endl;
            cout << T << endl;
        }
        outfile.close()
    }
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
    arma::mat A(L,L);
    double randomval;
    random_device rd;
    mt19937_64 gen(rd());
    uniform_real_distribution<double> RandomNumberGenerator(0.0,1.0); //Initialize a uniform PDF RNG
    for (int i = 0;i < L;i++)
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
    tuple<double,double> Init_Vals = Calc_Init_Vals(A,L);
    double E_init = get<0>(Init_Vals); double M_init = get<1>(Init_Vals);
    return make_tuple(E_init,M_init,A);
}



tuple<double,double> Calc_Init_Vals(arma::mat A, int L)
{
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
    double E_tot = 0; double E2_tot = 0;
    double M_tot = M; double M2_tot = M*M;
    double abs_M_tot = abs(M_tot); double abs_M2_tot = abs(M_tot)*abs(M_tot);

    random_device rand; mt19937_64 gen(rand()); uniform_int_distribution<int> RandintGen(0,L-1); uniform_real_distribution<double> RandFloatGen(0,1);

    int Row; int Column; double DeltaE; double DeltaM; double MetroReq; double Abs_M = 0;

    arma::vec DeltaEList(5);
    DeltaEList(0) = exp(8/Temp); DeltaEList(1) = exp(4/Temp);
    DeltaEList(2) = 1;
    DeltaEList(3) = exp(-4/Temp); DeltaEList(4) = exp(-8/Temp);
    for (int i = 1;i <= N;i++)
    {
        for (int j = 0;j < L*L;j++)
        {
            Row = RandintGen(gen); Column = RandintGen(gen);
            DeltaE = 0; DeltaM = 0;
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
        E_tot += E; E2_tot += E*E;
        M_tot += M; M2_tot += M*M;
        abs_M_tot += Abs_M; abs_M2_tot += Abs_M*Abs_M;
    }
    double mean_E = E_tot/(N+1);
    double mean_M = M_tot/(N+1);
    double mean_E2 = E2_tot/(N+1);
    double mean_M2 = M2_tot/(N+1);
    double absmean_M = abs_M_tot/(N+1);
    double absmean_M2 = abs_M2_tot/(N+1);

    double cv = (mean_E2-(mean_E*mean_E))/(Temp*Temp);
    double chi = (mean_M2-(mean_M*mean_M))/Temp;
    cout << "<E> = " << mean_E/(double)(L*L) << endl;
    cout << "<M> = " << mean_M/(double)(L*L) << endl;
    cout << "<|M|> = " << absmean_M/(double)(L*L) << endl;
    cout << "cv = " << cv/(double)(L*L) << endl << "chi = " << chi/(double)(L*L) << endl;
}

void flip_n_calc_equilibrium(int L, int N, arma::mat Lat, double E, double M, double Temp, string filename)
{
    ofstream outfile;
    outfile.open(filename+".txt");

    double E_tot = E; double M_tot = M; double Abs_M_tot = abs(M);
    outfile << "<E> - <M> - <accepted> - <|M|>" << endl;

    random_device rand; mt19937_64 gen(rand()); uniform_int_distribution<int> RandintGen(0,L-1); uniform_real_distribution<double> RandFloatGen(0,1);

    int Row; int Column; double DeltaE; double DeltaM; double MetroReq; double Abs_M;

    arma::vec DeltaEList(5);
    DeltaEList(0) = exp(8/Temp); DeltaEList(1) = exp(4/Temp);
    DeltaEList(2) = 1;
    DeltaEList(3) = exp(-4/Temp); DeltaEList(4) = exp(-8/Temp);
    for (int i = 1;i <= N;i++)
    {
        for (int j = 0;j < L*L;j++)
        {
            Row = RandintGen(gen); Column = RandintGen(gen);
            DeltaE = 0; DeltaM = 0;
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
        M_tot += M; Abs_M_tot += Abs_M;



        outfile << E_tot/(L*L*(i+1)) << " " << M_tot/(L*L*(i+1)) << " " << Abs_M_tot/(L*L*(i+1)) << endl;

    }
    outfile.close();
}

void flip_n_calc_probdist(int L, int N, arma::mat Lat, double E, double M, double Temp, string filename)
{
    ofstream outfile;
    outfile.open(filename+".txt");

    double E_tot = E; double E2_tot = E*E;

    random_device rand; mt19937_64 gen(rand()); uniform_int_distribution<int> RandintGen(0,L-1); uniform_real_distribution<double> RandFloatGen(0,1);

    int Row; int Column; double DeltaE; double MetroReq; double Abs_M; double E_mean;

    arma::vec DeltaEList(5);
    DeltaEList(0) = exp(8/Temp); DeltaEList(1) = exp(4/Temp);
    DeltaEList(2) = 1;
    DeltaEList(3) = exp(-4/Temp); DeltaEList(4) = exp(-8/Temp);
    for (int i = 1;i <= N;i++)
    {
        for (int j = 0;j < L*L;j++)
        {
            Row = RandintGen(gen); Column = RandintGen(gen);
            DeltaE = 0;
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
            E += DeltaE;
        }
        E_tot += E; E2_tot += E*E;
        if (i > 1000000)
        {
            outfile << E << endl;
        }

    }
    outfile << E2_tot/(N+1) << endl;
    outfile.close();
}

tuple<double,double,double,double> flip_n_calc_ordinary_output(int L, int N, arma::mat Lat, double E, double M, double Temp)
{
    double E_tot = 0; double E2_tot = 0;
    double M_tot = M; double M2_tot = M*M;
    double abs_M_tot = abs(M_tot); double abs_M2_tot = abs(M_tot)*abs(M_tot);

    random_device rand; mt19937_64 gen(rand()); uniform_int_distribution<int> RandintGen(0,L-1); uniform_real_distribution<double> RandFloatGen(0,1);

    int Row; int Column; double DeltaE; double DeltaM; double MetroReq; double Abs_M = 0;

    arma::vec DeltaEList(5);
    DeltaEList(0) = exp(8/Temp); DeltaEList(1) = exp(4/Temp);
    DeltaEList(2) = 1;
    DeltaEList(3) = exp(-4/Temp); DeltaEList(4) = exp(-8/Temp);
    for (int i = 1;i <= N;i++)
    {
        for (int j = 0;j < L*L;j++)
        {
            Row = RandintGen(gen); Column = RandintGen(gen);
            DeltaE = 0; DeltaM = 0;
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
        if(i > 1000000)
        {
            E_tot += E; E2_tot += E*E;
            M_tot += M; M2_tot += M*M;
            abs_M_tot += Abs_M; abs_M2_tot += Abs_M*Abs_M;
        }
    }
    double mean_E = E_tot/(N+1-1000000);
    double mean_M = M_tot/(N+1-1000000);
    double mean_E2 = E2_tot/(N+1-1000000);
    double mean_M2 = M2_tot/(N+1-1000000);
    double absmean_M = abs_M_tot/(N+1-1000000);
    double absmean_M2 = abs_M2_tot/(N+1-1000000);

    double cv = (mean_E2-(mean_E*mean_E))/(Temp*Temp);
    double chi = (mean_M2-(mean_M*mean_M))/Temp;

    return make_tuple(mean_E/(double)(L*L),absmean_M/(double)(L*L),cv/(double)(L*L),chi/(double)(L*L));
}
