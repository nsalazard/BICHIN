// Simulación de Bichines
#include <iostream>
#include <fstream>
#include <cmath>
#include "Random64.h"
using namespace std;

//---------- Constantes --------
const int P = 8;   // Numero de parámetros de los bichines
const int L = 200; // Espacio 2L*2L
const double K = 10;
const double TMAX = 5;
const int Ni = 10000;
const int Nfood = 5; // 200;
int Nlive = 2;       // 17;
double Energy = 0;
//--- ------ Clases ------------
class Bichin;
class Selection;

//---- interface e implementacion de clases ----
//---- clase Bichin ---
class Bichin
{
private:
    int t_live;
    double x = 0.0, y = 0.0, m = 0.0, E = 0.0, H, R = 0.0;
    int alive = 0;
    double moves[P];

public:
    void Move(double K, double prob);
    // void Feed(double food){E += food;};
    void Start(double x0, double y0, double E0, double m0, double R0, Crandom &ran64);
    double Getx(void) { return x; };      // inline
    double Gety(void) { return y; };      // inline
    double GetE(void) { return E; };      // inline
    double GetT(void) { return t_live; }; // inline
    void Genetic(Crandom &ran64);
    void Print(void);
    int Alive(void);
    friend class Selection;
    friend class Food;
};

void Bichin::Start(double x0, double y0, double E0, double m0, double R0, Crandom &ran64)
{
    R = R0;
    E = E0;
    x = x0;
    y = y0;
    t_live = 0;
    alive = 1;
    // Genetica
    Genetic(ran64);
}

void Bichin::Genetic(Crandom &ran64)
{
    int xx;
    double sum = 0.0;
    for (int ii = 0; ii < P; ii++)
    {
        xx = ran64.r() * 10.0;
        moves[ii] = xx;
        sum += moves[ii];
    }
    for (int ii = 0; ii < P; ii++)
    {
        moves[ii] = moves[ii] / sum;
    }
}

void Bichin::Move(double K, double prob)
{
    // double prob = 100*ran64.r();
    double min = 0.0, max = 0.0;
    for (int ii = 0; ii < P; ii++)
    {
        max = min + moves[ii];
        if (prob >= min && prob <= max)
        {
            x += (K * std::cos(ii * (M_PI / 4)));
            if (x >= L)
            {
                x = -L;
            }
            else if (x <= -L)
            {
                x = L;
            }

            y += (K * std::sin(ii * (M_PI / 4)));
            if (y >= L)
            {
                y = -L;
            }
            else if (y <= -L)
            {
                y = L;
            }
            E -= 1;
            Energy += 1;
            t_live += 1;
        }
        // Damos el min
        min += moves[ii];
    }
}

void Bichin::Print(void)
{
    cout << " , " << x << "+" << R << "*cos(t)," << y << "+" << R << "*sin(t)";
}

int Bichin::Alive(void)
{
    if (E > 0.0)
    {
        alive = 1;
    }
    else
    {
        alive = 0;
    }
    return alive;
}

//---------FOOD-------------------------------------
class Food
{
private:
    double x, y, E, R, val;

public:
    double Getx(void) { return x; };
    double Gety(void) { return y; };
    double GetE(void) { return E; };
    void Start(double x0, double y0, double E0, double R0, double val0);
    void Feed(Bichin &Bicho);
    void Print(void);
    void Blender(void);
    void Recharge(double newEnergy)
    {
        E += newEnergy;
        Energy = 0;
    };
    bool Alive(void) { return E > 0.0; };
    friend class Selection;
    friend class Bichin;
};
void Food::Start(double x0, double y0, double E0, double R0, double val0)
{
    x = x0;
    y = y0;
    R = R0;
    E = E0;
    val = val0;
}

void Food::Feed(Bichin &Bicho)
{
    double dis = sqrt(pow(Bicho.x - x, 2.0) + pow(Bicho.y - y, 2.0));
    if (dis < (R + Bicho.R) && E > 0.0)
    {
        // cout << Bicho.E << "\n";
        Bicho.E += val;
        E -= val;
        // cout << Bicho.E << "\n";
    }
}

void Food::Print(void)
{
    cout << " , " << x << "+" << R << "*cos(t)," << y << "+" << R << "*sin(t)";
}

//--- clase Selection----
class Selection
{
private:
public:
    void Birth(Bichin &BichoP, Bichin &BichoH, int prob, int prob2, Crandom &ran64);
    void Spread(Food *food, int N, double mu, double sigma, double Rfood, Crandom &ran64);
    void RechargeFood(Food *food, Crandom &ran64);
    friend class Food;
};

void Selection::Birth(Bichin &BichoP, Bichin &BichoH, int prob, int prob2, Crandom &ran64)
{
    Nlive += 1;
    BichoH.alive = 1;
    BichoP.E = BichoP.E / 2;
    BichoH.Start(BichoP.x, BichoP.y, BichoP.E, 1, BichoP.R, ran64);
    // Mutation
    for (int ii = 0; ii < P; ii++)
    {
        BichoH.moves[ii] = BichoP.moves[ii];
    }
    BichoH.moves[prob] -= 0.01;
    BichoH.moves[(prob2) % 7] += 0.01;
    // cout << BichoH.GetE() << "\n";
}

void Selection::Spread(Food *food, int N, double mu, double sigma, double Rfood, Crandom &ran64)
{
    int ix, iy;
    double E_inicial = 50;
    double val = 5;
    for (int ii = 0; ii < N; ii++)
    {
        // Escogemos celdas al azar segun una distribucion bigaussiana
        ix = (int)ran64.gauss(mu, sigma);
        iy = (int)ran64.gauss(mu, sigma);
        // Evitamos las fronteras del array
        if (ix < -L)
            ix = -L;
        if (iy < -L)
            iy = -L;
        if (ix > (L))
            ix = L;
        if (iy > L)
            iy = L;
        food[ii].Start(ix, iy, E_inicial, Rfood, val);
    }
}

void Selection::RechargeFood(Food *food, Crandom &ran64)
{
    // cout << "Comida recargada "<<( Energy) << "\n";
    int prob;
    while (Energy > 0)
    {
        prob = int(Nfood * ran64.r());
        // cout << "Comida antes"<< food[prob].E  << "\n";
        food[prob].E += 1;
        // cout << "Comida despues"<< food[prob].E  << "\n";
        Energy -= 1;
    }
};

//----------------- Funciones de Animacion ----------
void StartAnimacion(void)
{
    cout << "set terminal gif animate" << endl;
    cout << "set output 'Blender.gif'" << endl;
    cout << "unset key" << endl;
    cout << "set xrange[" << -L << ":" << L << "]" << endl;
    cout << "set yrange[" << -L << ":" << L << "]" << endl;
    cout << "set size ratio -1" << endl;
    cout << "set parametric" << endl;
    cout << "set trange [0:7]" << endl;
    cout << "set isosamples 12" << endl;
}
void InicieCuadro(void)
{
    cout << "plot 0,0 ";
}
void TermineCuadro(void)
{
    cout << endl;
}
//-----------  Programa Principal --------------
int main(void)
{
    Bichin Bichitos[Ni];
    Food food[Nfood];
    Selection Fate;
    Crandom ran64(1);
    double R = 4.0;
    int Ehijos = 100;   // Min Energy for reproduction
    double Thijos = 20; // Min Time for reproduction
    double Rfood = 2;
    double prob;
    double mu = 0.0, sigma = L / 4;
    int prob1, prob2;

    // Files Blender
    ofstream Bichin;
    ofstream Food;
    Bichin.open("Bichin.txt");
    Food.open("Food.txt");

    // Inicializar los bichines
    for (int jj = 0; jj < Nlive; jj++)
    {
        Bichitos[jj].Start(0.0, 0.0, 120.0, 1.0, R, ran64);
    }
    Fate.Spread(food, Nfood, mu, sigma, Rfood, ran64);
    // Animation
    // StartAnimacion(); //Dibujar

    for (int t = 0, tdibujo = 0; t < TMAX; t += 1)
    { // Move
        for (int ii = 0; ii < Nlive; ii++)
        {
            if (Bichitos[ii].Alive())
            { // Bichitos[ii].Alive()

                prob = ran64.r();
                Bichitos[ii].Move(K, prob);
                if (Bichitos[ii].Alive() == 1 && Bichitos[ii].GetE() > Ehijos && Bichitos[ii].GetT() > Thijos)
                {
                    prob1 = int(P * ran64.r());
                    prob2 = int(P * ran64.r());
                    Fate.Birth(Bichitos[ii], Bichitos[Nlive], prob1, prob2, ran64);
                }
                // Get some food
                for (int jj = 0; jj < Nfood; jj++)
                {
                    food[jj].Feed(Bichitos[ii]);
                }
            }
        }
        // if(tdibujo>tcuadro){
        Bichin << t << "\t";
        Food << t << "\t";
        // InicieCuadro();
        for (int ii = 0; ii < 4; ii++)
        {
            // Bichitos[ii].Print();
            Bichin << "\t" << Bichitos[ii].Getx() << "\t" << Bichitos[ii].Gety() << "\t" << Bichitos[ii].Alive();
        }
        for (int ii = 0; ii < Nfood; ii++)
        {
            // food[ii].Print();
            Food << "\t" << food[ii].Getx() << "\t" << food[ii].Gety() << "\t" << food[ii].Alive();
        }
        // TermineCuadro();
        Bichin << "\n";
        Food << "\n";
        // tdibujo=0;
        //}

        //---------------- PROOFS-----------

        // Recharge Food
        Fate.RechargeFood(food, ran64);
    }
    // Bichin.close();
    // Food.close();
    return 0;
}
