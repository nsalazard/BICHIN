// Simular el movimiento de N moleculas en un gas 2D
#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

//---------- Constantes --------
const int P = 8; //Numero de parámetros
const int L = 1000; //Espacio 2L*2L
const int K = 10;
const double TMAX = 100;
const int Ni = 1000;
int Nlive = 2;
//--- ------ Clases ------------
class Bichin;
class Selection;

//---- interface e implementacion de clases ----
//---- clase Bichin ---
class Bichin{
private:
  double x,y,m, E, H, R;
	double moves[P];
public:
  void Move(double K, Crandom ran64);
  void Feed(double food){E += food;};
	void Inicie(double x0,double y0,double E0,double m0, double R0);
  void Dibujese(void);
  double Getx(void){return x;}; //inline
  double Gety(void){return y;}; //inline
	double GetE(void){return E;}; //inline
	bool Alive(void){return E>0.0;};
  friend class Selection;
};

void Bichin::Inicie(double x0,double y0,double E0,double m0, double R0){
	R = R0;
	E = E0;
  x = x0; y = y0;
	moves[0] = 12.5; // East
	moves[1] = 12.5; // North-East
	moves[2] = 12.5; // North
	moves[3] = 15; // Nort-West
	moves[4] = 10; // West
	moves[5] = 15; // South-West
	moves[6] = 10; // South
	moves[7] = 12.5; // South-East
} 

void Bichin::Move(double K, Crandom ran64){
	double prob = 100*ran64.r();
	double min=0.0,max=0.0;
  for(int ii= 0; ii < (P-1); ii++ ){
		min += moves[ii]; max = min + moves[ii+1];
		if(prob >= min && prob <= max){
			x += K*std::cos(ii*(M_PI/4)); y += K*std::sin(ii*(M_PI/4));
			E-=1;
		}
	}
}

void Bichin::Dibujese(void){
  cout<<" , "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}
//--- clase Colisionador ----
class Selection{
private:
public:
  void Inicie(void); //Inicializa todos los individuos
  void CalculeFuerzas(Bichin * Grano,int Nlive, double dt);
  void Birth(Bichin & BichoP,Bichin & BichoH, double t, int prob);
};

void Selection::Birth(Bichin & BichoP,Bichin & BichoH, double t,int prob){
		Nlive +=1;
		BichoP.E = BichoP.E/2;
		BichoH.Inicie(0, 0, BichoP.E, 1, BichoP.R);
		//Mutation
	for(int ii=0; ii < P; ii++){
		BichoH.moves[ii] = BichoP.moves[ii];
	}
	BichoH.moves[prob] -=1; 
	//cout << BichoH.GetE() << "\n";
}


//---------FOOD-------------------------------------
class Food{
private:
  double x,y, E;
public:
  void Inicie(void);
};

//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output '1.gif'"<<endl;
  cout<<"unset key"<<endl;
  cout<<"set xrange["<<-L<<":"<<L<<"]"<<endl;
  cout<<"set yrange["<<-L<<":"<<L<<"]"<<endl;
  cout<<"set size ratio -1"<<endl;
  cout<<"set parametric"<<endl;
  cout<<"set trange [0:7]"<<endl;
  cout<<"set isosamples 12"<<endl;  
}
void InicieCuadro(void){
    cout<<"plot 0,0 ";
}
void TermineCuadro(void){
    cout<<endl;
}
//-----------  Programa Principal --------------  
int main(void){
  Bichin Bichitos[Ni];
	Selection Fate;
  Crandom ran64(1);
	double R = 10.0;
	int Ehijos = 100;
	double Thijos = 20;
	double prob1;
	int prob;
	
	InicieAnimacion(); //Dibujar
  
  //Inicializar los bichines
  Bichitos[0].Inicie(100, 100, 300, 1,R);
	Bichitos[1].Inicie(-100, -100, 100, 1,R);

  for(double t=0, tdibujo=0; t<TMAX ; t+=1){  //Move  
    for(int ii=0; ii<Nlive; ii++){
			if(Bichitos[ii].GetE() > 0){ //Bichitos[ii].Alive()
				prob1 = 100*ran64.r();
				Bichitos[ii].Move(K, prob1);
				if(Bichitos[ii].GetE()>Ehijos && t > Thijos){
				prob = int(P*ran64.r());
				Fate.Birth(Bichitos[ii], Bichitos[Nlive], t, prob);
					}
			}
		}
     //if(tdibujo>tcuadro){
      InicieCuadro();
			for(int ii=0;ii<Nlive;ii++) {Bichitos[ii].Dibujese();}
      TermineCuadro();
     // tdibujo=0;
    //}
  }     
  return 0;
}
