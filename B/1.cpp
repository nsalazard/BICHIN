// Simular el movimiento de N moleculas en un gas 2D
#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

//---------- Constantes --------
const int P = 8;
const int L = 1000;
const int Ni = 2;
const int K = 10;
const double TMAX = 100;
//--- ------ Clases ------------
class Bichin;
class Colisionador;

//---- interface e implementacion de clases ----
//---- clase Bichin ---
class Bichin{
private:
  double x,y,m, E, H;
	double moves[P];
public:

  void Move(double prob, double K);
  void Feed(double food){E += food;};
	void Inicie(double x0,double y0,double E0,double m0);
  void Dibujese(void);
	void Birth(double prob);
  double Getx(void){return x;}; //inline
  double Gety(void){return y;}; //inline
	double GetE(void){return y;}; //inline
  friend class Colisionador;
};

void Bichin::Inicie(double x0,double y0,double E0,double m0){
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

void Bichin::Move(double prob, double K){
	double min=0.0,max=0.0;
  for(int ii= 0; ii < (P-1); ii++ ){
		min += moves[ii]; max = min + moves[ii+1];
		if(prob >= min && prob <= max){
			x += K*std::cos(ii*(M_PI/4)); y += K*std::sin(ii*(M_PI/4));
			E-=1;
		}
	}
}

void Bichin::Birth(double prob){
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
  cout<<" , "<<x<<"+"<<5<<"*cos(t),"<<y<<"+"<<5<<"*sin(t)";
}
//--- clase Colisionador ----
class Food{
private:
  double x,y, E;
public:
  void Inicie(void);
  void CalculeFuerzas(Bichin * Grano,int Nlive, double dt);
  void CalculeFuerzaEntre(Bichin & Grano1,Bichin & Grano2
			  ,double & x_Cundall,double & s_old,double dt);
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
  Crandom ran64(1);
	
	InicieAnimacion(); //Dibujar
  
  //Inicializar los bichines
  Bichitos[0].Inicie(100, 100, 100, 1);
	Bichitos[1].Inicie(-100, -100, 100, 1);

  for(double t=0, tdibujo=0; t<TMAX ; t+=1){
    //Move
    for(int ii=0; ii<Ni; ii++){
			//if(Bichitos[ii].GetE()>0){
				double prob = 100*ran64.r();
				Bichitos[ii].Move(prob, K);
				//}
		}
     //if(tdibujo>tcuadro){
      InicieCuadro();
	for(int ii=0;ii<Ni;ii++) {Bichitos[ii].Dibujese();}
      TermineCuadro();
     // tdibujo=0;
    //}
  }     
  return 0;
}
