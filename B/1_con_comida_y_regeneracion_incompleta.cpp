// Simular el movimiento de N moleculas en un gas 2D
#include <iostream>
#include <cmath>
#include <vector>
#include "Random64.h"
using namespace std;

//---------- Constantes --------
const int P = 8; //Numero de parámetros
const int L = 500; //Espacio 2L*2L
const int K = 10; //
const double TMAX = 200;
const int Ni = 100; //Numero de bichines 
const int E_Total = 25000;  //Energia total del sistema
const int Efood = 50;    //Valor energetico de la comida
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
  friend class Food;
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
	    x += K*cos(ii*(M_PI/4)); y += K*sin(ii*(M_PI/4));
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
  double x,y,R,E;
  int murio;  //determina si la comida fue consumida o no
public:
  void Inicie(double x0,double y0,double E0, double R0, int murio); 
  void Feed(Bichin & Bicho); //Interacción bichines-comida
  void Dibujese(void);
  //void Regenere(int reg)
  int Getmurio(void){return murio;};
};

void Food::Inicie(double x0,double y0,double E0, double R0, int murio){
  murio = 0;  //inicialmente toda la comida esta disponible
  x=x0; y=y0;
  R=R0; E=E0;
}

void Food::Feed(Bichin & Bicho){
  double dis=sqrt(pow(Bicho.x - x,2.0) + pow(Bicho.y - y,2.0));
  if(dis < Bicho.R){
    Bicho.E += E;
    murio =1;
  }
}

void Food::Dibujese(void){
  cout<<" , "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}

//----------------- Funciones de Animacion ----------
void InicieAnimacion(void){
  cout<<"set terminal gif animate"<<endl; 
  cout<<"set output 'Generacion_comida.gif'"<<endl;
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
  Food Fruta[E_Total/Efood]; //Generacion de comida
  Crandom ran64(1);
  double R = 10.0;
  double Rf = 3.0; //Radio de la comida
  double E0_bichin = 300;   //Energia inicial de los bichines
  int Ehijos = 100;
  double Thijos = 20;
  double prob1;
  int prob;
  int reg;
  // (E_total - E_bichines - E_comida_viva)/E_food ----> numero de frutas que regenerar 
  
  InicieAnimacion(); //Dibujar
  
  //Inicializar los bichines
             //Inicie(x0,y0,E0, m0, R0)
  Bichitos[0].Inicie(100, 100, E0_bichin, 1,R);
  Bichitos[1].Inicie(-100, -100, E0_bichin, 1,R);

  //Inicializar la comida
    for(int j=0; j < E_Total/Efood ; j++){
      double rx= (2*ran64.r() -1)*L; double ry= (2*ran64.r() -1)*L;
      //       Inicie(x0,y0,E0, R0, murio);
      if(j < (E_Total - Nlive*E0_bichin)/Efood){
	int vivo = 0;
	 Fruta[j].Inicie(rx,ry,Efood,Rf,vivo);
      } else {
	int muerto = 1;
	Fruta[j].Inicie(rx,ry,Efood,Rf,muerto);
      }	
    }

  int en_comida = 0;
  int en_bichines = 0;

  
  for(double t=0, tdibujo=0; t<TMAX ; t+=1){  //Move
      for(int ii=0; ii<Nlive; ii++){
	if(Bichitos[ii].GetE() > 0){//Si el bichin esta vivo
	  en_bichines += Bichitos[ii].GetE();
	  prob1 = 100*ran64.r();
	  Bichitos[ii].Move(K, prob1);  //que se mueva
          for(int i=0; i<E_Total/Efood; i++){ //que se alimente
	    if(Fruta[i].Getmurio() == 0) {
	      if(ii==0){en_comida += Efood;} //solo se suma comida la 1ra vez
	      Fruta[i].Feed(Bichitos[ii]);
	    }

	    if(ii==1){
	      reg = (E_Total - en_bichines - en_comida)/Efood;
	      if(Fruta[i].Getmurio() == 1 && reg > 0){
		//guardar indices en un vector
	      }
	    }
	  }
	   if(Bichitos[ii].GetE()>Ehijos && t > Thijos){//y se reproduzca si puede
	    prob = int(P*ran64.r());
	    Fate.Birth(Bichitos[ii], Bichitos[Nlive], t, prob);
	    }
	}
      }

      /*Funcion regenere:
	Se regenera la comida necesaria (int reg) empleando los indices de la comida muerta
       */
     
    InicieCuadro();
    for(int i=0; i<E_Total/Efood; i++){  //Solo dibuje la comida que no ha sido consumida
      if(Fruta[i].Getmurio() == 0) {Fruta[i].Dibujese();}
    }
    for(int ii=0;ii<Nlive;ii++) {Bichitos[ii].Dibujese();}
    TermineCuadro();

    en_bichines = 0;
    en_comida = 0;
  }     
  return 0;
}
