// Simular el movimiento de N moleculas en un gas 2D
#include <iostream>
#include <cmath>
#include "Random64.h"
using namespace std;

//---------- Constantes --------
const int P = 8; //Numero de parámetros de los bichines
const int L = 400; //Espacio 2L*2L
const int K = 10;
const double TMAX = 400;
const int Ni = 10000;
const int Nfood = 2;
int Nlive = 2;
int Energy = 0;
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
  void Move(double K, double prob);
  //void Feed(double food){E += food;};
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
	//Genetica
	moves[0] = 12.5; // East
	moves[1] = 12.5; // North-East
	moves[2] = 12.5; // North
	moves[3] = 15; // Nort-West
	moves[4] = 10; // West
	moves[5] = 15; // South-West
	moves[6] = 10; // South
	moves[7] = 12.5; // South-East
} 

void Bichin::Move(double K, double prob){
	//double prob = 100*ran64.r();
	double min=0.0,max=0.0;
  for(int ii= 0; ii < (P-1); ii++ ){
		min += moves[ii]; max = min + moves[ii+1];
		if(prob >= min && prob <= max){
			x += K*std::cos(ii*(M_PI/4)); y += K*std::sin(ii*(M_PI/4));
			E-=1;
			Energy+=1;
		}
	}
}

void Bichin::Dibujese(void){
  cout<<" , "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}

//---------FOOD-------------------------------------
class Food{
private:
  double x,y, E, R, val;
public:
  void Inicie(double x0,double y0,double E0, double R0, double val0);
	void Feed(Bichin & Bicho);
	void Dibujese(void);
	void Recharge(double newEnergy){E += newEnergy; Energy =0;};
	friend class Selection;
  friend class Bichin;
};
void Food::Inicie(double x0,double y0,double E0, double R0, double val0){
	x = x0; y = y0;
	R = R0;
	E = E0;
  val = val0;
	}

void Food::Feed(Bichin & Bicho){
	double dis = sqrt(pow(Bicho.x-x,2.0)+pow(Bicho.y-y,2.0));
	if (dis < (R+ Bicho.R)){
		//cout << Bicho.E << "\n";
		Bicho.E += val;
		E -= val;
		//cout << Bicho.E << "\n";
	}
}

void Food::Dibujese(void){
  cout<<" , "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}


//--- clase Selection----
class Selection{
private:
public:
  void Birth(Bichin & BichoP,Bichin & BichoH, double t, int prob, int prob2);
	void RechargeFood(Food *food);
friend class Food;
};

void Selection::Birth(Bichin & BichoP,Bichin & BichoH, double t,int prob,int prob2){
		Nlive +=1;
		BichoP.E = BichoP.E/2;
		BichoH.Inicie(BichoP.x, BichoP.y, BichoP.E, 1, BichoP.R);
		//Mutation
	for(int ii=0; ii < P; ii++){
		BichoH.moves[ii] = BichoP.moves[ii];
	}
	BichoH.moves[prob] -=1;
	BichoH.moves[(prob2)%7] +=1;
	//cout << BichoH.GetE() << "\n";
}

void  Selection::RechargeFood(Food *food){
	double EnergyP= Energy/Nfood;
	for(int ii=0; ii < Nfood; ii++){
		food[ii].E += EnergyP;
	}
	Energy = 0;
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
	Food food[Nfood];
	Selection Fate;
  Crandom ran64(1);
	double R = 10.0;
	int Ehijos = 100; //Min Energy for reproduction
	double Thijos = 20; //Min Time for reproduction
	double Rfood = 10;
	double prob;
	int prob1,prob2;
	
  //Inicializar los bichines
  Bichitos[0].Inicie(0, 0, 300, 1,R);
	Bichitos[1].Inicie(0, 0, 100, 1,R);

	//Inicializar la comida
	food[0].Inicie(-50,0,10000,Rfood,5);
	food[1].Inicie(50,0, 10000, Rfood,9);

	InicieAnimacion(); //Dibujar

  for(double t=0, tdibujo=0; t<TMAX ; t+=1){  //Move  
    for(int ii=0; ii<Nlive; ii++){
			if(Bichitos[ii].Alive()){ //Bichitos[ii].Alive()
				prob = 100*ran64.r();
				Bichitos[ii].Move(K, prob);
				if(Bichitos[ii].GetE()>Ehijos && t > Thijos){
				prob1 = int(P*ran64.r());
				prob2 = int(P*ran64.r());
				Fate.Birth(Bichitos[ii], Bichitos[Nlive], t, prob1, prob2);
					}
				// Get some food
				for(int jj=0; jj< Nfood; jj++){
					food[jj].Feed(Bichitos[ii]);
					}
			}
		}
     //if(tdibujo>tcuadro){
      InicieCuadro();
			for(int ii=0;ii<Nlive;ii++){
				if(Bichitos[ii].Alive()){
					Bichitos[ii].Dibujese();
					}
				}
			for(int ii=0;ii<Nfood;ii++){
				food[ii].Dibujese();
				}
      TermineCuadro();
     // tdibujo=0;
    //}
		
	//Recharge Food
		Fate.RechargeFood(food);
		
  }     
  return 0;
}
