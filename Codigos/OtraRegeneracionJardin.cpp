 // Simular el movimiento de N moleculas en un gas 2D
#include <iostream>
#include <fstream>
#include <cmath>
#include "Random64.h"
using namespace std;

//---------- Constantes --------
const int P = 8; //Numero de parámetros de los bichines
const int L = 100; //Espacio 2L*2L
const double K = 10;
const double TMAX =1000;
const int Ni = 10000;
const int E_Total = 1200;  //Energia total del sistema
const int Efood = 10; 
int Nlive = 3;
int Energy = 0;
const int Nfood = 1200;
const double Thijos = 20; //Min Time for reproduction
double Rbichin = 5.0;
int E0_bichin=150;
int Ehijos = 150; //Min Energy for reproduction
double Rfood = 2.0 ;
//--- ------ Clases ------------
class Bichin; 
class Selection;

//---- interface e implementacion de clases ----
//---- clase Bichin ---
class Bichin{
private:
  double x,y,m, E, H, R,time;
  double moves[P];
public:
	int Move(double K, double prob);
	//void Feed(double food){E += food;};
  void Inicie(double x0,double y0,double E0,double m0, double R0, Crandom & ran64);
	void Dibujese(const char * NameFile);
	double Getx(void){return x;}; //inline
	double Gety(void){return y;}; //inline
	double GetE(void){return E;}; //inline
	void Genetic(Crandom & ran64);
	bool Alive(void){return E>0.0;};
        void Time(void);
        double GetTime(void){return time;};
	friend class Selection;
	friend class Food;
};
void Bichin::Time(void){
    time+=1;
    if(time>21){
      time=0;}}
void Bichin::Inicie(double x0,double y0,double E0,double m0, double R0, Crandom & ran64){
	R = R0;
	E = E0;
    x = x0; y = y0;
    time=0.0;
	//Genetica
	Genetic(ran64);
} 

void Bichin::Genetic(Crandom & ran64){
	int xx;
	double sum = 0.0;
	for(int ii=0; ii<P; ii++){
		xx = ran64.r()*10.0;
		moves[ii] = xx;
		sum += moves[ii];
	} 
	for(int ii=0; ii<P; ii++){
		moves[ii] = moves[ii]/sum;
	} 
	
	
/*for(int ii=0; ii<P; ii++){
		cout << moves[ii]<< endl;
	} */

	
} 


int Bichin::Move(double K, double prob){
	//double prob = 100*ran64.r();
	double paso=0;
	double min=0.0,max=0.0;
	for(int ii= 0; ii < (P-1); ii++ ){
		min += moves[ii]; max = min + moves[ii+1];
		if(prob >= min && prob <= max){
			x += (K*std::cos(ii*(M_PI/4)));
			if(x >= L) { x = -L;}
			else if(x <= -L) { x = L;}
			y += K*std::sin(ii*(M_PI/4));
			if(y >= L) { y = -L;}
			else if(y <= -L) { x = L;}
			E-=1;
			Energy+=1;
			return 1;

		}
	}
	return 0;
}

void Bichin::Dibujese(const char * NameFile){
  ofstream MyFile(NameFile,fstream::out | fstream::app);
  MyFile<<" , "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}

//---------FOOD-------------------------------------
class Food{
private:
	double x,y, E, R,val;
	int Alive;
public:
	void Inicie(double x0,double y0,double E0, double R0, double val0,int live);
	bool Feed(Bichin & Bicho);
	void Dibujese(const char * NameFile);
	void Recharge(double newEnergy){E += newEnergy; Energy =0;};
	void IsAlive(bool IsAlive){if(IsAlive){
		Alive=1;}
		else{Alive=0;}};
        int GetAlive(void){return Alive;}
	friend class Selection;
	friend class Bichin;
};
void Food::Inicie(double x0,double y0,double E0, double R0, double val0,int live){
	x = x0; y = y0;
	R = R0;
	E = E0;
	val=val0;
	Alive = live;
	}

bool Food::Feed(Bichin & Bicho){
	double dis = sqrt(pow(Bicho.x-x,2.0)+pow(Bicho.y-y,2.0));
	if (dis < (R+ Bicho.R)){
		
		Bicho.E += val;
		E -= val;
		
		Alive=0;
		return true;
	}
	else {return false;}
}

void Food::Dibujese(const char * NameFile){
  ofstream MyFile(NameFile,fstream::out | fstream::app);
  
  MyFile<<" , "<<x<<"+"<<R<<"*cos(t),"<<y<<"+"<<R<<"*sin(t)";
}


//--- clase Selection----
class Selection{
private:
public:
	void Birth(Bichin & BichoP,Bichin & BichoH, double t, int prob, int prob2, Crandom & ran64);
	void RechargeFood(Food *food,Bichin *Bicho);
	friend class Food;
};

void Selection::Birth(Bichin & BichoP,Bichin & BichoH, double t,int prob,int prob2, Crandom & ran64){
	Nlive +=1;
	BichoP.E = BichoP.E/2;
	BichoH.Inicie(BichoP.x, BichoP.y, BichoP.E, 1, BichoP.R,  ran64);
		//Mutation
	for(int ii=0; ii < P; ii++){
		BichoH.moves[ii] = BichoP.moves[ii];
	}
	BichoH.moves[prob] -=0.01;
	BichoH.moves[(prob2)%7] +=0.01;
	
}

void  Selection::RechargeFood(Food *food, Bichin *Bicho){
	double EnergyP= Energy/Nfood;
	for(int ii=0; ii < Nfood; ii++){
		food[ii].E += EnergyP;
	}
	Energy = 0;
	};




//----------------- Funciones de Animacion ----------
void InicieAnimacion(const char * NameFile){
  ofstream MyFile(NameFile,fstream::out | fstream::app);
  MyFile<<"set terminal gif animate"<<endl; 
  MyFile<<"set output '2.gif'"<<endl;
  MyFile<<"unset key"<<endl;
  MyFile<<"set xrange["<<-L<<":"<<L<<"]"<<endl;
  MyFile<<"set yrange["<<-L<<":"<<L<<"]"<<endl;
  MyFile<<"set size ratio -1"<<endl;
  MyFile<<"set parametric"<<endl;
  MyFile<<"set trange [0:7]"<<endl;
  MyFile<<"set isosamples 12"<<endl;  
}
void InicieCuadro(const char * NameFile){
	ofstream MyFile(NameFile,fstream::out | fstream::app);
    MyFile<<"plot 0,0 ";
}
void TermineCuadro(const char * NameFile){
	ofstream MyFile(NameFile,fstream::out | fstream::app);
    MyFile<<endl;
	MyFile.close();
}
//-----------  Programa Principal --------------  
int main(void){
	fstream my_file;
	
	Bichin Bichitos[Ni];
	Food food[Nfood];
	Selection Fate;
	Crandom ran64(1);
	double sobra=0;
	bool fed;
	int tcuadro=2;
	
	double prob;
	int prob1,prob2;
	if(Nlive*E0_bichin>E_Total){
	    cout<<"No hay suficiente energia en el sistema para iniciar a los bichines, cambiar parametros"<<endl; return 0;}
	else{
  //Inicializar los bichines
	for (int jj=0; jj<Nlive; jj++ ){
	  
	   
	 
	    Bichitos[jj].Inicie(0, 0, E0_bichin, 1,Rbichin, ran64);
	}

	
 
	//Inicializar la comida
    for(int j=0; j < Nfood ; j++){
      double rx= (2*ran64.r() -1)*L/4; double ry= (2*ran64.r() -1)*L/4-3*L/4;
      //       Inicie(x0,y0,E0, R0, murio);
      if(j < (E_Total - Nlive*E0_bichin)/Efood){
	int vivo = 1;
	food[j].Inicie(rx,ry,Efood,Rfood,Efood,vivo);
      } else {
	int muerto = 0;
	food[j].Inicie(rx,ry,Efood,Rfood,Efood,muerto);
      }	
    }
    int en_bichines=0;
    int en_comida=0;
    int regenerate=0;

	InicieAnimacion("Prueba.gp"); //Dibujar
	  
	int Epaso=0;//lleva cuenta de la energia que se pierde con cada paso de los bichitos   

  for(double t=0, tdibujo=0; t<TMAX ; t+=1,tdibujo+=1){  //Move
     
     for(int ii=0; ii<Nlive; ii++){
    
       if(Bichitos[ii].Alive()){ //Bichitos[ii].Alive()
	 Bichitos[ii].Time();
	 prob = ran64.r();
	 Epaso+=Bichitos[ii].Move(K, prob);
	 en_bichines += Bichitos[ii].GetE();
	 
	 for(int i=0; i<Nfood; i++){ //que se alimente
	   if(food[i].GetAlive()==1){
	     if(ii==0){ en_comida+=Efood;}
	     fed=food[i].Feed(Bichitos[ii]);
		 if(fed){
		 	double rx= (2*ran64.r() -1)*L/4; double ry= (2*ran64.r() -1)*L/4-3*L/4;
      //       Inicie(x0,y0,E0, R0, murio);
			int vivo = 1;
			food[i].Inicie(rx,ry,Efood,Rfood,Efood,vivo);}}}
	 if(Bichitos[ii].GetE()>Ehijos && Bichitos[ii].GetTime() > Thijos){
	   prob1 = int(P*ran64.r());
	   prob2 = int(P*ran64.r());
	   Fate.Birth(Bichitos[ii], Bichitos[Nlive], t, prob1, prob2, ran64);
					}
			       
	  
		}
      
		}
      /* regenerate=Epaso/Efood;
	   sobra+=Epaso%Efood;
	   cout<<"Time: "<<t<<"\t"<<"E Bicho:"<<"\t"<<en_bichines<<"\t"<<"E comida:"<<en_comida<<"\t"<<"Epaso:"<<Epaso<<"\t"<<sobra<<endl;
       if(regenerate>=1 or sobra>Efood){
		if(sobra>Efood){regenerate+=sobra/Efood;
		sobra=0;}
       while(regenerate>0){

		//cout<<regenerate<<endl;
       for(int i=0; i<Nfood; i++){ //Regenerar comida
	    
	   	if(food[i].GetAlive()==0){
			//cout<<"Iteration"<<"\t"<<i<<"\t"<<regenerate<<endl;
	   	double rx= (2*ran64.r() -1)*L; double ry= (2*ran64.r() -1)*L;
      //       Inicie(x0,y0,E0, R0, murio);
		int vivo = 1;
		food[i].Inicie(rx,ry,Efood,Rfood,Efood,vivo);
		regenerate-=1;

	 	}
		if(regenerate==0){break;}}}
       Epaso=0;
       }
       
       en_bichines=0;
       en_comida=0;
       */
     if(tdibujo>tcuadro){
       InicieCuadro("Prueba.gp");
       for(int ii=0;ii<Nlive;ii++){			  
		 	if(Bichitos[ii].Alive()){
	   		Bichitos[ii].Dibujese("Prueba.gp");
	  }
	}
       for(int ii=0;ii<Nfood;ii++){
			if(food[ii].GetAlive()){
	   		food[ii].Dibujese("Prueba.gp");}
				}
      TermineCuadro("Prueba.gp");
     
      tdibujo=0;
    }
		
	//Recharge Food
      // Fate.RechargeFood(food);
		
  }     
  
  return 0;
	}}
