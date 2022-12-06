#define _USE_MATH_DEFINES 
#include <iostream>
#include <cmath>
#include "Random64.h"
#include <fstream>
#include <exception>
#include <memory>
using namespace std;
ofstream salida;

//---------- Constantes --------
const int P = 8;            // Numero de parámetros de los bichines
const int L =700;           // Espacio 2L*2L
const double K = 10;        // Distancia recorrida en cada mov. por el bichin
const double TMAX = 1000;		// Tiempo de dibujo
const double E_gordo = 50;	// Energía a partir de la cual bichin no puede comer
const int Ni = 1000; 				// Numero maximo de bichines (?)
const int Nfood = 5000;  		// Numero de maximo comida | Nunca debe ser alcanzado.
int E_inicial = 10;  				// Energía inicial de la comida
int Nlive = 80;  						// Numero inicial de bichines
int Energy_bank = 0; 				// El banco temporal de energia
int Biome_energy=20000;			// Evita un bug con el colocamiento de la comida
														// Como buena practica  Biome_energy<Nfood*E_inicial; 
														// Podria funcionar incluso si esta condicion no se cumple pero se corre un riesgo.
int food_dis=1;	
double mu = 0.0, sigma = L / 4;  	//Parámetros distribución gaussiana de comida
//--- ------ Clases ------------
class Bichin;
class Selection;

//---- interface e implementacion de clases ----
//---- clase Bichin ---
class Bichin
	{
		private:
			int t_live;  //tiempo de vida(# de movimientos) del bichin
			double x, y, m, H, R;  
			int E=0;
			double theta=0;
			double moves[P];  //Genes del bichin

		public:
			void Move(double K, double prob) //Mueve bichin
				{    
					//SE RECIBE COMO ARGUMENTO VALOR ALEATORIO ENTRE 0 y 1
					//DE ACUERDO AL INTERVALO ENTRE 0 Y 1 AL QUE prob CORRESPONDA, SE INCREMENTA (x,y) en determinada dirección
					double min = 0.0, max = 0.0;
					for (int ii = 0; ii < P; ii++)
						{    
							if(ii>0) min += moves[ii-1]; 
							max = min + moves[ii];
							if (prob >= min && prob < max)
								{
									theta+=ii*M_PI/4;
								}
						}
						
					x += K * cos(theta);
					y += K * sin(theta);
					

					if (x>L)
						{x=-L;}
					else if(x<-L)
						{x=L;}
					if (y>L)
						{y=-L;}
					else if (y<-L)
						{y=L;}

					E -= 1;  //Disminuye la energía del bichin
					Energy_bank += 1;  //Aumenta la energía del sistema
					t_live += 1;  //Aumenta el tiempo de vida del bichin
				}
			void Start(double x0, double y0, int E0, double m0, double R0, Crandom &ran64) //Inicializa bichine
				{
					R = R0;
					E = E0;
					x = x0;
					y = y0;
					t_live = 0;
					
					Genetic(ran64);  //Inicializa su genética de forma aleatoria
				}
			double Getx(void) { return x; };      // inline
			double Gety(void) { return y; };      // inline
			double GetE(void) { return E; };      // inline
			double GetT(void) { return t_live; }; // inline
			void T_zero(void) {t_live=0;};
			void Genetic(Crandom &ran64)
				{   
					//LLena moves[] con porcentajes aleatorios tal que su suma de 1
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
					theta=2*M_PI*int(7*ran64.r());
				}
			void Print(void);
			void Blender(void);
			bool Alive(void) { return E > 0; }; //determina si vive
			friend class Selection;
			friend class Food;
	};

void Bichin::Print(void)
{
  salida << " , " << x << "+" << R << "*cos(t)," << y << "+" << R << "*sin(t)" << " lt -1";
}
void Bichin::Blender(void)
{
	salida << "\t" << x << "\t" << y;
}

//---------FOOD-------------------------------------
class Food
{
	private:
		double x, y, R;  //Posición, energía, radio y val de la comida
		int E=0;
	public:
		int GetE(void) { return E; };
		void Start(double x0, double y0, int E0, double R0)//Inicialice la comida
			{
				x = x0;
				y = y0;
				R = R0;
				E = E0; //Valor E de la comida que es cedido cuando se interacciona con el bichin 
			}
		void ReStart(int E0,Crandom &ran64)
			{	

				if(food_dis==0)
					{
						x =  2*L*ran64.r() - L;
						y =  2*L*ran64.r() - L;
					}
				if(food_dis==1)
					{
						x = (int)ran64.gauss(mu, sigma);
						y = (int)ran64.gauss(mu, sigma);		
					}
				E=E0;
			}
		void Feed(Bichin &Bicho)  //Alimente al bicho  
		{

			double dis = sqrt(pow(Bicho.x - x, 2.0) + pow(Bicho.y - y, 2.0));

			if (dis < (R + Bicho.R) && E > 0.0 && Bicho.E < E_gordo)\
				{ //Si el bicho y la comida se superponen y la energía de la comida es mayor a cero y la energía del bichin no es mayor a E_gordo
					Bicho.E += E; //Aumente la energía del bicho en val 
					E = 0; //Disminuya la energía del bicho en val
					// salida << Bicho.E << "\n";
					//cout<<"Comido \n";
				}
		}
		void Print(void);   //Imprime la comida
		void Blender(void);

		friend class Selection;
		friend class Bichin;
};

void Food::Print(void)
{
  salida << " , " << x << "+" << R << "*cos(t)," << y << "+" << R << "*sin(t)" << " lt 7";
}
void Food::Blender(void)
{
	salida << "\t" << x << "\t" << y;
}

//--- clase Selection----
class Selection
{
	private:
	public:
		void Birth(Bichin &BichoP, Bichin &BichoH, double t, int prob, int prob2, Crandom &ran64);  //Reproduce un bicho Padre en 2 bichos Hijos teniend en cuenta el tiempo de vida del padre. prob determina el gen que disminuye, prob 2 el que aumenta
		void Spread(Food *food, int N, double mu, double sigma, double Rfood, Crandom &ran64);  //Distribuya la comida con una gausiana
		void Uniform(Food *food, int N, double Rfood, Crandom &ran64);  //Distribuya la comida uniformemente
		void RechargeFood(Food *food, Crandom &ran64); //Recargue la comida de manera aleatoria en el ambiente
		int Biomass(Food *food, Bichin *Bichos)
			{
				int ii=0;
				int total_biomass=0;
				for(ii=0;ii<Ni;ii++)
					{total_biomass+=Bichos[ii].GetE();};
				for(ii=0;ii<Nfood;ii++)
					{total_biomass+=food[ii].GetE();};
				return total_biomass;

			}
		int food_Biomass(Food *food)
			{
				int ii=0;
				int total_biomass=0;
				for(ii=0;ii<Nfood;ii++)
					{total_biomass+=food[ii].GetE();};
				return total_biomass;

			}
		int Bichos_Biomass(Bichin *Bichos)
			{
				int ii=0;
				int total_biomass=0;
				for(ii=0;ii<Ni;ii++)
					{total_biomass+=Bichos[ii].GetE();};
				return total_biomass;
			}	
		void food_distribution(Food *food,int Nfood,double mux,double muy,double sigmax, double sigmay, double Rfood,Crandom &ran64, int dis )
			{
				int ix, iy;
				int  Ploted_energy=Biome_energy;
				if(dis==0)
					{
						for (int ii = 0; ii < Nfood; ii++)  //Se distribuyen N comidas
							{
								// Escogemos posición de la comida al azar segun una distribucion bigaussiana
								ix =  2*L*ran64.r() - L;
								iy =  2*L*ran64.r() - L;
								// Evitamos las fronteras del ambiente
								if (ix < -L)
									{ix = -L;}
								if (iy < -L)
									{iy = -L;}
								if (ix > L)
									{ix = L;}
								if (iy > L)
									{iy = L;}
								food[ii].Start(ix, iy, E_inicial, Rfood); //Se inicializa la comida con posiciones al azar, Energia E_inicial, radio R_food, porción val
								
								Ploted_energy-=E_inicial;
								if(Ploted_energy<=0)
									{break;};
							}
					}

				if(dis==1)
					{
						for (int ii = 0; ii < Nfood; ii++)  //Se distribuyen N comidas
							{
								// Escogemos posición de la comida al azar segun una distribucion bigaussiana
								ix = (int)ran64.gauss(mux, sigmax);
								iy = (int)ran64.gauss(muy, sigmay);
								// Evitamos las fronteras del ambiente
								if (ix < -L)
									ix = -L;
								if (iy < -L)
									iy = -L;
								if (ix > (L))
									ix = L;
								if (iy > L)
									iy = L;
								food[ii].Start(ix, iy, E_inicial, Rfood); //Se inicializa la comida con posiciones al azar, Energia E_inicial, radio R_food, porción val
								Ploted_energy-=E_inicial;
								if(Ploted_energy<=0)
									{break;};
							}
					}
				if(dis==2)
					{
						double place= ran64.r();

						for (int ii = 0; ii < Nfood; ii++)
							{
								place= ran64.r();
								if(place>0.5)
									{
										
									}
							}


					}
			}

		friend class Food;
};

void Selection::Birth(Bichin &BichoP, Bichin &BichoH, double t, int prob1, int prob2, Crandom &ran64)
{
	Nlive += 1; //Aumente el número de bichines vivos(que se dibujan)
	int E_original=BichoP.E;
	
	BichoP.E = int(BichoP.E / 2); //Divida la energía del padre en 2
	BichoH.Start(BichoP.x, BichoP.y, BichoP.E, 1, BichoP.R, ran64);//Inicialice al bichin hijo con la posición del padre, su energía (la mitad de la original), masa de 1, radio del padre y número aleatorio para determinar su genética(después se iguala a la del padre)

	//cout<<E_original<<" "<<BichoP.E+BichoH.E<<"-----\n";
	// Mutation
	for (int ii = 0; ii < P; ii++) //Iguale genética del Hijo a la del padre
	{
		BichoH.moves[ii] = BichoP.moves[ii];
	}
	BichoH.moves[prob1] -= 0.01;  //Disminuya el gen prob
	BichoH.moves[prob2] += 0.01;  //Aumente el gen prob2
	
}

void Selection::Uniform(Food *food, int N, double Rfood, Crandom &ran64)
{
	int ix, iy;
	int  Ploted_energy=Biome_energy;

	for (int ii = 0; ii < N; ii++)  //Se distribuyen N comidas
		{
			// Escogemos posición de la comida al azar segun una distribucion bigaussiana
			ix =  2*L*ran64.r() - L;
			iy =  2*L*ran64.r() - L;
			// Evitamos las fronteras del ambiente
			if (ix < -L)
				{ix = -L;}
			if (iy < -L)
				{iy = -L;}
			if (ix > L)
				{ix = L;}
			if (iy > L)
				{iy = L;}
			food[ii].Start(ix, iy, E_inicial, Rfood); //Se inicializa la comida con posiciones al azar, Energia E_inicial, radio R_food, porción val
			
			Ploted_energy-=E_inicial;
			if(Ploted_energy<=0)
				{break;};
		}
}

void Selection::Spread(Food *food, int N, double mu, double sigma, double Rfood, Crandom &ran64)
{
	int ix, iy;
	int Ploted_energy=Biome_energy;
	for (int ii = 0; ii < N; ii++)  //Se distribuyen N comidas
	{
		// Escogemos posición de la comida al azar segun una distribucion bigaussiana
		ix = (int)ran64.gauss(mu, sigma);
		iy = (int)ran64.gauss(mu, sigma);
		// Evitamos las fronteras del ambiente
		if (ix < -L)
			ix = -L;
		if (iy < -L)
			iy = -L;
		if (ix > (L))
			ix = L;
		if (iy > L)
			iy = L;
		food[ii].Start(ix, iy, E_inicial, Rfood); //Se inicializa la comida con posiciones al azar, Energia E_inicial, radio R_food, porción val
		Ploted_energy-=E_inicial;
		if(Ploted_energy<=0)
			{break;};
	}
}

void Selection::RechargeFood(Food *food, Crandom &ran64)
{
	int index=0;
	int Cn=0;
	while (Energy_bank>E_inicial)
		{	
			//cout<<"ciclo\n";
			index=int(Nfood*ran64.r());
			if(food[index].GetE()==0)
				{
					food[index].ReStart(E_inicial,ran64);
					Energy_bank-=E_inicial;
				}
			Cn+=1.;
			if(Cn*2>Nfood)
				{
					cout<<"Fallo en la distribucion de comida, re escribir Energia total \n";
					break;
					};
			
		}
};




//----------------- Funciones de Animacion ----------
void StartAnimacion(void)
{
	salida << "set terminal gif animate" << endl;
	salida << "set output '3.gif'" << endl;
	salida << "unset key" << endl;
	salida << "set xrange[" << -L << ":" << L << "]" << endl;
	salida << "set yrange[" << -L << ":" << L << "]" << endl;
	salida << "set size ratio -1" << endl;
	salida << "set parametric" << endl;
	salida << "set trange [0:7]" << endl;
	salida << "set isosamples 12" << endl;
}
void InicieCuadro(void)
{
	salida << "plot 0,0 ";
}
void TermineCuadro(void)
{
	salida << endl;
}
void StartBlender(int t)
{
	salida << t << "\t";
}
//-----------  Programa Principal --------------
int main(void)
	{	
		salida.open("console_out.gp");
		Bichin Bichitos[Ni]; 							//Array de bichines con numero maximo de bichines
		Food food[Nfood];  								//Array de food con numero maximo de food
		Selection Fate; 									//Nombre de clase Selection
		Crandom ran64(1);  								//Semilla del generador aleatorio
		double R = 5.0;  									//Radio del bichin
		int Ehijos = 20;   								// Min Energy for reproduction
		double Thijos = 20; 							// Min Time for reproduction
		double Rfood = 2;  								//Radio de la comida
		double prob;  										//variable auxiliar para mover bichines
		int prob1, prob2;  								//variables auxiliares para las mutaciones

		int total_bio=0,food_bio=0,Bichos_bio=0;

		int qq=0,nn=0;
		bool Blive=false;

		// Inicializar los bichines
		for (int jj = 0; jj < Nlive; jj++)
			{
				//Inicialice todos los bichines en el origen, 500 de energía, 1 de masa, radio R, ran64 para su genética

				//INICIALIZA EN POSICIONES ALEATORIAS
				double bix = 2*L*ran64.r() - L;
				double biy = 2*L*ran64.r() - L;
				Bichitos[jj].Start(bix, biy, E_inicial, 1, R, ran64);
			}


		if(food_dis==0)
			{Fate.Uniform(food, int(Nfood/2), Rfood, ran64);} //Distribuya Nfood(food maxima) con distribución uniforme}
		else if(food_dis==1)
			{Fate.Spread(food, int(Nfood/2),mu,sigma, Rfood, ran64);}
		else
			{
				cout<<"No se ha escogido una distribucion, linea 301"<<endl;
				return 0;
			};
		
		StartAnimacion(); // Dibujar


		for (int t = 0, tdibujo = 0; t < TMAX; t ++)
		{ 
			//total_bio=Fate.Biomass(food,Bichitos);
			//food_bio=Fate.food_Biomass(food);
			//Bichos_bio=Fate.Bichos_Biomass(Bichitos);
			//cout<<food_bio<<" "<<Bichos_bio<<" "<<Energy_bank<<"\n";
			//cout<<total_bio+Energy_bank<<"\n";
			
			for (int ii = 0,nn=0; ii < Ni; ii++)  						//Para todos los bichines vivos
				{	
					if (Bichitos[ii].Alive())  										//Si el bichin esta vivo
					{ nn++;
						
																											
						prob = ran64.r();  //Genere un número aleatorio
						Bichitos[ii].Move(K, prob); //Muevase con ese numero
				

						if (Bichitos[ii].GetE() > Ehijos && Bichitos[ii].GetT() > Thijos &&  int(Bichitos[ii].GetE())%2==0)  //Si el bichin cumple las 2 condiciones para reproducirse
							{
								prob1 = int(P * ran64.r());
								prob2 = int(P * ran64.r());

								Blive=false;
								for(qq=0;qq<Ni;qq++)
									{
										Blive=Bichitos[qq].Alive();
										if(!Blive){break;};
									}
									if(Blive){cout<<"Poblacion maxima \n";}

								Fate.Birth(Bichitos[ii], Bichitos[qq], t, prob1, prob2, ran64);   //Escoga a un nuevo bichin del array como hijo
							}


						for (int jj = 0; jj < Nfood; jj++)
							{   //Para toda la comida, revise si puede alimentarse con ella
								food[jj].Feed(Bichitos[ii]);
							}
						
						if(Bichitos[ii].GetE()<0)
							{Nlive-=1;}
						
						if(nn==Nlive)
							{break;};
					}
				}

			InicieCuadro(); //Dibuje los bichines vivos y la comida viva

			for (int ii = 0; ii < Nfood; ii++)
				{	
					if (food[ii].GetE() > 0)
					{	
						food[ii].Print();
					}
				}
			for (int ii = 0; ii < Ni; ii++)
				{
					if (Bichitos[ii].Alive())
					{
						Bichitos[ii].Print();
					}
				}
			//cout<<food[12].GetE()<<"-----\n";
			TermineCuadro();
			Fate.RechargeFood(food, ran64);
		}
		salida.close();
		return 0;
	}
