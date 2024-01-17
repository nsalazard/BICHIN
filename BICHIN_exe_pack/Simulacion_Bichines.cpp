#define _USE_MATH_DEFINES 
#include <iostream>
#include <cmath>
#include "Random64.h"
#include <fstream>
#include <exception>
#include <memory>
#include <unistd.h>
#include <string>

ofstream video_gp;
ofstream grafica;
ofstream Nodes;
ofstream Edges;
ofstream Hald;
ofstream gene_data;

//---------- Constantes --------
const int P = 8;            			// Numero de parámetros de los bichines
const int L =700;           			// Espacio 2L*2L
const double K = 10;        			// Distancia recorrida en cada mov. por el bichin
const double E_gordo = 50;				// Energía a partir de la cual bichin no puede comer
const int Ni = 4000; 					// Numero maximo de bichines (?)
const int Nfood = 10000;  				// Numero de maximo comida | Nunca debe ser alcanzado.
int E_inicial = 10;  					// Energía inicial de la comida
int Nlive = 100;  						// Numero inicial de bichines
int Energy_bank = 0; 					// El banco temporal de energia
int Biome_energy=60000;					// Evita un bug con el colocamiento de la comida
										// Como buena practica  Biome_energy<Nfood*E_inicial; 
										// Podria funcionar incluso si esta condicion no se cumple pero se corre un riesgo.
int food_dis=0;	
const double mu = 0.0, sigma = L / 4;  	// Parámetros distribución gaussiana de comida
const int genetic_start = 0; 			// Permite que todos los genes inicien con las mismas probabilidades 



int TMAX=1000;
int rand_seed=1;
int data_tick=100;
bool video=false;

void handleFlags(int argc, char *argv[]) {
	if(argc == 1) 
		{
			std::cout << "No flags were given. Using default values." << std::endl;
			return;
		}
	int opt;

	while((opt = getopt(argc, argv, "t:d:s:v:i:")) != -1)  
	{  
		switch(opt)  
		{  

			case 't': 
				TMAX = std::stoi(optarg); 
				break;  
			case 'd':  
				food_dis = std::stoi(optarg);
				break;  
			case 's':  
				rand_seed = std::stoi(optarg);
				break;
			case 'i':
				data_tick=std::stoi(optarg);
				break;
			case 'v':
				video=true;
				break;  
			
			case '?':  
				std::cout << "unknown option: " << optopt << std::endl;  
				break;  
		}  
	}  
}
void displayProgressBar(float progress, int barWidth = 50, std::string line="=", float total=1.0)
	{
		std::cout << "[";
		int pos = static_cast<int>(progress * barWidth);
		for (int i = 0; i < barWidth; ++i) 
			{
				if (i < pos)
					{std::cout << line;}
				else if (i == pos)
					{std::cout << ">";}
				else
					{std::cout << " ";}
			}
		std::cout << "] " << static_cast<int>(progress * 100.0)<< "%\r  "+std::to_string(total);
		std::cout.flush();
	}	
void write_config(std::string name)
	{
		std::ofstream config;
		config.open(name);
		config <<"TMAX,data_tick,food_dis,Biome_energy,E_incial,E_gordo\n";
		config <<TMAX<<","<<data_tick<<","<<food_dis<<","<<Biome_energy<<","<<E_inicial<<","<<E_gordo<<"\n";
		config.close();
	}




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
					if(genetic_start == 0){
						for (int ii = 0; ii < P; ii++)
							{moves[ii] = 0.125;}
					}

					else if(genetic_start ==1){
						moves[0]=1;
						for(int ii=1;ii<P;ii++)
							{moves[ii]=0;}

					}
					else//LLena moves[] con porcentajes aleatorios tal que su suma de 1	
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
								{moves[ii] = moves[ii] / sum;}
						}
					theta=M_PI/2*int(P*ran64.r());
				}
			void Print(void)
				{video_gp << " , " << x << "+" << R << "*cos(t)," << y << "+" << R << "*sin(t)" << " lt -1";}

			bool Alive(void) { return E > 0; }; //determina si vive
			int Main_gene(void)
				{
					int ii=0,imax=0;
					double max_gene=-1,gene=0;
					for(ii=0;ii<P;ii++)
						{
							gene=moves[ii];
							if(max_gene<gene )
								{
									imax=ii;
									max_gene=gene;
								}

						}
					return imax;
				}
			double Gene_value(int gene){return moves[gene];}
			friend class Selection;
			friend class Food;
	};


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
		void ReStart(int E0,Crandom &ran64,double mux,double muy,double sigmax,double sigmay)
			{	

				if(food_dis==0) //uniform
					{
						x =  2*L*ran64.r() - L;
						y =  2*L*ran64.r() - L;
					}
				if(food_dis==1)
					{
						x = (int)ran64.gauss(mu, sigma);
						y = (int)ran64.gauss(mu, sigma);		
					}
				if(food_dis==2)
					{
						double place= ran64.r();

								place= ran64.r();
								if(place>0.5)
									{
										x =  2*L*ran64.r() - L;
										y =  2*L*ran64.r() - L;
									}
								else	
									{	
										x=sigmax*(ran64.r()-0.5)+mux;
										y=sigmay*(ran64.r()-0.5)+muy;
									}
								if (x < -L)
									{x = -L;}
								if (y < -L)
									{y = -L;}
								if (x > L)
									{x = L;}
								if (y > L)
									{y = L;}
					}
				
				E=E0;
			}
		void Feed(Bichin &Bicho)  //Alimente al bicho  
			{
				double dis = sqrt(pow(Bicho.x - x, 2.0) + pow(Bicho.y - y, 2.0));
				if (dis < (R + Bicho.R) && E > 0.0 && Bicho.E < E_gordo)//Si el bicho y la comida se superponen y la energía de la comida es mayor a cero y la energía del bichin no es mayor a E_gordo
					{ 
						Bicho.E += E; //Aumente la energía del bicho en val 
						E = 0; //Disminuya la energía del bicho en val
						// video_gp << Bicho.E << "\n";
						//cout<<"Comido \n";
					}
			}
		void Print(void){video_gp << " , " << x << "+" << R << "*cos(t)," << y << "+" << R << "*sin(t)" << " lt 7";}
		void Blender(void);

		friend class Selection;
		friend class Bichin;
};

//--- clase Selection----
class Selection
{
	private:
	public:
		void Birth(Bichin &BichoP, Bichin &BichoH, double t,Crandom &ran64)  //Reproduce un bicho Padre en 2 bichos Hijos teniend en cuenta el tiempo de vida del padre. prob determina el gen que disminuye, prob 2 el que aumenta
			{
				Nlive += 1; //Aumente el número de bichines vivos(que se dibujan)
				int E_original=BichoP.E;

				int prob1= int(P * ran64.r());
				int prob2= int(P * ran64.r());
				
				BichoP.E = int(BichoP.E / 2); //Divida la energía del padre en 2
				BichoH.Start(BichoP.x, BichoP.y, BichoP.E, 1, BichoP.R, ran64);//Inicialice al bichin hijo con la posición del padre, su energía (la mitad de la original), masa de 1, radio del padre y número aleatorio para determinar su genética(después se iguala a la del padre)

				// Mutation
				for (int ii = 0; ii < P; ii++) //Iguale genética del Hijo a la del padre
					{BichoH.moves[ii] = BichoP.moves[ii];}

				if(BichoH.moves[prob1] > 0.01)//Disminuya el gen prob
					{BichoH.moves[prob1] -= 0.01;}

				else
					{
						for(int ii=0;ii<P; ii++)
							{
								prob1= int(P * ran64.r());
								if(BichoH.moves[prob1] > 0.01)
									{BichoH.moves[prob1] -= 0.01;break;}
							}	
					}

				BichoH.moves[prob2] += 0.01;  //Aumente el gen prob2
			}
		void Gaussian(Food *food, int N, double mu, double sigma, double Rfood, Crandom &ran64)  //Distribuya la comida con una gausiana
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
		void Uniform(Food *food, int N, double Rfood, Crandom &ran64)  //Distribuya la comida uniformemente
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
		void RechargeFood(Food *food, Crandom &ran64) //Recargue la comida de manera aleatoria en el ambiente
			{
				int index=0;
				int Cn=0;
				while (Energy_bank>E_inicial)
					{	
						//cout<<"ciclo\n";
						index=int(Nfood*ran64.r());
						if(food[index].GetE()==0)
							{
								food[index].ReStart(E_inicial,ran64,L/4,-L/4,sigma*2,sigma*2);
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
										ix =  2*L*ran64.r() - L;
										iy =  2*L*ran64.r() - L;
									}
								else	
									{	
										ix=sigmax*(ran64.r()-0.5)+mux;
										iy=sigmay*(ran64.r()-0.5)+muy;
									}
								if (ix < -L)
									{ix = -L;}
								if (iy < -L)
									{iy = -L;}
								if (ix > L)
									{ix = L;}
								if (iy > L)
									{iy = L;}
								food[ii].Start(ix, iy, E_inicial, Rfood);
								Ploted_energy-=E_inicial;
								if(Ploted_energy<=0)
									{break;};
							}


					}
			}
    	double Dist_Taxi(Bichin &Bicho1, Bichin &Bicho2)
			{
				double sum=0.0;
				for(int ii=0; ii<P;ii++) 
					{sum+= abs(Bicho1.moves[ii]-Bicho2.moves[ii]);}
				return sum;
    		}
		void Genetic_Nodes(Bichin *Bichos)
			{
				Nodes<<"id,gen_predominante\n";
				int ii=0,cc=0;
				
				string preferencia;
				for(ii=0;ii<Ni;ii++)
					{
						if(Bichos[ii].Alive())
							{	cc=Bichos[ii].Main_gene();
								preferencia=Direction(cc);
								Nodes<<"B"<<ii<<","<<preferencia<<"\n";
							}

					}
			}	
		void Genetic_Edges(Bichin *Bichos, double threshold)
			{
				int ii=0,jj=0,kk=0;
				double distancia=0;
				double peso=0;
				Edges<<"from,to,weight\n";

				for(ii=0;ii<Ni;ii++)
					{for(jj=ii+1;jj<Ni;jj++)
						{
							if(Bichos[ii].Alive() && Bichos[jj].Alive())
								{
									distancia=Dist_Taxi(Bichos[ii],Bichos[jj]);
									peso=1/distancia;
									if(distancia==0)
										{peso=0;}

									if(peso>threshold)
										{Edges<<"B"<<ii<<","<<"B"<<jj<<","<<peso<<"\n";}	
								}
						}
					}
			}
		double Prom_gene(Bichin *Bichos, int gene)
			{
				int oo=0;
				double prom=0;
				for(oo=0;oo<Nlive;oo++)
					{prom += Bichos[oo].Gene_value(gene);}

				prom=prom/Nlive;
				return prom;
				
			}
		double Std_gene(Bichin *Bichos, int gene)
			{
				double sd=0;// desviacion std
				double prom=0;
				double sum=0;
				double val=0;
				int oo=0;	//contador

				prom=Prom_gene(Bichos,gene);

				for(oo=0;oo<Nlive;oo++)
					{	
						val=Bichos[oo].Gene_value(gene)-prom;
						sd+= val*val;

					}
				sd=sqrt(sd/Nlive);
				return sd;

			}
		void Haldane(int time,Bichin *Bichos)
			{	
				// H es una unidad sin dimensiones para cuantificar tasas evolutivas, https://earth.geology.yale.edu/~ajs/1993/11.1993.17Gingerich.pdf
				int gene=0;
				Hald<<time<<" "<<Nlive<<" ";
				for(gene=0;gene<P;gene++)
					{Hald<<Prom_gene(Bichos,gene)<<" "<<Std_gene(Bichos,gene)<<" ";}

				Hald<<"\n";	
			}
		void Genetic_out(int time,Bichin *Bichos)
			{	
				string name1;
				name1="Genes_datos/genes_"+to_string(time)+".csv";
				gene_data.open(name1);
				//gene_data<<"name,PC1,PC2,PC\n";
				for(int ii=0;ii<Ni;ii++)
					{
						if(Bichos[ii].Alive())
							{	
								
								for(int jj=0;jj<P;jj++)
									{
										gene_data<<Bichos[ii].Gene_value(jj);
										if((jj<P-1))
											{gene_data<<",";}
									}
								gene_data<<"\n";
							}
					}
				gene_data.close();

				
			}

		string Direction(int dir)
			{	
				string pref;
				if(dir==0)
					{pref="F";}
				else if(dir==1)
					{pref="L1";}
				else if(dir==2)
					{pref="L2";}
				else if(dir==3)
					{pref="L3";}
				else if(dir==4)
					{pref="B";}
				else if(dir==5)
					{pref="R3";}
				else if(dir==6)
					{pref="R2";}
				else if(dir==7)
					{pref="R1";}
				else
					{pref="Error";
					cout<<dir<<"Error en la asignacion de direccion preferencial\n";}
				return pref;
			}		
			
		friend class Food;
};

//----------------- Funciones de Animacion ----------
void StartAnimacion(void)
	{
		video_gp << "set terminal gif animate" << endl;
		video_gp << "set output 'Bichin.gif'" << endl;
		video_gp << "unset key" << endl;
		video_gp << "set xrange[" << -L << ":" << L << "]" << endl;
		video_gp << "set yrange[" << -L << ":" << L << "]" << endl;
		video_gp << "set size ratio -1" << endl;
		video_gp << "set parametric" << endl;
		video_gp << "set trange [0:7]" << endl;
		video_gp << "set isosamples 12" << endl;
	}
void InicieCuadro(void)
	{
		video_gp << "plot 0,0 ";
	}
void TermineCuadro(void)
	{
		video_gp << endl;
	}

//-----------  Programa Principal --------------
int main(int argc, char **argv)
	{	

		handleFlags(argc, argv);

		write_config("config.csv");


		video_gp.open("console_out.gp");
		grafica.open("population.txt");
		Hald.open("Haldanes.txt");

		
		Bichin Bichitos[Ni]; 							//Array de bichines con numero maximo de bichines
		Food food[Nfood];  								//Array de food con numero maximo de food
		Selection Fate; 								//Nombre de clase Selection
		Crandom ran64(rand_seed);  						//Semilla del generador aleatorio
		double R = 5.0;  								//Radio del bichin
		int Ehijos = 20;   								//Min Energy for reproduction
		double Thijos = 40; 							//Min Time for reproduction
		double Rfood = 2;  								//Radio de la comida

		double prob;  									//variable auxiliar para mover bichines
		int prob1, prob2;  								//variables auxiliares para las mutaciones
		int gene=0;
    
		int total_bio=0,food_bio=0,Bichos_bio=0;

		int qq=0,nn=0,live_counter;
		bool Blive=false;
		string name,name2;

		for (int jj = 0; jj < Nlive; jj++)//Inicialice todos los bichines, ran64 para su genética
			{
				double bix = 2*L*ran64.r() - L;
				double biy = 2*L*ran64.r() - L;
				Bichitos[jj].Start(bix, biy, E_inicial, 1, R, ran64);//INICIALIZA EN POSICIONES ALEATORIAS
			}

		Fate.food_distribution(food,int(Nfood/2),L/4,-L/4,sigma*2,sigma*2,Rfood,ran64,food_dis);
		
		if(video)StartAnimacion(); // Dibujar

		for (int t = 0, tdibujo = 0; t < TMAX; t ++)
		{ 	

			
			//total_bio=Fate.Biomass(food,Bichitos);
			//food_bio=Fate.food_Biomass(food);
			//Bichos_bio=Fate.Bichos_Biomass(Bichitos);
			//cout<<food_bio<<" "<<Bichos_bio<<" "<<Energy_bank<<"\n";
			//cout<<total_bio+Energy_bank<<"\n";

			if(t%data_tick==0)
				{	
					//name="Nodos/Nodos"+to_string(t)+".csv";
					//Nodes.open(name);
					//Fate.Genetic_Nodes(Bichitos);
					//Nodes.close();
					
					//name2="Edges/Edges"+to_string(t)+".csv";
					//Edges.open(name2);
					//Fate.Genetic_Edges(Bichitos, 7);
					//Edges.close();
					
					//cout<<"std dev gen "<<gene<<": "<<Fate.Std_gene(Bichitos,gene)<<"\n";

					displayProgressBar(float(t)/TMAX,70,"=",t);
					Fate.Haldane(t,Bichitos);
					Fate.Genetic_out(t,Bichitos);
				}
			
			for (int ii = 0,nn=0; ii < Ni; ii++)  									//Para todos los bichines vivos
				{	
					if (Bichitos[ii].Alive())  										//Si el bichin esta vivo
					{ nn++;
						
																											
						prob = ran64.r();  //Genere un número aleatorio
						Bichitos[ii].Move(K, prob); //Muevase con ese numero
				

						if (Bichitos[ii].GetE() > Ehijos && Bichitos[ii].GetT() > Thijos &&  int(Bichitos[ii].GetE())%2==0)  //Si el bichin cumple las 2 condiciones para reproducirse
							{
								Blive=false;
								for(qq=0;qq<Ni;qq++)
									{
										Blive=Bichitos[qq].Alive();
										if(!Blive){break;};
									}
									if(Blive){cout<<"Poblacion maxima \n";}

								Fate.Birth(Bichitos[ii], Bichitos[qq], t, ran64);   //Escoga a un nuevo bichin del array como hijo
							}


						for (int jj = 0; jj < Nfood; jj++)//Para toda la comida, revise si puede alimentarse con ella
							{food[jj].Feed(Bichitos[ii]);}

						// if(Bichitos[ii].GetE()<10)
						// 	{cout<<Bichitos[ii].GetE()<<"----------"<<"\n";}

						if(Bichitos[ii].GetE()==0)//cout<<"muerte-----------------------------\n";
							{Nlive-=1;}
						
						if(nn==Nlive)
							{break;};
					}
				}

			if(video)InicieCuadro(); //Dibuje los bichines vivos y la comida viva

			for (int ii = 0; ii < Nfood; ii++)
				{ if (food[ii].GetE() > 0)
						{food[ii].Print();}
				}
			
			
			for (int ii = 0; ii < Ni; ii++)
				{
					if (Bichitos[ii].Alive())
						{	
							live_counter+=1;
							if(video)
								{Bichitos[ii].Print();}
							
						}
				}

			grafica<<t<<" "<<Nlive<<"\n";

			if(video)TermineCuadro();
			Fate.RechargeFood(food, ran64);

		}
		video_gp.close();
		grafica.close();
		Hald.close();
		return 0;
	}
