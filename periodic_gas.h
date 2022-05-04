#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <time.h>

using namespace std;

//Funzione segno
double sign (double a)
{
double b=0;

	if (a > 0) {b = b + 1.;}
	else if (a < 0) {b = b - 1.;};

return b;
};

//Funzione per il modulo di un vettore
double modulov (vector<double> vij)
{
	double d;
	d = sqrt(vij.at(0)*vij.at(0) + vij.at(1)*vij.at(1));

	return d; 
};



//Classe che contiene tutte le variabili utili di un atomo o una particella. Si può espandere con ad esempio un parametro che chiarisca se la particella è un protone o un neutrone ecc..
class atomo

	{
	double P_x; 	//posizione x
	double P_y;		//posizione y
	double V_x;		//velocita' x
	double V_y;		//velocita' y
	double F_x;		//forza sentita lungo x
	double F_y;		//forza sentita lungo y
	double Up;		//en potenziale atomo
	double M;		//massa

	public:

	atomo (double a, double b, double c, double d, double e, double f, double g, double h)
		{
		P_x = a;
		P_y = b;
		V_x = c;
		V_y = d;
		F_x = e;		
		F_y = f;
		Up = g;	
		M = h;
		}

//costruttore delegato utile per quando serve trovare il vettore rij tra due atomi i e j tenendo conto delle condizioni periodiche
	atomo(double a, double b) :  atomo(a,b,0,0,0,0,0,0) {
								P_x = a;
								P_y = b;
								}

	double x () {return P_x;}			void setx (double a) {P_x = a;}

	double y () {return P_y;}			void sety (double a) {P_y = a;}

	double vx () {return V_x;}			void setvx (double a) {V_x = a;}

	double vy () {return V_y;}			void setvy (double a) {V_y = a;}

	double fx () {return F_x;}			void setfx (double a) {F_x = a;}

	double fy () {return F_y;}			void setfy (double a) {F_y = a;}

	double U () {return Up;}			void setu (double a) {Up = a;}

	double EK () {
					double Ek = 0.5*M*(V_x*V_x + V_y*V_y) ;
					return Ek;
				  }


	};



//distanza tra atomi, utile solo nella prossima funzione
double dist (atomo i, atomo j)
{
	double d;
	d = sqrt(pow(i.x()-j.x() ,2.) + pow(i.y()-j.y() ,2.));

	return d; 
};



//Funzione fondamentale che dati due atomi i e j trova, tra j e le sue 8 repliche dovute alle condizioni periodiche, quello più vicino a i e anche il vettore congiungente rij. Note le due componenti del GIUSTO vettore congiungente si potrà trovare l'interazione tra i due atomi
// Vedi dispense pag. 22
// Vedi "periodic boundary conditions"
vector<double> vmin (atomo i, atomo j0)
{
//j0 e' l'atomo senza effetto pacman, j1...j8 quelli dovuti all'effetto pacman
vector <double> V;
vector <atomo> J;
double* D = new double[9];
int min=0;
											J.push_back(j0);
atomo j1 (j0.x() , j0.y()+1.);				J.push_back(j1);
atomo j2 (j0.x()+1. , j0.y()+1.);			J.push_back(j2);
atomo j3 (j0.x()+1. , j0.y());				J.push_back(j3);
atomo j4 (j0.x()+1. , j0.y()-1.);			J.push_back(j4);
atomo j5 (j0.x() , j0.y()-1.);				J.push_back(j5);
atomo j6 (j0.x()-1. , j0.y()-1.);			J.push_back(j6);
atomo j7 (j0.x()-1. , j0.y());				J.push_back(j7);
atomo j8 (j0.x()-1. , j0.y()+1.);			J.push_back(j8);

D[0] = (dist(i,j0));		
D[1] = (dist(i,j1));			
D[2] = (dist(i,j2));	
D[3] = (dist(i,j3));	
D[4] = (dist(i,j4));	
D[5] = (dist(i,j5));	
D[6] = (dist(i,j6));	
D[7] = (dist(i,j7));	
D[8] = (dist(i,j8));	

	for(int j=1; j<9; j++)
  	{
		if(D[j] < D[min]) 
   		min = j;
	}; 

atomo j = J.at(min);

//Vettore rij
V.push_back(i.x()-j.x());
V.push_back(i.y()-j.y());

return V;
};



//Il cuore del propgramma, questa è la funzione che applica una volta il propagatore Velocity-Verlet. Per passare dallo step temporale i a quello i+1, questo propagatore evolve prima le posizioni richiedendo di sapere le velocità al tempo i e le forza al tempo i. Fatto ciò serve fare la cosa scomoda, ovvero calcolare le forze al tempo i+1, che saranno diverse visto che le posizioni sono cambiate, e solo noto ciò si calcolano le velocità.
//Per maggiori dettagli vedi pag 22 delle dispense.
vector <atomo> Iter (double sigma, double epsilon, double m, double dt, vector<atomo> V0)
{
					//vector<atomo> V,Vf;
	
	//EVOLUZIONE TEMPORALE POSIZIONI

	for(int i=0;i<V0.size();i++)
	{
	atomo Xi = V0.at(i);

	double px = Xi.x() + Xi.vx()*dt + Xi.fx()*dt*dt/(2.*m); //(87)
	double py = Xi.y() + Xi.vy()*dt + Xi.fy()*dt*dt/(2.*m); //(87)

//Effetto pacman
	if		(px>1.) {px=px-1.;}
	else if (px<0.)	{px=px+1.;}
	else {px=px;};

	if		(py>1.)	{py=py-1.;}
	else if (py<0.)	{py=py+1.;}
	else {py=py;};

//Atomi con SOLO posizione aggiornata

	Xi.setx(px);
	Xi.sety(py);

	V0.at(i) = Xi;
 					//atomo X (px,py,Xi.vx(),Xi.vy(),Xi.fx(),Xi.fy(),Xi.U(),m);
					//V.push_back(X);
	};


	//EVOLUZIONE TEMPORALE VELOCITA'

//serve ricalcolare la forza! L'energia potenziale non serve per l'evoluzione temporale, serve per controllare che l'energia del sistema si conservi
	for(int i=0;i<V0.size();i++)
	{
	//atomo Xi = V.at(i);
	atomo Xi = V0.at(i);

	double ffx=0;
	double ffy=0;
	double Uf=0;

	for(int j=0;j<V0.size();j++)
		{
			//atomo Xj = V.at(j);
			atomo Xj = V0.at(j);

			vector<double> rij = vmin(Xi,Xj);
			double dij = modulov(rij);

//Considero solo le particelle entro 3 sigma da quella in esame. Per ditanze maggiori l'inetrazione è moolto debole, così i conti si alleggeriscono non poco. La distanza dev'essere > 0 altrimenti le particelle interagirebbero con loro stesse
			if(dij <= 3.*sigma && dij > sigma/100.)
			{
			double rx = rij.at(0);
			double ry = rij.at(1);

//sommo tutti i contruibuti di forza e potenziale di lennard-jones. Si veda pagina 21 delle dispense anche per il significato fisico si sigma ed epsilon
			ffx = ffx + 24.0*epsilon*(2.*pow(sigma/dij , 12.) - pow(sigma/dij , 6.))*(rx/pow(dij , 2.)); //(86)
			ffy = ffy + 24.0*epsilon*(2.*pow(sigma/dij , 12.) - pow(sigma/dij , 6.))*(ry/pow(dij , 2.)); //(86)
			Uf = Uf + 4.*epsilon*(pow(sigma/dij , 12.) - pow(sigma/dij , 6.)); //(84)
	
			}
			else
			{double null;};
		};

	double wx = Xi.vx() + dt*(ffx+Xi.fx())/(2.*m); //(90)
	double wy = Xi.vy() + dt*(ffy+Xi.fy())/(2.*m); //(90)

//Ora finalmente posso aggiornare tutte le variabili dell'atomo i-esimo. Le posizioni no perchè l'avevo fatto prima

 	//atomo XX (Xi.x(),Xi.y(),wx,wy,ffx,ffy,Uf,m);
	//Vf.push_back(XX);
	
	Xi.setvx(wx);
	Xi.setvy(wy);
	Xi.setfx(ffx);
	Xi.setfy(ffy);
	Xi.setu(Uf);

	V0.at(i) = Xi;
	
	};

//V0.clear();
//return Vf;

return V0;

};


//Funzione che dato un file di double conta quanti stanno nelle varie suddivisioni cosi' da fare i file di testo per gli istogrammi, in questo caso delle velocita'
vector<int> Histo_double(vector<double> Dati, double a, double b, int N)
	{
		double h=(b-a)/(N);
		vector<int> V(N);

	for(int i=0; i<Dati.size(); i++)
		{
		for(int l=0; l<N; l++)
			{
			if(Dati.at(i) > a+l*h  &&  Dati.at(i) < a+(l+1.)*h )
				{
				V.at(l) = V.at(l) + 1;
				}
				else { V.at(l) = V.at(l); };
			};
		};

	return V;
	};

//Classe LCG (linear congruential generator) che genera numeri casuali tra 0 e 1 (da cui poi si possono generare in modo da far seguire una qualsiasi distribuzione). La classe fa anche altre cose inuili, quella davvero importante è il Single_Number, infatti creo sempre le LCG con il costruttore delegato e le uso solo per quella funzione.
//Per maggiori dettagli vedere pagina 39 delle dispense (non lo shuffling, non lo uso)
class LCG								

{

vector <double> r;
int NNN;
const long int PARa = pow(7,5);
const long int PARm = pow(2,31)-1;
const int PARc=0;
long int xx0=0,xx1,xx2; 

public:

LCG() {}

double Single_Number () 
	{
	
	if(xx0==0)
	{
	xx0 = (PARa*clock()*clock()+PARc)%PARm;
	xx1 = (PARa*xx0+PARc)%PARm;
	}
	else
	{xx1 = xx2;};

	xx2 = (PARa*xx1+PARc)%PARm;
	return ((double) xx2/(double) PARm);
	}

};
