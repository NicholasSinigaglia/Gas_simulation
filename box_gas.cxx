#include <stdlib.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <stdio.h>
#include <time.h>

using namespace std;

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
	double U_c;		//en potenziale atomo interazione coulombiana
	double U_g;		//en potenziale gravitazionale
	double M;		//massa

	public:

	atomo (double a, double b, double c, double d, double e, double f, double g, double h, double i)
		{
		P_x = a;
		P_y = b;
		V_x = c;
		V_y = d;
		F_x = e;		
		F_y = f;
		U_c = g;		
		U_g = h;		
		M = i;
		}

//costruttore delegato utile per quando serve trovare il vettore rij tra due atomi i e j tenendo conto delle condizioni periodiche
	atomo(double a, double b) :  atomo(a,b,0,0,0,0,0,0,0) {
								P_x = a;
								P_y = b;
								}

	double x () {return P_x;}

	double y () {return P_y;}

	double vx () {return V_x;}

	double vy () {return V_y;}

	double fx () {return F_x;}

	double fy () {return F_y;}

	double EK () {
					double Ek = 0.5*M*(V_x*V_x + V_y*V_y) ;
					return Ek;
				  }

	double Uc () {return U_c;}
	double Ug () {return U_g;}

	};



//distanza tra atomi, utile solo nella prossima funzione
double dist (atomo i, atomo j)
{
	double d;
	d = sqrt(pow(i.x()-j.x() ,2.) + pow(i.y()-j.y() ,2.));

	return d; 
};



vector<double> Rij (atomo i, atomo j)
{

vector <double> V;

V.push_back(i.x()-j.x());
V.push_back(i.y()-j.y());

return V;
};



//Il cuore del propgramma, questa è la funzione che applica una volta il propagatore Velocity-Verlet. Per passare dallo step temporale i a quello i+1, questo propagatore evolve prima le posizioni richiedendo di sapere le velocità al tempo i e le forza al tempo i. Fatto ciò serve fare la cosa scomoda, ovvero calcolare le forze al tempo i+1, che saranno diverse visto che le posizioni sono cambiate, e solo noto ciò si calcolano le velocità.
//Per maggiori dettagli vedi pag 22 delle dispense.
vector <atomo> Iter (double sigma, double epsilon, double m, double dt, vector<atomo> V0, double g)
{
	vector<atomo> V,Vf;
	
	//EVOLUZIONE TEMPORALE POSIZIONI

	for(int i=0;i<V0.size();i++)
	{
	atomo Xi = V0.at(i);

	double px = Xi.x() + Xi.vx()*dt + Xi.fx()*dt*dt/(2.*m); //(87)
	double py = Xi.y() + Xi.vy()*dt + Xi.fy()*dt*dt/(2.*m); //(87)

//Effetto scatola
	double VX = Xi.vx();
	double VY = Xi.vy();

	if	(px>1.) {px = 2.-px;	VX = -Xi.vx();	VY = Xi.vy();}	else {double null;};
	if	(px<0.)	{px = -px;		VX = -Xi.vx();	VY = Xi.vy();}	else {double null;};
	if	(py>1.)	{py = 2.-py;	VX = Xi.vx();	VY = -Xi.vy();}	else {double null;};
	if	(py<0.)	{py = -py;		VX = Xi.vx();	VY = -Xi.vy();}	else {double null;};

//Atomi con SOLO posizione aggiornata
 	atomo X (px,py,VX,VY,Xi.fx(),Xi.fy(),Xi.Uc(),Xi.Ug(),m);

	V.push_back(X);
	};


	//EVOLUZIONE TEMPORALE VELOCITA'

//serve ricalcolare la forza! L'energia potenziale non serve per l'evoluzione temporale, serve per controllare che l'energia del sistema si conservi
	for(int i=0;i<V.size();i++)
	{
	atomo Xi = V.at(i);

	double ffx=0;				
	double ffy = -m*g;
	double UC = 0;
	double UG = m*g*Xi.y();


	for(int j=0;j<V.size();j++)
		{
			atomo Xj = V.at(j);

			vector<double> rij = Rij(Xi,Xj);
			double dij = modulov(rij);

//Considero solo le particelle entro 3 sigma da quella in esame. Per ditanze maggiori l'inetrazione è moolto debole, così i conti si alleggeriscono non poco. La distanza dev'essere > 0 altrimenti le particelle interagirebbero con loro stesse
			if(dij <= 3.*sigma && dij > sigma/100.)
			{
			double rx = rij.at(0);
			double ry = rij.at(1);

//sommo tutti i contruibuti di forza e potenziale di lennard-jones. Si veda pagina 21 delle dispense anche per il significato fisico si sigma ed epsilon
			ffx = ffx + 24.0*epsilon*(2.*pow(sigma/dij , 12.) - pow(sigma/dij , 6.))*(rx/pow(dij , 2.)); //(86)
			ffy = ffy + 24.0*epsilon*(2.*pow(sigma/dij , 12.) - pow(sigma/dij , 6.))*(ry/pow(dij , 2.)); //(86)
			UC = UC + 4.*epsilon*(pow(sigma/dij , 12.) - pow(sigma/dij , 6.)); //(84)
	
			}
			else
			{double null;};
		};
 
	double wx = Xi.vx() + dt*(ffx+Xi.fx())/(2.*m); //(90)
	double wy = Xi.vy() + dt*(ffy+Xi.fy())/(2.*m); //(90)

//Ora finalmente posso aggiornare tutte le variabili dell'atomo i-esimo. Le posizioni no perchè l'avevo fatto prima
 	atomo XX (Xi.x(),Xi.y(),wx,wy,ffx,ffy,UC,UG,m);

	Vf.push_back(XX);
	};

V0.clear();

return Vf;

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

LCG (int n) 
	{
	NNN = n;
	long int x0,x1;

	x0=clock();

	x1 = (PARa*x0+PARc)%PARm;
	for (int i=0;i<NNN;i++)
	{
	x0 = x1;
	x1 = (PARa*x0+PARc)%PARm;
	r.push_back(((double) x1/(double) PARm));	
	};	
	}

LCG() : LCG(0) {}

vector <double> Get_Vect () 
	{
	return r;
	}


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


void File_Numbers (string nome) 
	{
	ofstream fout;
	fout.open(nome);

	for(int i=0;i<NNN;i++)
	{
	fout << r.at(i) << endl;
	};

	//PER QUALCHE MOTIVO IGNOTO RIEMPIENDO UN ARRAY CON IL COSTRUTTORE IN OUTPUT QUI CI SONO SOLO ZERI
	fout.close();
	}


};


int main()

{	
	int N, n=1, nn=1, Nt, Ns=0, G; //n ed nn servono per la configurazione iniziale delle posizioni
	double sigma,epsilon,dt,m,g; 	//costanti
	vector<atomo> V0;			//conterrà configurazione iniziale
	LCG B;						//vedi header file
	int TEST,brow;

	cout << "inserire valore dell'accelerazione di gravita': " << endl;
	cin >> g;
	cout << "far partire tutte le particelle 'dal basso'? (1 si', altro no)" << endl;
	cin >> G;

	cout << "Si vuole tenere d'occhio la distribuzione del modulo delle velocita' tramite istogrammi? (1 si', altro no)" << endl;
	cin >> TEST;

	if(TEST==1)
	{
	cout << "Inserire numero di step in cui dividere l'intervallo [0,1.2], nel quale per costruzione dovrebbero stare circa tutte le velocita' (consigliato 22)" << endl;
	cin >> Ns;
	};

	cout << "Si vuole tenere d'occhio il moto della prima particella generata? (1 si', altro no)" << endl;
	cin >> brow;

	cout << "Inserire quanti atomi generare nella scatola 1x1" << endl;
	cin >> N;
	cout << "Inserire parametro sigma del modello di lennard jones" << endl;
	cin >> sigma;
	cout << "Inserire parametro epsilon del modello di lennard jones (il modulo, si tiene gia' conto del fatto che sia negativo)" << endl;
	cin >> epsilon;
	cout << "Inserire massa atomi" << endl;
	cin >> m;
	cout << "Inserire numero step temporali" << endl;
	cin >> Nt;
	cout << "Inserire larghezza step temporali" << endl;
	cin >> dt;

							// GENERAZIONE CONFIGURAZIONE INIZIALE RANDOMICA

	double r1,r2,r3,r4;

	if (G==1)
	{r1 = 0.5; 
	r2 =  0.5; 
	r3 = 0; 
	r4 = 0;}
 	else {
	r1 = 0.5; 
	r2 = 0.5; 
	r3 = B.Single_Number()-0.5; 
	r4 = B.Single_Number()-0.5;
	};

	atomo primo (r1,r2,r3,r4,0,0,0,0,m);

	V0.push_back(primo);

	while(n<N)
		{
		int TEST = 0;
		vector<double> D;

		if (G==1)
		{
		r1 = B.Single_Number(); 
		r2 = 0.4*B.Single_Number(); 
		r3 = 0; 
		r4 = 0;}
 		else {
		r1 = B.Single_Number(); 
		r2 = B.Single_Number();
		r3 = B.Single_Number()-0.5; 
		r4 = B.Single_Number()-0.5;
		};

		atomo X (r1,r2,r3,r4,0,0,0,0,m);


//Calcolo distranze tra atomo generato e tutti quelli già generati, motivo per il quale il primo è fuori dal ciclo
			for(int i=0; i<V0.size(); i++)
			{
			atomo Xi = V0.at(i);
			
			vector<double> ri = Rij(X, Xi);
			double di = modulov(ri);

			D.push_back(di);
			};

			for(int i=0; i<D.size(); i++)
			{
				if (sigma > D.at(i))
				{TEST++;}
				else if (sigma <= D.at(i)) {TEST=TEST;};
			};

			if (TEST==0)
			{
			cout << "Generato atomo " << n+1 << " in " << nn << " tentativi" << endl;
			V0.push_back(X);
			n++;
			nn=1;
			}
			else 
			{nn++;};

		D.clear();
		};

							// EVOLUZIONE TEMPORALE

	vector <atomo> V;
	
//cose per il tempo di esecuzione del programma
		clock_t start=0, end;
   		long double a=1000000;
		long double cpu_time_used=0;
		start = clock();

	//Calcolo l'interazione tra gli atomi nella configurazione iniziale fuori dal ciclo perché è molto più comodo

	vector<double> rij;
	double dij;

	for(int i=0;i<V0.size();i++) 
	{
	atomo Xi = V0.at(i);
	double f0x=0;				
	double f0y = -m*g;
	double UC0=0;
	double UG0 = m*g*Xi.y();

		for(int j=0;j<V0.size();j++)
			{
			atomo Xj = V0.at(j);

//vedi header file
			rij = Rij(Xi,Xj);
			dij = modulov(rij);

				if(dij <= 3.*sigma && dij > sigma/100.)
				{
				double rx = rij.at(0);
				double ry = rij.at(1);

				f0x = f0x + 24.0*epsilon*(2.*pow(sigma/dij , 12.) - pow(sigma/dij , 6.))*(rx/pow(dij , 2.)); //(86)
				f0y = f0y + 24.0*epsilon*(2.*pow(sigma/dij , 12.) - pow(sigma/dij , 6.))*(ry/pow(dij , 2.)); //(86)
				UC0 = UC0 + 4.*epsilon*(pow(sigma/dij , 12.) - pow(sigma/dij , 6.)); //(84)

				}
				else
				double null;

			};

//I nuovi atomi hanno stessa posizione e velocità della configurazione, ma sono state aggiornate le altre variabili che prima erano zero! Non è ancora stata fatta alcuna evoluzione temporale
	atomo X (Xi.x(),Xi.y(),Xi.vx(),Xi.vy(),f0x,f0y,UC0,UG0,m);

	V.push_back(X);
	};

	V0=V;
	V.clear();	

//Ora uso il propagatore di Velocity-Verlet, contenuto nella funzione Iter

	ofstream fout,eout,ekout,vout,bout;
	double Ng;

	cout << "Ogni quanti step temporali prendere frame per la gif?" << endl;
	cin >> Ng;

	fout.open("scatola.txt");
	eout.open("energie.txt");
	ekout.open("ek.txt");

	if(TEST==1)
	{
	vout.open("velocita.txt");
	};

	if(brow==1)
	{
	bout.open("motobrowniano.txt");
	};

//if TEST==1 sono le parti di codice dedicate alla gif con gli istogrammi delle velocità
	
						vector<double> modvel;
						double* medvel = new double[Ns];
						double e1 = 0.;
						double e2 = 1.2;	
						double h=(e2-e1)/(Ns-1.);

		for(int j=0; j<Nt; j++)
		{	
		V = Iter(sigma, epsilon, m, dt, V0, g);

			if((j+1.)/Ng == floor((j+1.)/Ng)) //solo in questo caso prendo giù i dati per la gif
			{

				for(int i=0; i<V.size(); i++)
					{
					atomo Xi = V.at(i);

						if(TEST==1)
						{
						double modv = sqrt(pow(Xi.vx(),2.) + pow(Xi.vy(),2.));
						modvel.push_back(modv);
						};

					fout << Xi.x() << "	" << Xi.y() << endl;	
					};

						if(TEST==1)
								{
								vector<int> Hvel = Histo_double(modvel,e1,e2,Ns);

								for(int i=0;i<Ns;i++)
								{
								vout << i*h << "	" << Hvel.at(i) << endl;
								double b = Ng/Nt;
								medvel[i] = medvel[i] + Hvel.at(i)*b; 
								};
								modvel.clear();
								};

					fout << endl;
					vout << endl;
					cout << "Passo " << j << endl;
			}

		double E=0,Ek=0;

//Calcolo l'energia ad ogni step per poi farci il grafico
		for(int i=0;i<V.size();i++) 
		{
		atomo Xi = V.at(i);
		E = E + Xi.EK() + Xi.Uc()/2. + Xi.Ug(); //(85)	
		Ek = Ek + Xi.EK(); 
		};

		eout << j << "	" << E << endl;
		ekout << j << "	" << Ek << endl;


	if(brow==1)
			{
			atomo A = V0.at(0);
			bout << A.x() << "	" << A.y() << endl;
			};

//ora che ho fatto i conti la mia configurazione vecchia viene aggiornata con quella nuova, in questo modo il programma non tiene nulla o quasi in memoria
		V0=V;
		V.clear();
		};


	fout.close();
	vout.close();
	eout.close();
	ekout.close();

	
			if(TEST==1)
			{
			ofstream vel;

			vel.open("velocitamedie.txt");

			for(int i=1;i<Ns;i++) 
				{
				vel << (i-0.5)*h << "	" << medvel[i] << endl;
				};

			vel.close();
			};


	end =clock();
	cpu_time_used = ((double) (end - start)) / a;

	cout << "Il tempo di calcolo del programma (nel computer di nicho) e' stato circa: " << cpu_time_used << " s" << endl;

	system("gnuplot -p -e \"call 'istruzioniscatola.txt'\""); 
	system("gnuplot -p -e \"call 'istruzioniscatolaE.txt'\""); 

	if(brow==1)
	{  
	system("gnuplot -p -e \"call 'istruzioniscatolaM.txt'\"");
	};


	if(TEST==1)
	{  
	system("gnuplot -p -e \"call 'istruzioniscatolaVm.txt'\"");
	system("gnuplot -p -e \"call 'istruzioniscatolaV.txt'\""); 
	};


	
	return 0;
}
