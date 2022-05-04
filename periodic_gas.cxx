#include "gas.h"

int main()

{	
	int N, n=1, nn=1, Nt, Ns=0; //n ed nn servono per la configurazione iniziale delle posizioni
	double sigma,epsilon,dt,m; 	//costanti
	vector<atomo> V0;			//conterrà configurazione iniziale
	LCG B;						//vedi header file
	int TEST,brow;

	cout << "Si vuole tenere d'occhio la distribuzione del modulo delle velocita' tramite istogramma? (1 si', altro no)" << endl;
	cin >> TEST;

	if(TEST==1)
	{
	cout << "Inserire numero di step in cui dividere l'intervallo [0,1.5], nel quale per costruzione dovrebbero stare circa tutte le velocita' (consigliato 30)" << endl;
	//Poichè le velocità sono generate in un certo modo, l'intervallo adatto a tenere d'occhio il loro modulo è [0,1.2]. Se si volesse far scegliere all'operatore anche come generare le v [0,1.2] non andrebbe più bene e andrebbe scritto generalizzato. Mi sembrava inutile farlo visto che porterebbe solo complicazioni (velocità bassa? per vedere evolvere il sistema serve dt grande, no bene; allora mi servono moltissimi passi temporali, no bene. velocità alte? serve dt piccolissimo, ecc...)
	cin >> Ns;
	};

cout << "Si vuole tenere d'occhio il moto della prima particella generata? (1 si', altro no)" << endl;
	cin >> brow;

	cout << "Inserire quanti atomi generare con condizioni periodiche" << endl;
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

	// posizioni casuali nel quadrato [0,1]x[0,1] 
	// velocita' distribuite in modo "uniforme" senza boltzmann, sia perche' sono poche particelle, sia perche' poi dovrebbero distribuirsi a quel modo

//primo atomo generato fuori dal ciclo per comodità. Nel mezzo della griglia per vederne il moto browniano

	double r1 = 0.5; 
	double r2 = 0.5; 
	//double r3 = B.Single_Number()-0.5; 
	//double r4 = B.Single_Number()-0.5;
	double r3 = sign(B.Single_Number()-0.5)*(B.Single_Number()/3.+0.2); 
	double r4 = sign(B.Single_Number()-0.5)*(B.Single_Number()/3.+0.2); 

	atomo primo (r1,r2,r3,r4,0,0,0,m);

	V0.push_back(primo);

	while(n<N)
		{
		int TEST = 0;
		vector<double> D;
	
		r1 = B.Single_Number(); 
		r2 = B.Single_Number();
		r3 = sign(B.Single_Number()-0.5)*(B.Single_Number()/3.+0.2); 
		r4 = sign(B.Single_Number()-0.5)*(B.Single_Number()/3.+0.2);

		atomo X (r1,r2,r3,r4,0,0,0,m);


//Calcolo distranze tra atomo generato e tutti quelli già generati, motivo per il quale il primo è fuori dal ciclo
			for(int i=0; i<V0.size(); i++)
			{
			atomo Xi = V0.at(i);
			
			vector<double> ri = vmin(X, Xi);
			double di = modulov(ri);

			D.push_back(di);
			};

//Se l'atomo generato fosse troppo vicino a uno o più atomi non andrebbe bene, non solo non è realistico, ma l'energia potenziale sarebbe pressochè infinita all'inizio della simulazione. La configurazione iniziale fatta in modo serio utilizzerebbe l'algoritmo di metropolis (vedi dispense pagina 42-43) e la minimizzazione dell'energia potenziale del sistema, ma è inutilmente complesso per questo giochino
			for(int i=0; i<D.size(); i++)
			{
				if (sigma > D.at(i))
				{TEST++;}
				else if (sigma <= D.at(i)) {TEST=TEST;};
			};


//Questo è utile per capire se la sigma scelta va bene o no, data una certa sigma possono starci solo un certo numero di atomi, circa. Ad esempio per sigma = 0.1 ce ne stanno 66-68 massimo. Se si impianta vuol dire che non c'è più posto, quindi rilanciatelo con dei parametri migliori
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
	double f0y=0;
	double U0=0;

		for(int j=0;j<V0.size();j++)
			{
			atomo Xj = V0.at(j);

//vedi header file
			rij = vmin(Xi,Xj);
			dij = modulov(rij);

				if(dij <= 3.*sigma && dij > sigma/100.)
				{
				double rx = rij.at(0);
				double ry = rij.at(1);

				f0x = f0x + 24.0*epsilon*(2.*pow(sigma/dij , 12.) - pow(sigma/dij , 6.))*(rx/pow(dij , 2.)); //(86)
				f0y = f0y + 24.0*epsilon*(2.*pow(sigma/dij , 12.) - pow(sigma/dij , 6.))*(ry/pow(dij , 2.)); //(86)
				U0 = U0 + 4.*epsilon*(pow(sigma/dij , 12.) - pow(sigma/dij , 6.)); //(84)

				}
				else
				double null;
			};

//I nuovi atomi hanno stessa posizione e velocità della configurazione, ma sono state aggiornate le altre variabili che prima erano zero! Non è ancora stata fatta alcuna evoluzione temporale
	atomo X (Xi.x(),Xi.y(),Xi.vx(),Xi.vy(),f0x,f0y,U0,m);

	V.push_back(X);
	};

	V0=V;
	V.clear();	

//Preparativi per i file da scrivere ecc ecc

	ofstream fout,eout,ekout,vout,bout;
	double Ng;

	cout << "Ogni quanti step temporali prendere frame per la gif?" << endl;
	cin >> Ng;

	fout.open("gas.txt");
	eout.open("energie.txt");
	ekout.open("ek.txt");

	//if(TEST==1)
	//{
	//vout.open("velocita.txt");
	//};

	if(brow==1)
	{
	bout.open("motobrowniano.txt");
	};

//if TEST==1 sono le parti di codice dedicate alla gif con gli istogrammi delle velocità
	
								vector<double> modvel;
								double* medvel = new double[Ns];
								double e1 = 0.;			//estremi per l'asse x dell'istogramma per le velocità
								double e2 = 1.5;		//si potrebbe anche farli inserire a mano ma non avevo voglia
								double h=(e2-e1)/(Ns-1.);

//Ora uso il propagatore di Velocity-Verlet, contenuto nella funzione Iter

		for(int j=0; j<Nt; j++)
		{	
		V = Iter(sigma, epsilon, m, dt, V0);

		if((j+1.)/Ng == floor((j+1.)/Ng)) //solo in questo caso prendo giù i dati per la gif del gas
		{

			for(int i=0; i<V.size(); i++)
				{
				atomo Xi = V.at(i);

				fout << Xi.x() << "	" << Xi.y() << endl;	
				};

				fout << endl;
				cout << "Passo " << j << endl;
			};
								if(TEST==1)
								{

								for(int i=0; i<V.size(); i++)
								{
								atomo Xi = V.at(i);

								double modv = sqrt(pow(Xi.vx(),2.) + pow(Xi.vy(),2.));
								modvel.push_back(modv);
								};

								vector<int> Hvel = Histo_double(modvel,e1,e2,Ns);

								for(int i=0;i<Ns;i++)
								{
								//vout << i*h << "	" << Hvel.at(i) << endl;
								double b = 1./Nt;
								medvel[i] = medvel[i] + Hvel.at(i)*b; 
								};

								modvel.clear();
								};
								//vout << endl;
		double E=0,Ek=0;

//Calcolo l'energia ad ogni step per poi farci il grafico
		for(int i=0;i<V.size();i++) 
		{
		atomo Xi = V.at(i);
		E = E + Xi.EK() + Xi.U()/2.; //(85)	
		Ek = Ek + Xi.EK(); 
		};

		eout << j << "	" << E << endl;
		ekout << j << "	" << Ek << endl;


	if(brow==1)
			{
			atomo A = V0.at(0);
			bout << A.x() << "	" << A.y() << endl;
			};

//ora che ho fatto i conti e tutto la mia configurazione vecchia V0 viene aggiornata con quella nuova, in questo modo il programma non tiene nulla o quasi in memoria
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

			for(int i=0;i<Ns;i++) 
				{
				vel << i*h << "	" << medvel[i] << endl;
				};

			vel.close();
			};


	end =clock();
	cpu_time_used = ((double) (end - start)) / a;

	cout << "Il tempo di calcolo del programma (nel computer di nicho) e' stato circa: " << cpu_time_used << " s" << endl;

	system("gnuplot -p -e \"call 'istruzionigas.txt'\""); 
	system("gnuplot -p -e \"call 'istruzionigasE.txt'\""); 

	if(brow==1)
	{  
	system("gnuplot -p -e \"call 'istruzionigasM.txt'\"");
	};


	if(TEST==1)
	{  
	system("gnuplot -p -e \"call 'istruzionigasVm.txt'\"");
	//system("gnuplot -p -e \"call 'istruzionigasV.txt'\""); 
	};


	
	return 0;
}
