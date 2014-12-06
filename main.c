#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>

#include "initcond.h"

// DEFINIZIONE DI COSTANTI
#define N 4
#define ETA 0.55
#define PI M_PI //3.14159265358979

// VARIABILI GLOBALI
double delta_t; // intervallo fra due passi
double sigma; // diametro dei dishi

FILE *cfPtr; // puntatore a file data
FILE *cfPtrEn; // puntatore a file energy



int main( void ){
	
	// variabili per l'inizializzazione
	double vx[ N ], vy[ N ], rx[ N ], ry[ N ]; // vettori posizione e velocità
	sigma = sqrt( ( 4.0 * ETA )/( PI * N ) ); // diametro dei dischi
	
	// assegnazione delle condizioni iniziali (e stampa opzionale)
	CI( rx, ry, vx, vy, N );
	printCI( rx, ry, vx, vy, N );
	
	return 0;
}
