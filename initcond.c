/*
 *************************************************************************************************************************************************************************
**************************************************************************************************************************************************************************
*
* In questo file sono definiti i passi per l'inizializzazione dei vettori posizione e velocità del sistema di dischi rigidi. La funzione CI( ... ) prende in ingresso  
* gli indirizzi di quattro vettori definiti nel main (rx, ry, vx, vy), in modo da poterli modificare. Gli indirizzi vengono passati semplicemente scrivendo il loro nome  
*
* e.g.: CI( rx, ry, vx, vy ) se nel main si vogliono passare i vettori tipo double  rx, ry, vx, vy.
*
* Nel main deve essere definita la lunghezza N  di tali vettori.
* 
**************************************************************************************************************************************************************************
**************************************************************************************************************************************************************************
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include "initcond.h"

/* funzione principale di initcond.h: genera le condizioni iniziali del sistema */
void CI( double *rx, double *ry, double *vx, double *vy, int N ){
	
	/* DEFINIZIONE VARIABILI DI CI( ... ) */
	int n = sqrt( N );
	double step = 1/(double)n; // distanza fra i centri di dischi nella configurazione iniziale
	double vcmx=0, vcmy=0; // velocità del centro di massa
	int i;
	
	srand( time( NULL ) ); // imposto il seme della funzione random

	/* assegno coordinate dei punti */
	for( i=0; i<N; i++ ){
		rx[i] = ( ( i % n ) * step );
	}	
	for( i=0; i<N; i++ ){
		ry[i] = ( i / n ) * step;
	}

	/* assegno le velocità */
	for( i=0; i<N; i++ ){
		vx[i] = speed_gen();
		vy[i] = speed_gen();
	}
	
	/* calcolo la velocità del centro di massa */
	for( i=0; i<N; i++ ){
		vcmx += (double)vx[i] / N;
		vcmy += (double)vy[i] / N;
	}
	
	/* riassegno le velocità con la condizione vcm=0 */
	for( i=0; i<N; i++ ){
		vx[i] = vx[i] - vcmx;
		vy[i] = vy[i] - vcmy;
	}
	
	
	/* controllo che il centro di massa sia fermo */
	vcmx = 0, vcmy = 0;
	for( i=0; i<N; i++ ){
		vcmx += (double)vx[i] / N;
		vcmy += (double)vy[i] / N;
	}
	if( vcmx > pow( 10, -10 ) || vcmy > pow( 10, -10 ) ) abort( ); // se il cdm no è fermo, interrompo l'esecuzione del programma
	// printf( "Velcità centro di massa:\t(%e,%e).\n", vcmx, vcmy );
	

}


/* funzione che genera le velocità casualmente */
double speed_gen( void ){
	
	return ( rand()/(double)RAND_MAX ) * pow( ( -1 ), ( rand() % 2 ) );

}

void printCI( double *rx, double *ry, double *vx, double *vy, int N ){
	
	int i;
	
	for( i=0; i<N; i++ ){
		printf( "posizioni:\t(%f,%f)", rx[i], ry[i] );
		printf( "\t velocità:\t(%f,%f)", vx[i], vy[i] );
		printf( "\n\n" );
	}
		
}
