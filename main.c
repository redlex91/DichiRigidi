
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <float.h>
#include <limits.h>

// DEFINIZIONE DI COSTANTI
#define N 400
#define NP 17
#define ETA 0.45
#define PI M_PI //3.14159265358979

// VARIABILI GLOBALI
double delta_t; // intervallo fra due passi
FILE *cfPtr;
FILE *cfPtrEn;
double sigma;
int count = 0;
float press[ NP ];

// PROTOTIPI DELLE FUNZIONI
double speed_gen( void );
double coll_time_search( double r1x, double r1y, double v1x, double v1y, double r2x, double r2y, double v2x, double v2y, double s );
double v_min_search( double v[ 9 ]);
double amongst_copies_search( double r1x, double r1y, double v1x, double v1y, double r2x, double r2y, double v2x, double v2y, double s );
void evolve( double rx[ N ], double ry[ N ], double vx[ N ], double vy[ N ], double delta );
void speed_refresh( double rx[ N ], double ry[ N ], double vx[ N ], double vy[ N ], int p1, int p2 ); // p1 e p2 sono risp. le posizioni delle due particelle che collidono
double energy_compute( double vx[ N ], double vy[ N ] );
double distance( double r1x, double r1y, double r2x, double r2y );
int pressure( int nc );

int main( void ){
	 

	
	pressure( 500 );
	

	
	return 0;
}

int pressure( int nc ){

	/******************************************************* DICHIARAZIONE DELLE VARIABILI ***************************************************************************
	******************************************************************************************************************************************************************
	*****************************************************************************************************************************************************************/
	int n = sqrt( N );

	double Ei, Ef = 0;
	// variabili per l'inizializzazione
	double vx[ N ], vy[ N ], rx[ N ], ry[ N ];
	sigma = sqrt( ( 4.0 * ETA )/( PI * N ) ); // diametro dei dischi
	// double v0x[N], v0y[N], r0x[N], r0y[N]; // posizioni e velocità iniziali
	double Mtime[N][N]; // matrice dei tempi
	double step = 1/(double)n; // distanza fra i centri di dischi nella configurazione iniziale
	double vcmx=0, vcmy=0; // velocità del centro di massa
	double ctime, min;
	//double r0x[ N ], r0y[ N ]; // coordinate dei punti durante l'urto, per calcolo velocità
	// variabili per l'evoluzione
	double t0 = 0; // tempo iniziale al passo i-esimo
	
	// puntatore a file velocità
	FILE * cfPtrVel;
	cfPtrVel = fopen( "speed.dat", "a" );
	
	// variabili per il calcolo della pressione
	double delta_v = 0, v0x, v0y, intT = 0, P_mod = 0;
	FILE *cfPtrPr;

	// indici per cicli
	int i, j, k, l, m;
	int P1, P2;

	// gestione file
	cfPtr = fopen( "data.dat", "w" );
	cfPtrEn = fopen( "energy.dat", "w" );
	cfPtrPr = fopen( "pressure.dat", "w" );

	srand( time( NULL ) );// seme della funzione random
	
	/* printf( "\n#%g#\n", sigma ); */

	// variabile contatore di pressure( )
	
	
	// imposto la differenza di velocità delle particelle che collidono a 0
	
	/***************************************************************** ASSEGNAZIONE DELLE COORDINATE *****************************************************************
	******************************************************************************************************************************************************************
	*****************************************************************************************************************************************************************/

	// coordinate x dei punti
	for( i=0; i<N; i++ ){
		rx[i] = ( ( i % n ) * step );
	}	
	// coordinate y dei punti
	for( i=0; i<N; i++ ){
		ry[i] = ( i / n ) * step;
	}

	// assegno le velocità a caso
	for( i=0; i<N; i++ ){
		vx[i] = speed_gen();
		vy[i] = speed_gen();
	}
	
	// calcolo la velocità del centro di massa
	for( i=0; i<N; i++ ){
		vcmx += (double)vx[i] / N;
		vcmy += (double)vy[i] / N;
	}
	// riassegno le velocità con la condizione vcm=0
	for( i=0; i<N; i++ ){
		vx[i] = vx[i] - vcmx;
		vy[i] = vy[i] - vcmy;
	}
	
	vcmx=0, vcmy=0;
	// controllo che il centro di massa sia fermo
		for( i=0; i<N; i++ ){
		vcmx += (double)vx[i] / N;
		vcmy += (double)vy[i] / N;
	}
	
	fprintf( cfPtrEn, "Passo\t\t\tEnergia\n\n" );
	
	// scalo le velocità in modo da normalizzare l'energia del sistema
	Ei = energy_compute( vx, vy );
	
	for( i=0; i<N; i++ ){
		
			vx[ i ] = vx[ i ] / sqrt( Ei );
			vy[ i ] = vy[ i ] / sqrt( Ei );
	}
	
	
	// inizia il ciclo di evoluzione del sistema
	while( count <= nc ){

	fprintf( cfPtrEn, "%d\t\t\t%f\n", count, energy_compute( vx, vy ) );
	/*
	printf( "\n Esempio= %f\n" , amongst_copies_search( r0x[1], r0y[1], v0x[1], v0y[1], r0x[0], r0y[0], v0x[0], v0y[0], sigma ) );
	// controlli vari
	printf( "%f,%f\n\n", vcmx, vcmy );*/
	/*for( i=0; i<N; i++ ){
		fprintf( cfPtr, "posizioni:\t(%f,%f)", rx[i], ry[i] );
		fprintf( cfPtr, "\t velocità:\t(%f,%f)", vx[i], vy[i] );
		fprintf( cfPtr,  "\n\n" );
	}*/
			
	/**************************************************************** DETERMINAZIONE DEI TEMPI DI COLLISIONE *********************************************************
	******************************************************************************************************************************************************************
	*****************************************************************************************************************************************************************/		
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			ctime = amongst_copies_search( rx[i], ry[i], vx[i], vy[i], rx[j], ry[j], vx[j], vy[j], sigma );
			if( i==j ) Mtime[i][j] = 0;
			else if( ctime == FLT_MAX ) Mtime[i][j] = INFINITY;
			else Mtime[i][j] = ctime;
		}
	}
	
	
	/*fprintf( cfPtr, "\n Tabella\n -------------------------------------------\n" );
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			fprintf( cfPtr, "\t%.2f", Mtime[i][j] );
		}
		fprintf( cfPtr, "\n" );
	}
	fprintf( cfPtr, "\n -------------------------------------------" );
	*/
	
	// determino il tempo di collisione più piccolo


	/***************************** algoritmo di scansione completa della matrice dei tempi
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			if( Mtime[i][j] < min && i != j ){
				min = Mtime[ i ][ j ];
				P1 = i;
				P2 = j;
			}
		}
	}
	printf( "\n min = %g.\n", min );
	*/


	min = INFINITY;
	for( i=0; i<N; i++ ){
		for( j=i; j<N; j++ ){
			if( Mtime[i][j] < min && i != j ){
				min = Mtime[ i ][ j ];
				P1 = i;	
				P2 = j;
			}
		}
	}

	fprintf( cfPtr, "\nIl disco %d collide col disco %d al tempo t=%E.\n", P1, P2, min );
	
	
	if( cfPtr == NULL ) printf(" Memoria non disponibile!\n ");
	else{
		fprintf( cfPtr, "\t x1\t y1\t x2\t y2 \t vx1\t vy1\t vx2\t vy2\n" );
		fprintf( cfPtr, "\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\n", rx[ P1 ], ry[ P1 ], rx[ P2 ], ry[ P2 ], vx[ P1 ], vy[ P1 ], vx[ P2 ], vy[ P2 ] );
		fprintf( cfPtr, "\t Calcolo la distanza fra le particelle che collidono:\t %f $ %f.\n", sqrt( pow( rx[ P1 ] - rx[ P2 ], 2 ) + pow( ry[ P1 ] - ry[ P2 ], 2 ) ) - sigma, distance( rx[ P1 ], ry[ P1 ], rx[ P2 ], ry[ P2 ] ) - sigma );
		fprintf( cfPtr, "\n\n\n" );
		} // chiudo else - scrittura file
	//fclose( cfPtr );
		
	delta_t = min;
	
	/*
	// determino la posizione della particella 2 coinvolta nell'urto
	if( abs( rx[ P1 ] - rx[ P2 ]) > 0.5 )
		rx[ P2 ] += floor( 0.5 +  rx[ P1 ] - rx[ P2 ] );
 
	if( abs( ry[ P1 ] - ry[ P2 ]) > 0.5 )
		ry[ P2 ] += floor( 0.5 +  ry[ P1 ] - ry[ P2 ] );
	

	//cfPtr = fopen( "data.dat", "a" );
	if( cfPtr == NULL ) printf(" Memoria non disponibile!\n ");
	else{
		fprintf( cfPtr, "\t x1\t y1\t x2\t y2 \t vx1\t vy1\t vx2\t vy2\n" );
		fprintf( cfPtr, "\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\n", rx[ P1 ], ry[ P1 ], rx[ P2 ], ry[ P2 ], vx[ P1 ], vy[ P1 ], vx[ P2 ], vy[ P2 ] );
		fprintf( cfPtr, "\t Calcolo la distanza fra le particelle che collidono dopo l'aggiustamento:\t %f.\n",  sqrt( pow( rx[ P1 ] - rx[ P2 ], 2 ) + pow(ry[ P1 ] - ry[ P2 ], 2 ) ) - sigma );
		fprintf( cfPtr, "\t Energia iniziale del sistema: %f.\n", energy_compute( vx, vy ) );

		fprintf( cfPtr, "\n\n\n" );
	}*/
	
	/************************************************************* INIZIO EVOLUZIONE DEL SISTEMA *********************************************************************
	******************************************************************************************************************************************************************
	******************************************************************************************************************************************************************
	*****************************************************************************************************************************************************************/
	evolve( rx, ry, vx, vy, delta_t );


/*	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			if( sqrt( pow((rx[i] - rx[j]), 2)+ pow(( ry[i] - ry[ j ] ),2) ) - sigma <0 ) printf("\nPasso %d: %f\n",count,  sqrt( pow((rx[i] - rx[j]), 2)+ pow(( ry[i] - ry[ j ] ),2) ) - sigma );
		}
	}	
*/	

	if( cfPtr == NULL ) printf(" Memoria non disponibile!\n ");
	else{
		fprintf( cfPtr, "\t\t--- EVOLUZIONE ---\n\n" );
		fprintf( cfPtr, "\t Calcolo la distanza fra le particelle che collidono dopo l'evoluzione:\t %f $ %f.\n",  sqrt( pow( rx[ P1 ] - rx[ P2 ], 2 ) + pow(ry[ P1 ] - ry[ P2 ], 2 ) ) - sigma, distance( rx[ P1 ], ry[ P1 ], rx[ P2 ], ry[ P2 ] ) - sigma );
		fprintf( cfPtr, "\t Energia iniziale del sistema: %f.\n", Ei = energy_compute( vx, vy ) );
		
		fprintf( cfPtr, "\n\n\n" );
	}

	// memorizzo le posizioni delle particelle durante l'urto
	/*for( k = 0; k < N; k++ ){
		r0x[ k ] = rx[ k ];
		r0y[ k ] = ry[ k ];
	}*/
	
	// memorizzo le velocità delle particelle che collidono prima dell'urto per il calcolo di delta_v (per la pressione)
	v0x = vx[ P1 ];
	v0y = vy[ P1 ];

	/************************************************************* AGGIORNAMENTO DELLE VELOCITA **********************************************************************
	******************************************************************************************************************************************************************
	******************************************************************************************************************************************************************
	*****************************************************************************************************************************************************************/
	speed_refresh( rx, ry, vx, vy,  P1, P2 );	
	
	if( cfPtr == NULL ) printf(" Memoria non disponibile!\n ");
	else{
		fprintf( cfPtr, "\t x1\t y1\t x2\t y2 \t vx1\t vy1\t vx2\t vy2\n" );
		fprintf( cfPtr, "\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\n", rx[ P1 ], ry[ P1 ], rx[ P2 ], ry[ P2 ], vx[ P1 ], vy[ P1 ], vx[ P2 ], vy[ P2 ] );
		fprintf( cfPtr, "\t Calcolo la distanza fra le particelle che collidono dopo l'evoluzione:\t %f $$ %f.\n",  sqrt( pow( rx[ P1 ] - rx[ P2 ], 2 ) + pow( ry[ P1 ] - ry[ P2 ], 2 ) ) - sigma, distance( rx[ P1 ], ry[ P1 ], rx[ P2 ], ry[ P2 ] ) - sigma );
	
		Ef = energy_compute( vx, vy );
		fprintf( cfPtr, "\t Energia iniziale - Energia finale del sistema: %f.\n",  fabs( Ef - Ei ) );
		} // chiudo else - scrittura file
	
	
	// calcolo la pressione
	delta_v += sqrt( pow( ( v0x - vx[ P1 ] ), 2 ) + pow( ( v0y - vy[ P2 ] ), 2 ) );
	intT += delta_t;
	P_mod = 1.0 + ( sigma )/( 2.0 * Ef * intT ) * delta_v;
	fprintf( cfPtrPr, "\n%f\t%f\t%f ", delta_v, intT,  P_mod );
	
	// ---> passo successivo
	
	if( count == nc ){
		for( i=0; i<N; i++ ){
			fprintf( cfPtrVel, "\n%f\t%f", vx[ i ], vy[ i ] );
		}
	}
	
	count++;
	
	/*
	for( i=0; i<N; i++ ){
		for( j=i+1; j<N; j++ ){
			if( sqrt( pow((rx[i] - rx[j]), 2)+ pow(( ry[i] - ry[ j ] ),2) ) - sigma < FLT_EPSILON ) printf("\nPasso %d: %e\n",count,  sqrt( pow((rx[i] - rx[j]), 2)+ pow(( ry[i] - ry[ j ] ),2) ) - sigma );
		}
	}*/
	
	
	} // chiudo il while di pressure( )
	
	fclose( cfPtr );
	fclose( cfPtrEn );
	fclose( cfPtrPr );
	for( i=0; i<N; i++ ){
		for( j=0; j<N; j++ ){
			if( ( sqrt( pow((rx[i] - rx[j]), 2) + pow(( ry[i] - ry[ j ] ),2) )  - sigma && i!=j ) < - pow( 10, -10 ) ) printf("\nPasso %d: sovrapposizione: %f\nIndici: %d, %d",count,  sqrt( pow((rx[i] - rx[j]), 2)+ pow(( ry[i] - ry[ j ] ),2) ) - sigma, i, j );
			
		}
	}	

	return 0;
}

// DICHIARAZIONE DELLE FUNZIONI

/************************************************************************************************************************************************************************
 ************************************************************************************************************************************************************************
 *                                                   funzioni per la determinazione dei tempi di collisione
 ************************************************************************************************************************************************************************
 ************************************************************************************************************************************************************************
 */
// genera le velocità casualmente
double speed_gen( void ){
	
	return ( rand()/(double)RAND_MAX ) * pow( ( -1 ), ( rand() % 2 ) );

}


// determina il tempo di collisione fra due dischi
double coll_time_search( double r1x, double r1y, double v1x, double v1y, double r2x, double r2y, double v2x, double v2y, double s ){
	double dr_x, dr_y, dv_x, dv_y, dot_prod, sq_norm_r, sq_norm_v;
	int l, m;
	
	// calcoli intermedi
	dr_x = r1x - r2x;
	dr_y = r1y - r2y;
	dv_x = v1x - v2x;
	dv_y = v1y - v2y;
	
	dot_prod = dr_x * dv_x + dr_y * dv_y;
	sq_norm_r = dr_x * dr_x + dr_y * dr_y;
	sq_norm_v = dv_x * dv_x + dv_y * dv_y;
	
	// uscita della funzione
	
	/*fprintf( cfPtr, "\n Cerca: %g, %g, %g, %e, %g\n", dot_prod, sq_norm_r, sq_norm_v, s * s, sqrt(sq_norm_r) - s );
	fprintf( cfPtr, " tempo: \%g\n", (- dot_prod - sqrt( dot_prod * dot_prod - sq_norm_v * ( sq_norm_r - s * s ) ) ) / sq_norm_v );*/
	
	
	
	if( dot_prod < 0 && ( dot_prod * dot_prod - sq_norm_v * ( sq_norm_r - s * s ) >=0 ) ){
 
		return ( - dot_prod - sqrt( dot_prod * dot_prod - sq_norm_v * ( sq_norm_r - s * s ) ) ) / sq_norm_v;
 	
	}
	else return FLT_MAX;
	
}

// determina il minimo di un vettore
double v_min_search( double v[ 9 ]){
	int k = 0;
	double min = v[0];
	
	for( k=1; k<9; k++ ){
		if( min > v[k] ) min = v[k];
	}
	
	return min;
}

// cerca il tempo di collisione di un disco con le 9 repliche di un secondo
double amongst_copies_search( double r1x, double r1y, double v1x, double v1y, double r2x, double r2y, double v2x, double v2y, double s ){
	int p, q, k=0;
	double vec[ 9 ];
	
	for( p=-1; p<=1; p++){
		for( q=-1; q<=1; q++ ){
			vec[ k ] = coll_time_search( r1x, r1y, v1x, v1y, r2x + p, r2y + q, v2x, v2y, s );
			k++;
		}
	}
	
	/*
	k--;
	printf( "\n\n" );
	while( k>=0 ){
		printf( "%d___ %e ___\n", k, vec[k] );
		k--;
	}
	printf( "\n\n" );
	*/
	
	return v_min_search( vec );
}

/*************************************************************************************************************************************************************************
 *************************************************************************************************************************************************************************
                                                           funzioni per l'evoluzione del sistema
 *************************************************************************************************************************************************************************
 ************************************************************************************************************************************************************************/

// evoluzione del sistema
void evolve( double rx[ N ], double ry[ N ], double v0x[ N ], double v0y[ N ], double delta ){
	int k; // contatore del ciclo

	for( k = 0; k < N; k++ ){
		rx[ k ] += v0x[ k ] * delta;
		ry[ k ] += v0y[ k ] * delta;

		// proietto le posizioni dentro alla scatola
		rx[ k ] -= floor( rx[ k ] );
		ry[ k ] -= floor( ry[ k ] ); 
	}
}

// aggiornamento delle velocità
void speed_refresh( double rx[ N ], double ry[ N ], double vx[ N ], double vy[ N ], int p1, int p2){
	double vers_rx, vers_ry, dot_prod, sq_norm;
	double x1 = rx[ p1 ], x2 = rx[ p2 ], y1 = ry[ p1 ], y2 = ry[ p2 ]; 

	//determino la posizione della particella 2 coinvolta nell'urto
	/*if( fabs( x1 - x2 ) > 0.5 )
		x2 = ( x1 - x2 ) / fabs( x1 - x2 ) * ( 1 - fabs( x1 - x2 ) );
		//rx[ p2 ] += floor( 0.5 +  rx[ p1 ] - rx[ p2 ] );
 
	if( fabs( y1 - y2 ) > 0.5 )
		y2 = ( y1 - y2 )/ fabs( y1 - y2 ) * ( 1 - fabs( y1 - y2 ) );
		//ry[ p2 ] += floor( 0.5 +  ry[ p1 ] - ry[ p2 ] );*/

	if( fabs( x1 - x2 ) > 0.5 )
		x2 += floor( 0.5 +  x1 - x2 );
	if( fabs( y1 - y2 ) > 0.5 )
		y2 += floor( 0.5 +  y1 - y2 );

	// qualche calcolo preliminare	
	sq_norm = ( x2 - x1 ) * ( x2 -x1 ) + ( y2 -y1 ) * ( y2 -y1 );
	vers_rx = ( x2 - x1 ) / sqrt( sq_norm );
	vers_ry = ( y2 -y1 ) / sqrt( sq_norm );
	dot_prod = ( vx[ p1 ] - vx[ p2 ] ) * vers_rx + ( vy[ p1 ] -vy[ p2 ] ) * vers_ry;
	
	// assegno le nuove velocità
	vx[ p1 ] -= dot_prod * vers_rx;
	vy[ p1 ] -= dot_prod * vers_ry;
	vx[ p2 ] += dot_prod * vers_rx;
	vy[ p2 ] += dot_prod * vers_ry;
}

/*************************************************************************************************************************************************************************
 *************************************************************************************************************************************************************************
                                                            funzioni per calcoli ricorrenti
 *************************************************************************************************************************************************************************
 ************************************************************************************************************************************************************************/

// calcolo dell'energia totale del sistema
double energy_compute( double vx[ N ], double vy[ N ] ){

	double sum=0;
	int k;

	for( k=0; k<N; k++ )	sum += 1./2 * ( pow( vx[ k ], 2 ) + pow( vy[ k ], 2 ) );
	return sum;
}

// calcolo delle distanza fra due particelle tenendo conto delle condizioni di bordo
double distance( double r1x, double r1y, double r2x, double r2y ){
	
	double x1 = r1x, x2 = r2x, y1 = r1y, y2 = r2y;
	
	/*if( fabs( x1 - x2 ) > 0.5 )
		x2 = ( x1 - x2 )/fabs( x1 - x2 ) * ( 1 - fabs( x1 - x2 ) );
	if( fabs( y1 - y2 ) > 0.5 )
		y2 = ( y1 - y2  )/ fabs( y1 - y2 ) * ( 1 -fabs( y1 - y2 ) );*/
		
	if( fabs( x1 - x2 ) > 0.5 )
		x2 += floor( 0.5 +  x1 - x2 );
	if( fabs( y1 - y2 ) > 0.5 )
		y2 += floor( 0.5 +  y1 - y2 );
	
	return sqrt( pow( (x1 - x2 ), 2 ) + pow( (y1 - y2 ), 2 ) );
}
