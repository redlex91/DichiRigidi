/*
 *************************************************************************************************************************************************************************
**************************************************************************************************************************************************************************
*
* In questo header sono definiti i passi per l'inizializzazione dei vettori posizione e velocità del sistema di dischi rigidi. La funzione CI( ... ) prende in ingresso  
* gli indirizzi di quattro vettori definiti nel main (rx, ry, vx, vy), in modo da poterli modificare. Gli indirizzi vengono passati semplicemente scrivendo il loro nome  
*
* e.g.: CI( rx, ry, vx, vy ) se nel main si vogliono passare i vettori tipo double  rx, ry, vx, vy.
*
* Nel main deve essere definita la lunghezza N  di tali vettori.
* 
**************************************************************************************************************************************************************************
**************************************************************************************************************************************************************************
*/

double speed_gen( void );
void CI( double rx[ N ], double ry[ N ], double vx[ N ], double vy[ N ] );

// genera le velocità casualmente
double speed_gen( void ){
	return ( rand()/(double)RAND_MAX ) * pow( ( -1 ), ( rand() % 2 ) );
}
