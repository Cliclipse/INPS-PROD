#include "Hermit.h"



void Poly::calcHermit( int n , arma::rowvec zVals ) {

    nombreValeur = zVals.n_elem;
    hermit = arma::mat(n, nombreValeur, arma::fill::zeros);


    hermit.row(0).fill(1);
    if (ordreMax > 1) {
        hermit.row(1) = 2 * zVals;

        for (int ligne = 2 ; ligne < ordreMax + 1 ; ligne++ ){
        hermit.row(ligne) = 2.0 * (zVals % hermit.row(ligne - 1)) - 2.0 * (ligne - 1) * hermit.row(ligne - 2);
        }
    }
}


void Poly::calcLaguerre( int n , int m , arma::rowvec zVals ) {

    nombreValeur = zVals.n_elem;
    laguerre = arma::cube(n, m , nombreValeur, arma::fill::zeros);


    arma::mat slice0 = arma::mat(n, nombreValeur, arma::fill::0); // voir si y'a pas moyen de remplir cette slice directement de 1
    slice.fill(1);

    arma::vec vecTemporaire();
    for (int i =0 ; i<m ; i++)
    {
           vecTemporaire(i) = 1 + i;
    }


    arma::mat slice = arma::mat(n, nombreValeur, arma::fill::0); // voir si y'a pas moyen de remplir cette slice directement de 1
    slice.fill(1);


    //je fais un vect de 1+m 
    //pui pour chaje fais un slice avec ce vecteur ajouter à chaque ligne 



    laguerre.mat(0).row(1) = 1 + m - zVals 

    if (ordreMax > 1) {
        hermit.row(1) = 2 * zVals;

        for (int ligne = 2 ; ligne < ordreMax + 1 ; ligne++){
        hermit.row(ligne) = 2.0 * (zVals % hermit.row(ligne - 1)) - 2.0 * (ligne - 1) * hermit.row(ligne - 2);
        }
    }



}
