#include "Poly.h"




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


    // =================Partie complétion slice 0=========================
    nombreValeur = zVals.n_elem;
    laguerre = arma::cube(m , nombreValeur,n, arma::fill::zeros);

    arma::mat slice0 = arma::mat(m, nombreValeur, arma::fill::zeros); // voir si y'a pas moyen de remplir cette slice directement de 1
    slice0.fill(1);



    laguerre.slice(0)= slice0;

    // ===============================================================


 
    // =================Partie complétion slice 1=========================
   //je fais un vect de 1+m 
    //pui pour chaque je  fais un slice avec ce vecteur ajouter à chaqu, évite double boucle


    arma::vec vecTemporaire(m);
    for (int i =0 ; i<m ; i++)
    {
           vecTemporaire(i) = 1 + i;
    }
    

    arma::mat slice1 = arma::mat(m, nombreValeur);
    for (int nu = 0; nu < nombreValeur; nu++) {
        arma::vec tempo(m);
        tempo.fill(nu);
        slice1.col(nu) = vecTemporaire + tempo;     }
    laguerre.slice(1) = slice1;

       
    // =========================== Partie complétion autres Slices ============================

    arma::mat mat1 = arma::mat(m, nombreValeur, arma::fill::ones); // pour éviter de les récréer à chaque itération
    arma::mat mat2 = arma::mat(m, nombreValeur);
    mat2.fill(2);


    matFact1 = arma::mat(m, nombreValeur); // Ca aussi c'est indep de n donc autant les créer là, en vrai go en faire des attributs privé ensuite, libère la fct
    matFact2 = arma::mat(m, nombreValeur);


    arma::vec vecTemporaireColonne(m);
    for (int i =0 ; i<m ; i++)
    {
           vecTemporaireColonne(i) = i - 1;
    }

    
    for (int nu =0 ; nu<nombreValeur; nu++)
    {
        arma::vec nuVec(m);
        nuVec.fill(nu);
        matFact1.col(nu) = vecTemporaireColonne - nuVec;
        matFact2.col(nu) = vecTemporaireColonne;

    }

 
    for (int N = 2 ; N < n ; N++){
        laguerre.slice(N) = Poly::calcSliceN(N , laguerre.slice(N-1) , laguerre.slice(N-1) , mat1 , mat2);
    }

}
    // ===============================================================

    





void Poly::printMatrix(arma::mat mat) {
    int rows = mat.n_rows;
    int cols = mat.n_cols;

    std::cout << "[\n";
    for (int i = 0; i < rows; i++) {
        std::cout << "  [";
        for (int j = 0; j < cols; j++) {
            std::cout << mat(i, j) << " ";
        }
        std::cout << "]\n";
    }
    std::cout << "]\n";
}

arma::mat Poly::calcSliceN(int n, arma::mat sliceNMoins1, arma::mat sliceNMoins2 , arma::mat mat1 , arma::mat mat2)
{

    arma::mat sliceN = arma::mat(sliceNMoins1.n_rows, sliceNMoins1.n_cols , arma::fill::zeros);

    sliceN = (mat2 + (matFact1/ n)) % sliceNMoins1 - (mat1 + (matFact2/ n) % sliceNMoins2);

    printMatrix(sliceN);
    return (sliceN);
}
