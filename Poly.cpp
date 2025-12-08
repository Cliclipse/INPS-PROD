#include "Poly.h"

#include <stdexcept>


void Poly::calcHermite(int n, const arma::vec & zVals)
{
    if (n < 0)
    {
        throw std::invalid_argument("Hermite order must be non-negative");
    }

    ordreMax = n;
    nombreValeur = zVals.n_elem;
    arma::rowvec zRow = zVals.t();

    // Allocate the Hermite table and seed the constant row
    hermiteTable = arma::mat(n + 1, nombreValeur, arma::fill::zeros);
    hermiteTable.row(0).ones();

    if (n == 0)
    {
        return;
    }

    // Base recurrence row for n = 1
    hermiteTable.row(1) = 2.0 * zRow;

    // Fill higher orders via the standard Hermite recurrence H_n = 2 x H_{n-1} - 2(n-1) H_{n-2}
    for (int ligne = 2; ligne <= n; ++ligne)
    {
        hermiteTable.row(ligne) = 2.0 * (zRow % hermiteTable.row(ligne - 1)) - 2.0 * (ligne - 1) * hermiteTable.row(ligne - 2);
    }
}


void Poly::calcLaguerre(int mCount, int nCount, const arma::vec& zVals)
{
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

    


    if (mCount <= 0 || nCount <= 0)
    {
        laguerreTable.reset();
        return;
    }

    // Allocate the cube: rows=m alphas, cols=n orders, slices=x samples
    laguerreTable = arma::cube(mCount, nCount, nombreValeur, arma::fill::zeros);

    for (int alpha = 0; alpha < mCount; ++alpha)
    {
        const double alpha_d = static_cast<double>(alpha);

        for (arma::uword idx = 0; idx < nombreValeur; ++idx)
        {
            const double x = zVals(idx);

            // L_0^alpha(x) = 1 for all x
            laguerreTable(alpha, 0, idx) = 1.0;

            if (nCount == 1)
            {
                continue;
            }

            // L_1^alpha(x) = 1 + alpha - x
            laguerreTable(alpha, 1, idx) = 1.0 + alpha_d - x;

            // Higher orders use the generalized Laguerre recurrence
            for (int order = 2; order < nCount; ++order)
            {
                const double order_d = static_cast<double>(order);
                const double coeffA = 2.0 * order_d - 1.0 + alpha_d - x;
                const double coeffB = order_d - 1.0 + alpha_d;
                laguerreTable(alpha, order, idx) = (coeffA * laguerreTable(alpha, order - 1, idx) - coeffB * laguerreTable(alpha, order - 2, idx)) / order_d;
            }
        }
    }
}




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
    arma::mat result(sliceNMoins1.n_rows, sliceNMoins1.n_cols, arma::fill::zeros);
    
    result = sliceNMoins1;

    return result;
}

arma::vec Poly::getHermiteRow(int order) const
{
    if (hermiteTable.is_empty())
    {
        throw std::logic_error("Hermite table is empty; call calcHermite first");
    }

    return hermiteTable.row(order).t();
}

arma::vec Poly::hermite(int order) const
{
    return getHermiteRow(order);
}

arma::vec Poly::laguerre(int mOrder, int nOrder) const
{
    // Extract the tube (all x values) for the requested (m,n) pair
    return laguerreTable.tube(mOrder, nOrder);
}
