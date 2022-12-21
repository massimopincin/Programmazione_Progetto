#include <iostream>
#include <exception>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
#include <boost/numeric/ublas/triangular.hpp>
#include <boost/numeric/ublas/lu.hpp>
#include <boost/qvm/mat_operations.hpp>

using namespace boost::numeric::ublas;



enum class SolutionType
{
    BOUNDED,  // trovata una soluzione ottima
    UNBOUNDED // l'insieme delle soluzioni non è limitato superiormente
};

/**
 * Funzione per ottenere una matrice diagonale da un vettore
 *
 * La matrice diagonale avrà in posizione (i, i) il valore vector(i)
 *
 * @param vector è il vettore i cui valori devono formare la diagonale della matrice
 * @return è la matrica diagonale
 */
matrix<double> diagonale(vector<double> &vector);

/**
 * Funzione per ottenere il quadrato dell'inversa di una matrice
 *
 * @param input è la matrice su cui si vuole compiere l'inversione
 * @param inverse è la matrice inversa
 * @return è una condizione di verità che è vera se la matrice è invertibile
 */
bool invertMatrix(const matrix<double> &input, matrix<double> &inverse);

/**
 * Versione del algoritmo di Karmakar con l'utilizzo delle trasformazioni affini
 * Non è un algoritmo di complessità polinomiale
 *
 * @param A matrice di Ax <= b
 * @param b
 * @param c vettore di cTx (cT è il vettore riga di c)
 * @param x0
 * @param repetitions è il numero di iterazioni consentite per ottenere il risultato
 * @return è il vettore massimo
 */

SolutionType affineScaling(
    matrix<double> &A,
    vector<double> &b,
    vector<double> &c,
    vector<double> &x0,
    unsigned int &repetitions);


matrix<double> diagonale(vector<double> &vector)
{
    matrix<double> m(vector.size(), vector.size());

    for (size_t i{0}; i < vector.size(); ++i){
        for (size_t j{0}; j < vector.size(); ++j){
            m(i,j) = 0;
        }
    }
    for (size_t i{0}; i < vector.size(); ++i)
    {
        m(i, i) = vector(i);
    }

    return m;
}

bool invertMatrix(const matrix<double> &input, matrix<double> &inverse)
{
    typedef permutation_matrix<std::size_t> pmatrix;
    matrix<double> A(input);             // copia di lavoro della matrice input
    pmatrix pm(A.size1());               // crea una matrice di permutazione per la fattorizzazione LU
    const int res = lu_factorize(A, pm); // esegue la fattorizzazione LU
    if (res != 0)
    {
        return false;
    }
    inverse.assign(identity_matrix<double>(A.size1())); // crea la matrice identità
    lu_substitute(A, pm, inverse);                      // backsubstitute to get the inverse
    return true;
}

SolutionType affineScaling(matrix<double> &A, vector<double> &b, vector<double> &c, vector<double> &x0, unsigned int &repetitions)
{
    unsigned int k{0};
    vector<double> v[repetitions];
    vector<double> x[repetitions];
    matrix<double> Dv(b.size(), b.size());
    matrix<double> Dv_inversa(b.size(), b.size());
    matrix<double> mInv(b.size(), b.size());
    vector<double> hx(b.size());
    vector<double> hv(b.size());
    // double alpha;
    double min{0};
    double tmp;
    // double gamma{1};

    x[0] = x0;

    while (k < repetitions)
    {
        
        std::cout << "\n\nSTARTING ITERATION n. " << k << std::endl;

        v[k] = b - prod(A, x[k]);

        

        Dv = diagonale(v[k]);

        std::cout << "x[k]: " << x[k] << " Dv: " << Dv << std::endl;

        if (!invertMatrix(Dv, Dv_inversa))
        {
            std::cout << "Matrice Dv non invertibile --> Dv: " << Dv << std::endl;
            throw("Matrice Dv non invertibile");
        }

        std::cout << "Dv_inversa: " << Dv_inversa << std::endl;

        Dv = prod(Dv_inversa, Dv_inversa);
        Dv = prod(trans(A), Dv);
        Dv = prod(Dv, A);
        if (!invertMatrix(Dv, Dv_inversa))
        {
            std::cout << "Matrice A^tDv^{-2}A non invertibile --> A^tDv^{-2}A: " << Dv << std::endl;
            throw("Matrice A^tDv^{-2}A non invertibile");
        }
        hx = prod(Dv_inversa, c);
        hv = prod(-A, hx);
        for (unsigned int i{0}; i < hv.size(); ++i)
        {
            if (hv(i) > 0)
            {
                return SolutionType::UNBOUNDED;
            }
        }
        for (unsigned int i{0}; i < b.size() - 1; ++i)
        {
            std::cout << "hv: " << hv << std::endl;
            if (hv(i) < 0)
            {

                    min = -v[k](i) / hv(i);
                    std::cout << "MIN: " << min << std::endl;

                tmp = -v[k](i + 1) / hv(i + 1);
                if (min > tmp)
                {
                    min = tmp;
                }
            }
        }


        x[k + 1] = x[k] + min * hx;

        ++k;
    }
    double res{0};
    for (unsigned int i{0}; i < c.size(); ++i) {
        res += c(i)*x[repetitions - 1](i);
        std::cout << "\ni: " << i << "  res: " << res << "  c: " << c << "  x[repetitions-1](i): " << x[repetitions-1](i);
    }
    std::cout << "\n\nIl massimo è: " << res << std::endl;

    return SolutionType::BOUNDED;
}

int main()
{
    using namespace boost::numeric::ublas;

    matrix<double> A{identity_matrix<double>(2)};
    vector<double> b(2, 1);
    b[0] = 2;
    vector<double> c(2, 2);
    vector<double> x0(2, 0.1);
    x0[0] = 1;
    unsigned int repetitions{24};


    if (affineScaling(A, b, c, x0, repetitions) == SolutionType::BOUNDED)
    {
        std::cout << "BOUNDED solution" << std::endl;
    }
    else
    {
        std::cout << "UNBOUNDED solution" << std::endl;
    }

    return 0;
}