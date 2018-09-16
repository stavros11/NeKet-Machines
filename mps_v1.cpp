#include <C:/cygwin64/home/Stavros/eigen/Eigen/Dense>
#include <complex>
//#include <fstream>
#include <vector>
#include <random>
#include <iostream> // temporarily for tests
#include "lookup.hpp"
//#include "Utils/random_utils.hpp"

void RandomGaussian(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &par, double sigma) {
  std::default_random_engine generator(1);
  std::normal_distribution<double> distribution(0, sigma);
  for (int i = 0; i < par.size(); i++) {
    par(i) =
        std::complex<double>(distribution(generator), distribution(generator));
  }
};

/**
  Abstract class for Machines.
  This class prototypes the methods needed
  by a class satisfying the Machine concept.
*/

template <typename T>
class MPS {
    using VectorType = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    using MatrixType = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;

    // Number of sites
    int N = 6;
    // Physical dimension
    int d = 2;
    // Bond dimension
    int D = 3;
    // Number of variational parameters
    int npar_ = N*d*D*D;

    // MPS Matrices
    std::vector<MatrixType> W;
    // Local lookup matrices
    std::vector<MatrixType> loc_lt;
    
 public:
  using StateType = T;
  using LookupType = Lookup<T>;

  // Auxiliary function that defines the matrices
  void Init() {
      MatrixType x=MatrixType::Zero(D, D);
      for (int i=0; i<N*d; i++){
          W.push_back(x);
      }
  };

  /**
  Member function returning the number of variational parameters.
  @return Number of variational parameters in the Machine.
  OK */ 
  int Npar() {return npar_;};

  /**
  Member function returning the current set of parameters in the machine.
  @return Current set of variational parameters of the machine.
  OK */
  VectorType GetParameters() {
    int k = 0;
    VectorType pars(npar_);

    for (int p=0; p<N*d; p++){
        for (int i=0; i<D; i++){
            for (int j=0; j<D; j++){
                pars(k) = W[p](i,j);
                k++;
            }
        }
    }
    return pars;
  };

  /**
  Member function setting the current set of parameters in the machine.
  OK*/
  void SetParameters(const VectorType &pars) {
    int k = 0;

    for (int p=0; p<N*d; p++){
        for (int i=0; i<D; i++){
            for (int j=0; j<D; j++){
                W[p](i,j) = pars(k);
                k++;
            }
        }
    }
  };

  /**
  Member function providing a random initialization of the parameters.
  @param seed is the seed of the random number generator.
  @param sigma is the variance of the gaussian.
  OK*/
 void InitRandomPars(int seed, double sigma) {
    VectorType pars(npar_);

    ::RandomGaussian(pars, sigma);
    SetParameters(pars);
  };

  /**
  Member function returning the number of visible units.
  @return Number of visible units in the Machine.
  OK*/
  int Nvisible() {return N;};

  /**
  Member function initializing the look-up tables.
  If needed, a Machine can make use of look-up tables
  to speed up some critical functions. For example,
  to speed up the calculation of wave-function ratios.
  The state of a look-up table depends on the visible units.
  This function should initialize the look-up tables
  making sure that memory in the table is also being allocated.
  @param v a constant reference to the visible configuration.
  @param lt a reference to the look-up table to be initialized.
  */
  void InitLookup(const Eigen::VectorXd &v, LookupType &lt) {
      int site;
      // Initializes local lookup too!

      // First (left) site
      _InitLookup_check(lt, 0);
      lt.M(0) = W[0 + (int)v(0)];
      loc_lt.push_back(W[0 + (int)v(0)]);

      // Last (right) site
      _InitLookup_check(lt, 1);
      lt.M(1) = W[N - 1 + (int)v(N-1)];
      loc_lt.push_back(W[N - 1 + (int)v(N-1)]);

      // Rest sites
      for (int i=2; i<2*N; i+=2) {
        _InitLookup_check(lt, i);
        site = i / 2;
        lt.M(i) = lt.M(i-2) * W[site + (int)v(site)];
        loc_lt.push_back(lt.M(i));
        
        _InitLookup_check(lt, i+1);
        site = N - 1 - site;
        lt.M(i+1) = W[site + (int)v(site)] * lt.M(i-1);
        loc_lt.push_back(lt.M(i+1));
      }
  };

  // Auxiliary function
  void _InitLookup_check(LookupType &lt, int i){
      if (lt.MatrixSize() == i){
          lt.AddMatrix(D, D);
      }
      else {
          lt.M(i).resize(D, D);
      }
  };

  /**
  Member function computing the logarithm of the wave function for a given
  visible vector. Given the current set of parameters, this function should
  comput the value of the logarithm of the wave function from scratch.
  @param v a constant reference to a visible configuration.
  @return Logarithm of the wave function.
  OK*/
 T LogVal(const Eigen::VectorXd &v) {
      MatrixType p = W[0 + (int)v(0)];
      for (int i=1; i<N; i++) { p *= W[i + (int)v(i)]; }
      return p.trace(); // A log is needed here
  };

  /**
  Member function computing the logarithm of the wave function for a given
  visible vector. Given the current set of parameters, this function should
  comput the value of the logarithm of the wave function using the information
  provided in the look-up table, to speed up the computation.
  @param v a constant reference to a visible configuration.
  @param lt a constant eference to the look-up table.
  @return Logarithm of the wave function.
  */
  T LogVal(const Eigen::VectorXd &v, const LookupType &lt) { return lt.M(2*N-2).trace(); }; // A log is needed here

  /**
  Member function computing the difference between the logarithm of the
  wave-function computed at different values of the visible units (v, and a set
  of v').
  @param v a constant reference to the current visible configuration.
  @param tochange a constant reference to a vector containing the indeces of the
  units to be modified.
  @param newconf a constant reference to a vector containing the new values of
  the visible units: here for each v', newconf(i)=v'(tochange(i)), where v' is
  the new visible state.
  @return A vector containing, for each v', log(Psi(v')) - log(Psi(v))
  */
  VectorType LogValDiff(
      const Eigen::VectorXd &v, const std::vector<std::vector<int>> &tochange,
      const std::vector<std::vector<double>> &newconf) {
    // This function has some mistakes! It assumes that newconf has the whole state
    // but probably it only has the changed positions
    const std::size_t nconn = tochange.size();
    int site = 0;
    std::size_t nchange;
    VectorType logvaldiffs(nconn);
    MatrixType current_prod = W[(int)v(0)];
    MatrixType new_prods(D, D);

    for (std::size_t k=0; k<nconn; k++) {
        nchange = tochange[k].size();

        if (tochange[k][0] == 0) {
            new_prods = W[(int)newconf[k][0]];
        }
        else {
            new_prods = W[(int)v(0)];
            for (site=1; site<tochange[k][0]; site++) {
                new_prods *= W[site + (int)v(site)];
                current_prod *= W[site + (int)v(site)];
            }
            site = tochange[k][0];
            new_prods *= W[site + (int)newconf[k][0]];
            current_prod *= W[site + (int)v(site)];
        }

        for (std::size_t i=1; i<nchange; i++){
            for (site=tochange[k][i-1]+1; site<tochange[k][i]; site++) {
                new_prods *= W[site + (int)v(site)];
                current_prod *= W[site + (int)v(site)];
            }
            site = tochange[k][i];
            new_prods *= W[site + (int)newconf[k][i]];
            current_prod *= W[site + (int)v(site)];
        }

        for (site=tochange[k][nchange-1]+1; site<N; site++) {
            new_prods *= W[site + (int)v(site)];
            current_prod *= W[site + (int)v(site)];
        }

        logvaldiffs(k) = new_prods.trace() / current_prod.trace(); // A log is needed here
    }

    return logvaldiffs;
  };

  /**
  Member function computing the difference between the logarithm of the
  wave-function computed at different values of the visible units (v, and a
  single v'). This version uses the look-up tables to speed-up the calculation.
  @param v a constant reference to the current visible configuration.
  @param tochange a constant reference to a vector containing the indeces of the
  units to be modified.
  @param newconf a constant reference to a vector containing the new values of
  the visible units: here newconf(i)=v'(tochange(i)), where v' is the new
  visible state.
  @param lt a constant eference to the look-up table.
  @return The value of log(Psi(v')) - log(Psi(v))
  */
  T LogValDiff(const Eigen::VectorXd &v, const std::vector<int> &toflip,
               const std::vector<double> &newconf,
               const LookupType &lt) {
    // Assumes that vector toflip is in ascending order

    const std::size_t nconn = toflip.size();
    MatrixType new_prod;

    if (toflip[0] == 0) {
        new_prod = W[(int)newconf[0]];
    }
    else {
        new_prod = lt.M(2 * (toflip[0] - 1)) * W[toflip[0] + (int)newconf[0]];
    }

    for (std::size_t k=1; k<nconn; k++) {
        for (int site=toflip[k-1]+1; site<toflip[k]; site++) {
            new_prod *= W[site + (int)v(site)];
        }
        new_prod *= W[toflip[k] + (int)newconf[k]];
    }

    if (toflip[nconn-1] == N-1) {
        return new_prod.trace() / LogVal(v, lt); // A log is needed here
    }
    else {
        return (new_prod * lt.M(2*(N - toflip[nconn-1])-3)).trace() / LogVal(v, lt); // A log is needed here
    }
 };

 /**
  Member function computing the derivative of the logarithm of the wave function
  for a given visible vector.
  @param v a constant reference to a visible configuration.
  @return Derivatives of the logarithm of the wave function with respect to the
  set of parameters.
  */
  VectorType DerLog(const Eigen::VectorXd &v) {
    int k = 0, Dsq = D*D;
    VectorType der(npar_);

    for (int spin=0; spin<d; spin++) {
        if ((int)v(0) != spin) {
            der.segment(k, Dsq) = VectorType::Zero(Dsq);
        }
        else {
            der.segment(k, Dsq) = Eigen::Map<VectorType>((loc_lt[2*N - 3]).transpose().data(), Dsq);
        }
        k += Dsq;
    }

    for (int site=1; site<N-1; site++) {
        for (int spin=0; spin<d; spin++) {
            if ((int)v(site) != spin) {
                der.segment(k, Dsq) = VectorType::Zero(Dsq);
            }
            else {
                der.segment(k, Dsq) = Eigen::Map<VectorType>((loc_lt[2*(N - site) - 3] * loc_lt[2 * site - 2]).transpose().data(), Dsq);
            }
            k += Dsq;
        }
    }

    for (int spin=0; spin<d; spin++) {
        if ((int)v(N-1) != spin) {
            der.segment(k, Dsq) = VectorType::Zero(Dsq);
        }
        else {
            der.segment(k, Dsq) = Eigen::Map<VectorType>((loc_lt[2*N - 4]).transpose().data(), Dsq);
        }
        k += Dsq;
    }

      return der;
  };

};

int main() {
    using namespace std;
    double sigma=0.1;
    MPS<complex<double>> x;
    Lookup<complex<double>> ltable;
    Eigen::VectorXd chain(x.Nvisible()), chain2(x.Nvisible());

    x.Init();
    x.InitRandomPars(1, sigma);
    cout << "Machine Initialized with " << x.Nvisible() << " units and " << x.Npar() << " parameters." << endl;

    chain.setZero(x.Nvisible()); chain2.setZero(x.Nvisible());
    chain(0)=1; chain(4)=1; chain(5)=1;
    chain2(0)=1; chain2(3)=1; chain2(5)=1;
    cout << "Chains updated." << endl;

    x.InitLookup(chain, ltable);
    cout << "Look up table initialized with " << ltable.MatrixSize() << " matrices." << endl;
    
    // Test LogVals
    complex<double> psi = x.LogVal(chain), psi2 = x.LogVal(chain2);

    cout << endl;
    cout << "Psi (full calculation) = " << psi << endl;
    cout << "Psi (look up) = " << x.LogVal(chain, ltable) << endl;
    cout << "Psi (look up right) =" << ltable.M(2*x.Nvisible()-1).trace() << endl;
    cout << endl;

    // Test DerLog
    //cout << "DerLog = " << x.DerLog(chain) << endl;

    // Test LogValDiff
    vector<int> change_ind;
    vector<double> upd_conf;
    vector<vector<int>> change_ind_v;
    vector<vector<double>> upd_conf_v;

    change_ind.push_back(3); upd_conf.push_back(1);
    change_ind.push_back(4); upd_conf.push_back(0);

    change_ind_v.push_back(change_ind); upd_conf_v.push_back(upd_conf);

    cout << "LogValDiff (full calculation) = " << x.LogValDiff(chain, change_ind_v, upd_conf_v) << endl;
    cout << "LogValDiff (with look up) = " << x.LogValDiff(chain, change_ind, upd_conf, ltable) << endl;
    cout << "Ratio from LogVal calculation = " << psi2 / psi << endl;

    return 0;
}