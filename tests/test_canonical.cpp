#include <iostream>
#include <C:/cygwin64/home/Stavros/eigen/Eigen/Dense>
#include <random>
#include <complex>
#include <fstream>

using VectorType = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1>;
using MatrixType = Eigen::Matrix<std::complex<double>, Eigen::Dynamic, Eigen::Dynamic>;

void RandomGaussian(Eigen::Matrix<std::complex<double>, Eigen::Dynamic, 1> &par, double sigma) {
  std::default_random_engine generator(1);
  std::normal_distribution<double> distribution(0, sigma);
  for (int i = 0; i < par.size(); i++) {
    par(i) =
        std::complex<double>(distribution(generator), distribution(generator));
  }
};

class MPS {
    // Number of sites
	int N_=6;
	// Physical dimension
	int d_=2;
	// Bond dimension (stored in vector of size N_ + 1)
	std::vector<int> D_;
	//Bond dimension given by user
	int Duser_=4;
	// Number of variational parameters
	int npar_;

    // MPS Matrices (stored as [N * d, D, D]
	std::vector<MatrixType> W_;

  public:
    void Initialize(double sigma) {
            // Initialize D vector
            D_.push_back(1);
            D_.push_back(d_);
			for (int i = 2; i < N_-1; i++) {
				D_.push_back(Duser_);
			}
            D_.push_back(d_);
			D_.push_back(1);

            // Initialize npar_
            npar_ = 0;
			for (int i = 0; i < N_; i++) {
				npar_ += D_[i] * D_[i + 1];
			}
			npar_ *= d_;

			VectorType pars(npar_);
			::RandomGaussian(pars, sigma);
			SetParameters(pars);
		};

    void SetParameters(const VectorType &pars) {
			int k = 0;

			for (int site = 0; site < N_; site++) {
				for (int spin = 0; spin < d_; spin++) {
                    W_.push_back(MatrixType::Zero(D_[site], D_[site+1]));
					for (int i = 0; i < D_[site]; i++) {
						for (int j = 0; j < D_[site + 1]; j++) {
							W_[site * d_ + spin](i, j) = pars(k);
							k++;
						}
					}
				}
			}
		};
    
    int Nsites() {
        return N_;
    }

    int physdim() {
        return d_;
    }

    std::vector<int> dims() {
        return D_;
    }

    std::complex<double> Wval(int site, int spin, int row, int col) {
        return W_[site * d_ + spin](row, col);
    }

   		 // #################################### //
		// ### Functions for canonical form ### //
		// #################################### //
		void normalize2canonical() {
			MatrixType SdotV;
            Eigen::JacobiSVD<MatrixType> svd;
			double last_norm=0.0;

			// Do SVD for site 0
            svd = JSVD(list2matLeft());
			// Update W for site 0
			mat2list(0, svd.matrixU());
            
			// Repeat for the rest sites
			for (int site=1; site<N_-1; site++) {
				SdotV = svd.singularValues().asDiagonal() * svd.matrixV().adjoint();
				svd = JSVD(list2mat(site, SdotV));
				mat2list(site, svd.matrixU());
			}
			
			// Normalize final state
			SdotV = svd.singularValues().asDiagonal() * svd.matrixV().adjoint();
			for (int i=d_*(N_ -1); i<d_*N_; i++) {
				W_[i] = SdotV * W_[i];
				last_norm += std::real((W_[i].conjugate().cwiseProduct(W_[i]).sum()));
			}
			for (int i=d_*(N_ -1); i<d_*N_; i++) {
				W_[i] *= 1.0 / std::sqrt(last_norm);
			}
		};

        inline Eigen::JacobiSVD<MatrixType> JSVD(const MatrixType x) {
            using namespace Eigen;
            JacobiSVD<MatrixType> svd(x, ComputeThinU | ComputeThinV);
            return svd;
        };

		inline MatrixType list2matLeft() {
			MatrixType mat(d_ * D_[0], D_[1]);
			for (int i=0; i<d_; i++) {
				mat.block(i * D_[0], 0, D_[0], D_[1]) = W_[i];
			}
			return mat;
		}

		inline MatrixType list2mat(const int site, const MatrixType SdotV) {
			MatrixType mat(d_ * D_[site], D_[site+1]);
			for (int i=0; i<d_; i++) {
				mat.block(i * D_[site], 0, D_[site], D_[site+1]) = SdotV * W_[site * d_ + i];
			}
			return mat;
		}

		inline void mat2list(const int site, const MatrixType mat) {
			for (int i=0; i<d_; i++) {
				W_[site * d_ + i] = mat.block(i * D_[site], 0, D_[site], D_[site+1]);
			}
		}
		// ############################################ //
		// ### End of functions for canonical form ### //
		// ########################################## //
};

int main() {
    using namespace std;
	int ver;
	double sigma;
    MPS x;

	cin >> sigma;
	x.Initialize(sigma);

    std::vector<int> dims = x.dims();
    cout << "Bond dimension vector: ";
    for (int i=0; i<x.Nsites()+1; i++) {
        cout << dims[i] << " ";
    }
    cout << endl;

	cin >> ver;

    // Save W to files
    fstream myfile;
    myfile.open("Winit" + to_string(ver) + ".txt", ios::out|ios::binary);
    for (int site=0; site<x.Nsites(); site++) {
        for (int spin=0; spin<x.physdim(); spin++) {
            for (int i=0; i<dims[site]; i++) {
                for (int j=0; j<dims[site+1]; j++) {
                    myfile << x.Wval(site, spin, i, j) << " " << endl << endl;
                }
            }
        }
    }
    myfile.close();

    // Save W to files after svd
    x.normalize2canonical();
    myfile.open("Wsvd" + to_string(ver) + ".txt", ios::out|ios::binary);
    for (int site=0; site<x.Nsites(); site++) {
        for (int spin=0; spin<x.physdim(); spin++) {
            for (int i=0; i<dims[site]; i++) {
                for (int j=0; j<dims[site+1]; j++) {
                    myfile << x.Wval(site, spin, i, j) << " " << endl << endl;
                }
            }
        }
    }
    myfile.close();

    return 0;
};