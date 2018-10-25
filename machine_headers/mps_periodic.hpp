// Copyright 2018 The Simons Foundation, Inc. - All Rights Reserved.
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//    http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.
//
// by S. Efthymiou, October 2018

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "Lookup/lookup.hpp"
#include "Utils/all_utils.hpp"

#ifndef NETKET_MPS_PERIODIC_HPP
#define NETKET_MPS_PERIODIC_HPP

namespace netket {

template <typename T, bool diag>
class MPSPeriodic : public AbstractMachine<T> {
  using VectorType = typename AbstractMachine<T>::VectorType;
  using MatrixType = typename AbstractMachine<T>::MatrixType;

  // Number of sites
  int N_;
  // Physical dimension
  int d_;
  // Bond dimension
  int D_;
  // Second matrix dimension (D for normal MPS, 1 for diagonal)
  int Dsec_;
  // D squared
  int Dsq_;
  // Number of variational parameters
  int npar_;
  // Period of translational symmetry (has to be a divisor of N)
  int symperiod_;

  // Used for tree look up
  int Nleaves_;
  // Map from site that changes to the corresponding "leaves"
  // Shape (N_, nr of leaves for the corresponding site)
  std::vector<std::vector<int>> leaves_of_site_;
  // Contractions needed to produce each leaf
  // Shape (total leaves, 2)
  std::vector<std::vector<int>> leaf_contractions_;

  // MPS Matrices (stored as [symperiod, d, D, D] or [symperiod, d, D, 1])
  std::vector<std::vector<MatrixType>> W_;

  // Map from Hilbert states to MPS indices
  std::map<double, int> confindex_;
  // Identity Matrix
  MatrixType identity_mat_;

  const Hilbert &hilbert_;

 public:
  using StateType = T;
  using LookupType = Lookup<T>;

  // constructor as a machine
  explicit MPSPeriodic(const Hilbert &hilbert, const json &pars)
      : N_(hilbert.Size()), d_(hilbert.LocalSize()), hilbert_(hilbert) {
    from_json(pars);
  }

  inline MatrixType prod(const MatrixType &m1, const MatrixType &m2) const {
    if (diag) {
      return m1.cwiseProduct(m2);
    }

    return m1 * m2;
  }

  inline T trace(const MatrixType &m) const {
    if (diag) {
      return m.sum();
    }
    return m.trace();
  }

  inline void setparamsident(MatrixType &m, const VectorType &pars) const {
    if (diag) {
      for (int i = 0; i < D_; i++) {
        m(i, 0) = T(1, 0) + pars(i);
      }
    } else {
      for (int i = 0; i < D_; i++) {
        for (int j = 0; j < D_; j++) {
          m(i, j) = pars(i * D_ + j);
          if (i == j) {
            m(i, j) += T(1, 0);
          }
        }
      }
    }
  }

  // Auxiliary function that defines the matrices
  void Init() {
    // Initialize parameters
    std::vector<MatrixType> pushback_vec;
    if (diag) {
      Dsec_ = 1;
      identity_mat_ = MatrixType::Ones(D_, 1);
    } else {
      Dsec_ = D_;
      identity_mat_ = MatrixType::Identity(D_, D_);
    }
    MatrixType init_mat = MatrixType::Zero(D_, Dsec_);
    Dsq_ = D_ * Dsec_;
    npar_ = symperiod_ * d_ * Dsq_;

    for (int site = 0; site < symperiod_; site++) {
      W_.push_back(pushback_vec);
      for (int spin = 0; spin < d_; spin++) {
        W_[site].push_back(init_mat);
      }
    }

    // Initialize tree parameters
    InitTree();

    InfoMessage() << "Number of leaves " << Nleaves_ << std::endl;
    InfoMessage() << "Leaves of site: " << std::endl;
    for (int i = 0; i < N_; i++) {
      for (std::size_t k = 0; k < leaves_of_site_[i].size(); k++) {
        std::cout << leaves_of_site_[i][k] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;

    InfoMessage() << "Leaf contractions: " << std::endl;
    for (int l = 0; l < Nleaves_; l++) {
      for (int i = 0; i < 2; i++) {
        std::cout << leaf_contractions_[l][i] << " ";
      }
      std::cout << std::endl;
    }

    // Machine creation messages
    if (diag) {
      InfoMessage() << "Periodic diagonal MPS machine with " << N_
                    << " sites created" << std::endl;
    } else {
      InfoMessage() << "Periodic MPS machine with " << N_ << " sites created"
                    << std::endl;
    }
    InfoMessage() << "Physical dimension d = " << d_
                  << " and bond dimension D = " << D_ << std::endl;
    if (symperiod_ < N_) {
      InfoMessage() << "Translation invariance is used. Number of "
                       "variational parameters is "
                    << npar_ << " instead of " << npar_ * N_ / symperiod_
                    << std::endl;
    } else {
      InfoMessage() << "Number of variational parameters is " << npar_
                    << std::endl;
    }
    // Initialize map from Hilbert space states to MPS indices
    auto localstates = hilbert_.LocalStates();
    for (int i = 0; i < d_; i++) {
      confindex_[localstates[i]] = i;
    }
  }

  void InitTree() {
    // Initializes vectors used for tree look up tables
    // leaves_of_site_ and leaf_contractions_

    std::vector<int> two_vector(2), empty_vector, level_start;
    std::vector<std::vector<int>> above;
    int level = 1, available_ind = 0, available_level = 0, Nlevels = 0;
    bool available = false;

    level_start.push_back(0);
    level_start.push_back(N_);

    int length = level_start[level] - level_start[level - 1];
    while (length > 1 or available) {
      above.push_back(empty_vector);
      // Iterate level-1
      for (int i = 0; i < length - 1; i += 2) {
        // Construct above for level-1
        above.back().push_back(level_start[level] + i / 2);
        above.back().push_back(level_start[level] + i / 2);
        // Construct leaf_contractions_ for level
        two_vector[0] = level_start[level - 1] + i;
        two_vector[1] = level_start[level - 1] + i + 1;
        leaf_contractions_.push_back(two_vector);
      }
      if (length % 2 == 1) {
        if (available) {
          // Connect the two odd leaves
          above.back().push_back(level_start[level] + (length - 1) / 2);
          above[available_level].push_back(level_start[level] +
                                           (length - 1) / 2);

          two_vector[0] = level_start[level] - 1;
          two_vector[1] = available_ind;
          leaf_contractions_.push_back(two_vector);

          level_start.push_back(level_start.back() + length / 2 + 1);
          available = false;
        } else {
          available = true;
          available_ind = level_start[level] - 1;
          available_level = level - 1;
          level_start.push_back(level_start.back() + length / 2);
        }
      } else {
        level_start.push_back(level_start.back() + length / 2);
      }
      level++;
      length = level_start[level] - level_start[level - 1];
    }
    Nleaves_ = level_start.back();

    // Flatten above vector
    std::vector<int> flat_above;
    for (std::size_t l = 0; l < above.size(); l++) {
      for (std::size_t k = 0; k < above[l].size(); k++) {
        flat_above.push_back(above[l][k]);
      }
    }

    // Create leaves_of_site_ from above vector
    for (int i = 0; i < N_; i++) {
      std::vector<int> leaves;
      leaves.push_back(flat_above[i]);
      while (flat_above[leaves.back()] < Nleaves_ - 1) {
        leaves.push_back(flat_above[leaves.back()]);
      }
      leaves.push_back(Nleaves_ - 1);
      leaves_of_site_.push_back(leaves);
      leaves.clear();
    }
    Nleaves_ += -N_;
  }

  int Npar() const override { return npar_; }

  VectorType GetParameters() override {
    int k = 0;
    VectorType pars(npar_);

    for (int site = 0; site < symperiod_; site++) {
      for (int spin = 0; spin < d_; spin++) {
        for (int i = 0; i < D_; i++) {
          for (int j = 0; j < Dsec_; j++) {
            pars(k) = W_[site][spin](i, j);
            k++;
          }
        }
      }
    }
    return pars;
  }

  void SetParameters(const VectorType &pars) override {
    int k = 0;

    for (int site = 0; site < symperiod_; site++) {
      for (int spin = 0; spin < d_; spin++) {
        for (int i = 0; i < D_; i++) {
          for (int j = 0; j < Dsec_; j++) {
            W_[site][spin](i, j) = pars(k);
            k++;
          }
        }
      }
    }
  }

  // Auxiliary function used for setting initial random parameters and adding
  // identities in every matrix
  void SetParametersIdentity(const VectorType &pars) {
    int k = 0;
    for (int site = 0; site < symperiod_; site++) {
      for (int spin = 0; spin < d_; spin++) {
        setparamsident(W_[site][spin], pars.segment(k, Dsq_));
        k += Dsq_;
      }
    }
  }

  void InitRandomPars(int seed, double sigma) override {
    VectorType pars(npar_);

    netket::RandomGaussian(pars, seed, sigma);
    SetParametersIdentity(pars);
  }

  int Nvisible() const override { return N_; }

  void InitLookup(const Eigen::VectorXd &v, LookupType &lt) override {
    InfoMessage() << "InitLookup called!" << std::endl;

    for (int k = 0; k < Nleaves_; k++) {
      std::vector<MatrixType *> m(2);
      for (int i = 0; i < 2; i++) {
        if (leaf_contractions_[k][i] < N_) {
          m[i] = &(W_[leaf_contractions_[k][i] % symperiod_]
                     [confindex_[v(leaf_contractions_[k][i])]]);
        } else {
          m[i] = &(lt.M(leaf_contractions_[k][i] - N_));
        }
      }

      InfoMessage() << "InitLookup check " << std::endl;
      InfoMessage() << (*(m[0])).size() << std::endl;
      InfoMessage() << (*(m[1])).size() << std::endl;

      _InitLookup_check(lt, k);
      lt.M(k) = prod(*(m[0]), *(m[1]));

      InfoMessage() << "InitLookup check " << k << " again" << std::endl;
    }

    InfoMessage() << "InitLookup ended!" << std::endl;
  }

  /**
  // Old lookups
  void InitLookup(const Eigen::VectorXd &v, LookupType &lt) override {
    // We need 2 * L matrix lookups for each string, where L is the MPS length

    // First (left) site
    _InitLookup_check(lt, 0);
    lt.M(0) = W_[0][confindex_[v(0)]];

    // Last (right) site
    _InitLookup_check(lt, 1);
    lt.M(1) = W_[(N_ - 1) % symperiod_][confindex_[v(N_ - 1)]];

    // Rest sites
    for (int i = 2; i < 2 * N_; i += 2) {
      _InitLookup_check(lt, i);
      int site = i / 2;
      lt.M(i) = prod(lt.M(i - 2), W_[(site %
  symperiod_)][confindex_[v(site)]]);

      _InitLookup_check(lt, i + 1);
      site = N_ - 1 - site;
      lt.M(i + 1) =
          prod(W_[site % symperiod_][confindex_[v(site)]], lt.M(i - 1));
    }
  } */

  // Auxiliary function
  inline void _InitLookup_check(LookupType &lt, int i) {
    if (lt.MatrixSize() == i) {
      lt.AddMatrix(D_, Dsec_);
    } else {
      lt.M(i).resize(D_, Dsec_);
    }
  }

  // Auxiliary function for sorting indeces
  // (copied from stackexchange - original answer by Lukasz Wiklendt)
  inline std::vector<std::size_t> sort_indeces(const std::vector<int> &v) {
    // initialize original index locations
    std::vector<std::size_t> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    // sort indexes based on comparing values in v
    std::sort(idx.begin(), idx.end(),
              [&v](std::size_t i1, std::size_t i2) { return v[i1] < v[i2]; });
    return idx;
  }

  void UpdateLookup(const Eigen::VectorXd &v, const std::vector<int> &tochange,
                    const std::vector<double> &newconf,
                    LookupType &lt) override {
    std::size_t nchange = tochange.size();
    if (nchange <= 0) {
      return;
    }

    InfoMessage() << "UpdatedLookup called!" << std::endl;

    std::vector<std::size_t> sorted_ind = sort_indeces(tochange);
    std::vector<int> updated_leaves_index;
    std::map<int, MatrixType *> ltpM;
    int site;

    // Add tochanges in the map
    for (std::size_t k = 0; k < nchange; k++) {
      site = tochange[sorted_ind[k]];
      ltpM[site] = &(W_[site % symperiod_][confindex_[newconf[sorted_ind[k]]]]);
    }

    // Add leaves for the first site
    site = tochange[sorted_ind[0]];
    for (std::size_t l = 0; l < leaves_of_site_[site].size(); l++) {
      updated_leaves_index.push_back(leaves_of_site_[site][l] - N_);
      ltpM[leaves_of_site_[site][l]] = &(lt.M(leaves_of_site_[site][l] - N_));
      // Add matrices needed for contraction in the map
      // if they do not already exist
      for (int i = 0; i < 2; i++) {
        int leaf_nr = leaf_contractions_[leaves_of_site_[site][l] - N_][i];
        if (ltpM.count(leaf_nr) == 0) {
          if (leaf_nr < N_) {
            ltpM[leaf_nr] = &(W_[leaf_nr % symperiod_][confindex_[v(leaf_nr)]]);
          } else {
            ltpM[leaf_nr] = &(lt.M(leaf_nr - N_));
          }
        }
      }
    }

    InfoMessage() << "UpdateLookup check 1!" << std::endl;

    // Add leaves for the rest sites that flip
    for (std::size_t k = 1; k < nchange; k++) {
      site = tochange[sorted_ind[k]];
      for (std::size_t l = 0; l < leaves_of_site_[site].size(); l++) {
        std::size_t m = 0;
        while (updated_leaves_index[m] < leaves_of_site_[site][l]) {
          m++;
        }
        if (updated_leaves_index[m] != leaves_of_site_[site][l]) {
          updated_leaves_index.insert(updated_leaves_index.begin() + m,
                                      leaves_of_site_[site][l] - N_);
          ltpM[leaves_of_site_[site][l]] =
              &(lt.M(leaves_of_site_[site][l] - N_));

          // Add matrices needed for contraction in the map
          for (int i = 0; i < 2; i++) {
            int leaf_nr = leaf_contractions_[leaves_of_site_[site][l] - N_][i];
            if (ltpM.count(leaf_nr) == 0) {
              if (leaf_nr < N_) {
                ltpM[leaf_nr] =
                    &(W_[leaf_nr & symperiod_][confindex_[v(leaf_nr)]]);
              } else {
                ltpM[leaf_nr] = &(lt.M(leaf_nr - N_));
              }
            }
          }
        }
      }
    }

    InfoMessage() << "UpdateLookup check 2!" << std::endl;

    // Do updates
    for (std::size_t k = 0; k < updated_leaves_index.size(); k++) {
      lt.M(updated_leaves_index[k]) =
          prod(*ltpM[leaf_contractions_[updated_leaves_index[k]][0]],
               *ltpM[leaf_contractions_[updated_leaves_index[k]][1]]);
    }

    InfoMessage() << "UpdateLookup ended!" << std::endl;
  }

  /**
  // Old lookups
  void UpdateLookup(const Eigen::VectorXd &v, const std::vector<int>
  &tochange, const std::vector<double> &newconf, LookupType &lt) override {
    std::size_t nchange = tochange.size();
    if (nchange <= 0) {
      return;
    }
    std::vector<std::size_t> sorted_ind = sort_indeces(tochange);
    int site = tochange[sorted_ind[0]];

    // Update left (site++)
    if (site == 0) {
      lt.M(0) = W_[0][confindex_[newconf[sorted_ind[0]]]];
    } else {
      lt.M(2 * site) =
          prod(lt.M(2 * (site - 1)),
               W_[site % symperiod_][confindex_[newconf[sorted_ind[0]]]]);
    }

    for (std::size_t k = 1; k < nchange; k++) {
      for (site = tochange[sorted_ind[k - 1]] + 1;
           site < tochange[sorted_ind[k]]; site++) {
        lt.M(2 * site) = prod(lt.M(2 * (site - 1)),
                              W_[site % symperiod_][confindex_[v(site)]]);
      }
      site = tochange[sorted_ind[k]];
      lt.M(2 * site) =
          prod(lt.M(2 * (site - 1)),
               W_[site % symperiod_][confindex_[newconf[sorted_ind[k]]]]);
    }

    for (site = tochange[sorted_ind[nchange - 1]] + 1; site < N_; site++) {
      lt.M(2 * site) = prod(lt.M(2 * (site - 1)),
                            W_[site % symperiod_][confindex_[v(site)]]);
    }

    // Update right (site--)
    site = tochange[sorted_ind[nchange - 1]];
    if (site == N_ - 1) {
      lt.M(1) = W_[(N_ - 1) % symperiod_]
                  [confindex_[newconf[sorted_ind[nchange - 1]]]];
    } else {
      lt.M(2 * (N_ - site) - 1) = prod(
          W_[site % symperiod_][confindex_[newconf[sorted_ind[nchange - 1]]]],
          lt.M(2 * (N_ - site) - 3));
    }

    for (std::size_t k = 0; k < nchange - 1; k++) {
      for (site = tochange[sorted_ind[nchange - 1 - k]] - 1;
           site > tochange[sorted_ind[nchange - 2 - k]]; site--) {
        lt.M(2 * (N_ - site) - 1) =
            prod(W_[site % symperiod_][confindex_[v(site)]],
                 lt.M(2 * (N_ - site) - 3));
      }
      site = tochange[sorted_ind[nchange - 2 - k]];
      lt.M(2 * (N_ - site) - 1) =
          prod(W_[site % symperiod_]
                 [confindex_[newconf[sorted_ind[nchange - 2 - k]]]],
               lt.M(2 * (N_ - site) - 3));
    }

    for (site = tochange[sorted_ind[0]] - 1; site >= 0; site--) {
      lt.M(2 * (N_ - site) - 1) =
          prod(W_[site % symperiod_][confindex_[v(site)]],
               lt.M(2 * (N_ - site) - 3));
    }
  } */

  // Auxiliary function that calculates contractions from site1 to site2
  inline MatrixType mps_contraction(const Eigen::VectorXd &v, const int &site1,
                                    const int &site2) {
    MatrixType c = identity_mat_;
    for (int site = site1; site < site2; site++) {
      c = prod(c, W_[site % symperiod_][confindex_[v(site)]]);
    }
    return c;
  }

  T LogVal(const Eigen::VectorXd &v) override {
    return std::log(trace(mps_contraction(v, 0, N_)));
  }

  T LogVal(const Eigen::VectorXd & /* v */, const LookupType &lt) override {
    // return std::log(trace(lt.M(2 * N_ - 2)));
    return std::log(trace(lt.M(Nleaves_ - 1)));
  }

  VectorType LogValDiff(
      const Eigen::VectorXd &v, const std::vector<std::vector<int>> &tochange,
      const std::vector<std::vector<double>> &newconf) override {
    const std::size_t nconn = tochange.size();

    std::vector<std::size_t> sorted_ind;
    VectorType logvaldiffs = VectorType::Zero(nconn);
    StateType current_psi = trace(mps_contraction(v, 0, N_));
    MatrixType new_prods(D_, Dsec_);

    for (std::size_t k = 0; k < nconn; k++) {
      std::size_t nchange = tochange[k].size();
      if (nchange > 0) {
        sorted_ind = sort_indeces(tochange[k]);
        int site = tochange[k][sorted_ind[0]];

        if (site == 0) {
          new_prods = W_[0][confindex_[newconf[k][sorted_ind[0]]]];
        } else {
          new_prods = prod(
              mps_contraction(v, 0, site),
              W_[site % symperiod_][confindex_[newconf[k][sorted_ind[0]]]]);
        }

        for (std::size_t i = 1; i < nchange; i++) {
          site = tochange[k][sorted_ind[i]];
          new_prods = prod(
              new_prods,
              prod(mps_contraction(v, tochange[k][sorted_ind[i - 1]] + 1, site),
                   W_[site % symperiod_]
                     [confindex_[newconf[k][sorted_ind[i]]]]));
        }
        site = tochange[k][sorted_ind[nchange - 1]];
        if (site < N_ - 1) {
          new_prods = prod(new_prods, mps_contraction(v, site + 1, N_));
        }
        logvaldiffs(k) = std::log(trace(new_prods) / current_psi);
      }
    }
    return logvaldiffs;
  }

  T LogValDiff(const Eigen::VectorXd &v, const std::vector<int> &toflip,
               const std::vector<double> &newconf,
               const LookupType &lt) override {
    // Assumes number of levels > 1 (?)
    std::size_t nflip = toflip.size();
    if (nflip <= 0) {
      return T(0, 0);
    }

    InfoMessage() << "LogValDiff called!" << std::endl;

    MatrixType empty_matrix(D_, Dsec_);
    std::vector<std::size_t> sorted_ind = sort_indeces(toflip);
    std::vector<MatrixType> updated_leaves;
    std::vector<int> updated_leaves_index;
    std::vector<const MatrixType *> two_vector(2);
    std::vector<std::vector<const MatrixType *>> updated_leaves_contractions;
    // Map from leaf number to the corresponding matrix (changed or not)
    std::map<int, MatrixType *> ltpM;
    int site;

    // Add toflips in the map
    for (std::size_t k = 0; k < nflip; k++) {
      site = toflip[sorted_ind[k]];
      ltpM[site] = &(W_[site][confindex_[newconf[sorted_ind[k]]]]);
    }

    InfoMessage() << "LogValDiff check 1!" << std::endl;

    // Add leaves for the first site
    site = toflip[sorted_ind[0]];
    for (std::size_t l = 0; l < leaves_of_site_[site].size(); l++) {
      updated_leaves.push_back(empty_matrix);
      updated_leaves_index.push_back(leaves_of_site_[site][l] - N_);
      ltpM[leaves_of_site_[site][l]] = &(updated_leaves.back());
      // Add matrices needed for contraction in the map
      // if they do not already exist
      for (int i = 0; i < 2; i++) {
        int leaf_nr = leaf_contractions_[leaves_of_site_[site][l] - N_][i];
        if (ltpM.count(leaf_nr) == 0) {
          if (leaf_nr < N_) {
            two_vector[i] = &(W_[leaf_nr % symperiod_][confindex_[v(leaf_nr)]]);
          } else {
            two_vector[i] = &(lt.M(leaf_nr - N_));
          }
        } else {
          two_vector[i] = ltpM[leaf_nr];
        }
      }
      updated_leaves_contractions.push_back(two_vector);
    }

    InfoMessage() << "LogValDiff check 2!" << std::endl;

    // Add leaves for the rest sites that flip
    for (std::size_t k = 1; k < nflip; k++) {
      site = toflip[sorted_ind[k]];
      for (std::size_t l = 0; l < leaves_of_site_[site].size(); l++) {
        std::size_t m = 0;
        while (updated_leaves_index[m] < leaves_of_site_[site][l]) {
          m++;
        }
        if (updated_leaves_index[m] != leaves_of_site_[site][l]) {
          updated_leaves_index.insert(updated_leaves_index.begin() + m,
                                      leaves_of_site_[site][l] - N_);
          updated_leaves.insert(updated_leaves.begin() + m, empty_matrix);
          ltpM[leaves_of_site_[site][l]] = &(updated_leaves[m]);

          // Add matrices needed for contraction in the map
          for (int i = 0; i < 2; i++) {
            int leaf_nr = leaf_contractions_[leaves_of_site_[site][l] - N_][i];
            if (ltpM.count(leaf_nr) == 0) {
              if (leaf_nr < N_) {
                two_vector[i] =
                    &(W_[leaf_nr % symperiod_][confindex_[v(leaf_nr)]]);
              } else {
                two_vector[i] = &(lt.M(leaf_nr - N_));
              }
            } else {
              two_vector[i] = ltpM[leaf_nr];
            }
          }
          updated_leaves_contractions.insert(
              updated_leaves_contractions.begin() + m, two_vector);
        }
      }
    }

    InfoMessage() << "LogValDiff check 3!" << std::endl;

    // Update products
    for (std::size_t k = 0; k < updated_leaves.size(); k++) {
      updated_leaves[k] = prod(*(updated_leaves_contractions[k][0]),
                               *(updated_leaves_contractions[k][1]));
    }

    InfoMessage() << "LogValDiff ended!" << std::endl;

    return std::log(trace(updated_leaves.back()) / trace(lt.M(Nleaves_ - 1)));
  }

  /**
  // Old look up
  T LogValDiff(const Eigen::VectorXd &v, const std::vector<int> &toflip,
               const std::vector<double> &newconf,
               const LookupType &lt) override {
    const std::size_t nflip = toflip.size();
    if (nflip <= 0) {
      return T(0, 0);
    }
    MatrixType new_prod;
    std::vector<std::size_t> sorted_ind = sort_indeces(toflip);
    int site = toflip[sorted_ind[0]];

    if (site == 0) {
      new_prod = W_[0][confindex_[newconf[sorted_ind[0]]]];
    } else {
      new_prod =
          prod(lt.M(2 * (site - 1)),
               W_[site % symperiod_][confindex_[newconf[sorted_ind[0]]]]);
    }

    for (std::size_t k = 1; k < nflip; k++) {
      site = toflip[sorted_ind[k]];
      new_prod =
          prod(new_prod,
               prod(mps_contraction(v, toflip[sorted_ind[k - 1]] + 1, site),
                    W_[site %
  symperiod_][confindex_[newconf[sorted_ind[k]]]]));
    }

    site = toflip[sorted_ind[nflip - 1]];
    if (site < N_ - 1) {
      new_prod = prod(new_prod, lt.M(2 * (N_ - site) - 3));
    }

    return std::log(trace(new_prod) / trace(lt.M(2 * N_ - 2)));
  } */

  // Derivative with full calculation
  VectorType DerLog(const Eigen::VectorXd &v) override {
    MatrixType temp_product(D_, Dsec_);
    std::vector<MatrixType> left_prods, right_prods;
    VectorType der = VectorType::Zero(npar_);

    // Calculate products
    left_prods.push_back(W_[0][confindex_[v(0)]]);
    right_prods.push_back(W_[(N_ - 1) % symperiod_][confindex_[v(N_ - 1)]]);
    for (int site = 1; site < N_ - 1; site++) {
      left_prods.push_back(prod(left_prods[site - 1],
                                W_[site % symperiod_][confindex_[v(site)]]));
      right_prods.push_back(
          prod(W_[(N_ - 1 - site) % symperiod_][confindex_[v(N_ - 1 - site)]],
               right_prods[site - 1]));
    }
    left_prods.push_back(prod(
        left_prods[N_ - 2], W_[(N_ - 1) % symperiod_][confindex_[v(N_ - 1)]]));
    right_prods.push_back(prod(W_[0][confindex_[v(0)]], right_prods[N_ - 2]));

    der.segment(confindex_[v(0)] * Dsq_, Dsq_) +=
        Eigen::Map<VectorType>(right_prods[N_ - 2].transpose().data(), Dsq_);
    for (int site = 1; site < N_ - 1; site++) {
      temp_product = prod(right_prods[N_ - site - 2], left_prods[site - 1]);
      der.segment((d_ * (site % symperiod_) + confindex_[v(site)]) * Dsq_,
                  Dsq_) +=
          Eigen::Map<VectorType>(temp_product.transpose().data(), Dsq_);
    }
    der.segment((d_ * ((N_ - 1) % symperiod_) + confindex_[v(N_ - 1)]) * Dsq_,
                Dsq_) +=
        Eigen::Map<VectorType>(left_prods[N_ - 2].transpose().data(), Dsq_);

    return der / trace(left_prods[N_ - 1]);
  }

  const Hilbert &GetHilbert() const { return hilbert_; }

  // Json functions
  void to_json(json &j) const override {
    j["Machine"]["Name"] = "MPSperiodic";
    j["Machine"]["Length"] = N_;
    j["Machine"]["BondDim"] = D_;
    j["Machine"]["PhysDim"] = d_;
    j["Machine"]["Diagonal"] = diag;
    j["Machine"]["SymmetryPeriod"] = symperiod_;
    for (int i = 0; i < symperiod_; i++) {
      for (int k = 0; k < d_; k++) {
        j["Machine"]["W" + std::to_string(d_ * i + k)] = W_[i][k];
      }
    }
  }

  void from_json(const json &pars) override {
    if (pars.at("Machine").at("Name") != "MPSperiodic") {
      throw InvalidInputError("Error while constructing MPS from Json input");
    }

    if (FieldExists(pars["Machine"], "Length")) {
      N_ = pars["Machine"]["Length"];
    }
    if (N_ != hilbert_.Size()) {
      throw InvalidInputError(
          "Number of spins is incompatible with given Hilbert space");
    }

    if (FieldExists(pars["Machine"], "PhysDim")) {
      d_ = pars["Machine"]["PhysDim"];
    }
    if (d_ != hilbert_.LocalSize()) {
      throw InvalidInputError(
          "Number of spins is incompatible with given Hilbert space");
    }

    if (FieldExists(pars["Machine"], "BondDim")) {
      D_ = pars["Machine"]["BondDim"];
    } else {
      throw InvalidInputError("Unspecified bond dimension");
    }

    if (FieldExists(pars["Machine"], "SymmetryPeriod")) {
      symperiod_ = pars["Machine"]["SymmetryPeriod"];
    } else {
      // Default is symperiod = N, resp. no translational symmetry
      symperiod_ = N_;
    }

    Init();

    // Loading parameters, if defined in the input
    from_jsonWeights(pars["Machine"]);
  }

  inline void from_jsonWeights(const json &pars) {
    for (int i = 0; i < symperiod_; i++) {
      for (int k = 0; k < d_; k++) {
        if (FieldExists(pars, "W" + std::to_string(d_ * i + k))) {
          W_[i][k] = pars["W" + std::to_string(d_ * i + k)];
        }
      }
    }
  }
};

}  // namespace netket

#endif
