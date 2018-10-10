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
#include "mps_periodic_general.hpp"

#ifndef NETKET_MPS_PERIODIC_HPP
#define NETKET_MPS_PERIODIC_HPP

namespace netket {

template <typename T>
class MPSPeriodic
    : public MPSPeriodicGeneral<T, typename AbstractMachine<T>::MatrixType,
                                MPSPeriodic<T>> {
  using MatrixType = typename AbstractMachine<T>::MatrixType;
  using VectorType = typename AbstractMachine<T>::VectorType;

  int D_;

 public:
  using LookupType = Lookup<T>;

  // constructor as a machine
  explicit MPSPeriodic(const Hilbert &hilbert, const json &pars)
      : MPSPeriodicGeneral(hilbert, pars) {}

  int flattened_dim(const int &D) {
    D_ = D;
    return D * D;
  }

  inline MatrixType product(const MatrixType &a,
                            const MatrixType &b) const override {
    return a * b;
  }

  inline MatrixType empty_initialization(const bool &ones) const {
    if (ones) {
      return MatrixType::Identity(D_, D_);
    }
    MatrixType mat(D_, D_);
    return mat;
  }

  inline void print_creation_message(const int &N) const {
    InfoMessage() << "Periodic MPS machine with " << N << " sites created"
                  << std::endl;
  }

  inline VectorType X2Vec(const MatrixType &X) const override {
    VectorType vec(D_ * D_);
    int k = 0;
    for (int i = 0; i < D_; i++) {
      for (int j = 0; j < D_; j++) {
        vec(k) = X(i, j);
        k++;
      }
    }
    return vec;
  }

  inline void Vec2X(MatrixType &X, const VectorType &vec,
                    const bool ident) const override {
    int k = 0;
    for (int i = 0; i < D_; i++) {
      for (int j = 0; j < D_; j++) {
        X(i, j) = vec(k);
        if (ident and i == j) {
          X(i, j) += T(1, 0);
        }
        k++;
      }
    }
  }

  inline void InitLookup_check(LookupType &lt, const int i) const override {
    if (lt.MatrixSize() == i) {
      lt.AddMatrix(D_, D_);
    } else {
      lt.M(i).resize(D_, D_);
    }
  }

  inline MatrixType *PltP(LookupType &lt, const int i) const override {
    return &(lt.M(i));
  }

  inline MatrixType ltP(const LookupType &lt, const int i) const override {
    return lt.M(i);
  }
};

}  // namespace netket

#endif
