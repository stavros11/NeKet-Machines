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

#ifndef NETKET_MPS_DIAGONAL_HPP
#define NETKET_MPS_DIAGONAL_HPP

namespace netket {

template <typename T>
class MPSDiagonal
    : public MPSPeriodicGeneral<T, typename AbstractMachine<T>::VectorType,
                                MPSDiagonal<T>> {
  using VectorType = typename AbstractMachine<T>::VectorType;

  int D_;

 public:
  using LookupType = Lookup<T>;

  // constructor as a machine
  explicit MPSDiagonal(const Hilbert &hilbert, const json &pars)
      : MPSPeriodicGeneral(hilbert, pars) {}

  int flattened_dim(const int &D) {
    D_ = D;
    return D;
  }

  inline VectorType product(const VectorType &a,
                            const VectorType &b) const override {
    return a.cwiseProduct(b);
  }

  inline VectorType empty_initialization(const bool &ones) const {
    if (ones) {
      return VectorType::Ones(D_);
    }
    VectorType mat(D_);
    return mat;
  }

  inline void print_creation_message(const int &N) const {
    InfoMessage() << "Periodic diagonal MPS machine with " << N
                  << " sites created" << std::endl;
  }

  inline VectorType X2Vec(const VectorType &X) const override { return X; }

  inline void Vec2X(VectorType &X, const VectorType &vec,
                    const bool ident) const override {
    if (ident) {
      X = vec + VectorType::Ones(D_);
    } else {
      X = vec;
    }
  }

  inline void InitLookup_check(LookupType &lt, const int i) const override {
    if (lt.VectorSize() == i) {
      lt.AddVector(D_);
    } else {
      lt.V(i).resize(D_);
    }
  }

  inline VectorType *PltP(LookupType &lt, const int i) const override {
    return &(lt.V(i));
  }

  inline VectorType ltP(const LookupType &lt, const int i) const override {
    return lt.V(i);
  }
};

}  // namespace netket

#endif
