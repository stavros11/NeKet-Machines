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
#include "abstract_mps.hpp"
#include "mps_diagonal.hpp"
#include "mps_periodic.hpp"

#ifndef NETKET_SBS_HPP
#define NETKET_SBS_HPP

namespace netket {

template <typename T>
class SBS : public AbstractMachine<T> {
  using VectorType = typename AbstractMachine<T>::VectorType;
  using MatrixType = typename AbstractMachine<T>::MatrixType;
  using Ptype = std::unique_ptr<AbstractMPS<T>>;

  // Number of sites
  int N_;
  // Physical dimension
  int d_;
  // Number of strings
  int M_;
  // Bond dimension of each string (allow different dimensions among strings)
  std::vector<int> Dstr_;
  // Length of each string (allow different lengths)
  std::vector<int> Lstr_;
  // Cumulative sum of Lstr used for lookup start indices
  std::vector<int> Lstr_cumsum_;

  // Number of total SBS variational parameters
  int npar_;

  // Map from Hilbert states to MPS indices
  std::map<double, int> confindex_;

  // Map from site numbering in v to numbering to each string
  // for each site in v we get n 2-component vectors: (string number, index
  // within string), where n is the number of strings that the site is in
  std::vector<std::vector<std::vector<int>>> site2string_;
  // Map from strings to sites
  std::vector<std::vector<int>> string2site_;

  // Vector that stores the MPS object for each string
  std::vector<Ptype> strings_;

  const Hilbert &hilbert_;

 public:
  using StateType = T;
  using LookupType = Lookup<T>;

  // constructor
  explicit SBS(const Hilbert &hilbert, const json &pars)
      : N_(hilbert.Size()), d_(hilbert.LocalSize()), hilbert_(hilbert) {
    from_json(pars);
  }

  // Auxiliary function that defines the matrices
  void Init(const int Mdiag) {
    std::vector<std::vector<int>> pushback_vec2;
    std::vector<int> two_component_vec;
    two_component_vec.push_back(0);
    two_component_vec.push_back(0);

    // Initialize map from Hilbert space states to MPS indices
    auto localstates = hilbert_.LocalStates();
    for (int i = 0; i < d_; i++) {
      confindex_[localstates[i]] = i;
    }

    // Find site2string vector by reversing string2site (defined in json)
    for (int site = 0; site < N_; site++) {
      site2string_.push_back(pushback_vec2);
      for (int i = 0; i < M_; i++) {
        for (int pos = 0; pos < Lstr_[i]; pos++) {
          if (site == string2site_[i][pos]) {
            site2string_[site].push_back(two_component_vec);
            site2string_[site].back()[0] = i;
            site2string_[site].back()[1] = pos;
          }
        }
      }
    }

    // Machine creation messages
    InfoMessage() << "SBS machine with " << N_ << " sites created" << std::endl;
    InfoMessage() << "Physical dimension d = " << d_ << std::endl;
    InfoMessage() << "Number of strings M = " << M_ << std::endl;
    InfoMessage() << "Number of diagonal strings Mdiag = " << Mdiag
                  << std::endl;
    InfoMessage() << "Bond dimensions of strings are: ";
    for (int i = 0; i < M_; i++) {
      std::cout << Dstr_[i] << " ";
    }
    std::cout << std::endl;
    for (int i = 0; i < M_; i++) {
      std::cout << Lstr_cumsum_[i] << " ";
    }
    std::cout << std::endl;
    InfoMessage() << "The number of variational parameters is " << npar_
                  << std::endl;
  };

  int Npar() const override { return npar_; };

  VectorType GetParameters() override {
    int seg_init = 0;
    VectorType pars(npar_);

    for (int i = 0; i < M_; i++) {
      pars.segment(seg_init, strings_[i]->Npar()) =
          strings_[i]->GetParameters();
      seg_init += strings_[i]->Npar();
    }
    return pars;
  };

  void SetParameters(const VectorType &pars) override {
    int seg_init = 0;

    for (int i = 0; i < M_; i++) {
      strings_[i]->SetParameters(pars.segment(seg_init, strings_[i]->Npar()));
      seg_init += strings_[i]->Npar();
    }
  };

  // Auxiliary function used for setting initial random parameters and adding
  // identities in every matrix
  inline void SetParametersIdentity(const VectorType &pars) {
    int seg_init = 0;

    for (int i = 0; i < M_; i++) {
      strings_[i]->SetParametersIdentity(
          pars.segment(seg_init, strings_[i]->Npar()));
      seg_init += strings_[i]->Npar();
    }
  }

  void InitRandomPars(int seed, double sigma) override {
    VectorType pars(npar_);

    netket::RandomGaussian(pars, seed, sigma);
    SetParametersIdentity(pars);
  };

  int Nvisible() const override { return N_; };

  void InitLookup(const Eigen::VectorXd &v, LookupType &lt) override {
    for (int i = 0; i < M_; i++) {
      strings_[i]->InitLookup(extract(v, i), lt, Lstr_cumsum_[i]);
    }
  };

  void UpdateLookup(const Eigen::VectorXd &v, const std::vector<int> &tochange,
                    const std::vector<double> &newconf,
                    LookupType &lt) override {
    if (tochange.size() > 0) {
      std::vector<std::map<int, std::vector<int>>> string_lists =
          tochange4string(tochange, newconf);

      for (auto const &ent : string_lists[0]) {
        int i = ent.first;
        strings_[i]->UpdateLookup(extract(v, i), string_lists[0][i],
                                  string_lists[1][i], lt, Lstr_cumsum_[i]);
      }
    }
  };

  // Auxiliary function that takes v and returns the visible vector for the
  // corresponding string
  inline std::vector<int> extract(const Eigen::VectorXd &v, const int string) {
    std::vector<int> x;
    for (int i = 0; i < Lstr_[string]; i++) {
      x.push_back(confindex_[v(string2site_[string][i])]);
    }
    return x;
  }

  T LogVal(const Eigen::VectorXd &v) override {
    T s = T(0, 0);
    for (int i = 0; i < M_; i++) {
      s += strings_[i]->LogVal(extract(v, i));
    }
    return s;
  };

  T LogVal(const Eigen::VectorXd &v, const LookupType &lt) override {
    T s = T(0, 0);
    for (int i = 0; i < M_; i++) {
      s += strings_[i]->LogVal(lt, Lstr_cumsum_[i]);
    }
    return s;
  };

  VectorType LogValDiff(
      const Eigen::VectorXd &v, const std::vector<std::vector<int>> &tochange,
      const std::vector<std::vector<double>> &newconf) override {
    const std::size_t nconn = tochange.size();
    if (nconn <= 0) {
      return VectorType::Zero(nconn);
    }
    VectorType logvaldiffs = VectorType::Zero(nconn);
    std::vector<std::map<int, std::vector<int>>> string_lists;

    for (std::size_t k = 0; k < nconn; k++) {
      std::size_t nchange = tochange[k].size();
      if (nchange > 0) {
        string_lists = tochange4string(tochange[k], newconf[k]);
        for (auto const &ent : string_lists[0]) {
          int i = ent.first;
          logvaldiffs(k) += strings_[i]->LogValDiff(
              extract(v, i), string_lists[0][i], string_lists[1][i]);
        }
      }
    }
    return logvaldiffs;
  };

  T LogValDiff(const Eigen::VectorXd &v, const std::vector<int> &toflip,
               const std::vector<double> &newconf,
               const LookupType &lt) override {
    // InfoMessage() << "LogValDiff lookup called" << std::endl;

    std::size_t nflip = toflip.size();
    if (nflip <= 0) {
      return T(0, 0);
    }

    T result = T(0, 0);
    std::vector<std::map<int, std::vector<int>>> string_lists =
        tochange4string(toflip, newconf);

    for (auto const &ent : string_lists[0]) {
      int i = ent.first;
      if (string_lists[0][i].size() > 1) {
        result +=
            strings_[i]->LogValDiff(extract(v, i), string_lists[0][i],
                                    string_lists[1][i], lt, Lstr_cumsum_[i]);
      } else if (string_lists[0][i].size() > 0) {
        result += strings_[i]->FastLogValDiff(
            string_lists[0][i], string_lists[1][i], lt, Lstr_cumsum_[i]);
      }
    }
    return result;
  };

  // Auxiliary function that gives the tochange and newconf maps for each string
  inline std::vector<std::map<int, std::vector<int>>> tochange4string(
      const std::vector<int> &tochange, const std::vector<double> &newconf) {
    std::vector<std::map<int, std::vector<int>>> results;
    std::map<int, std::vector<int>> string_tochange, string_newconf;

    for (std::size_t j = 0; j < tochange.size(); j++) {
      int site = tochange[j];
      // For each site in toflip find the strings it belongs. Save the
      // corresponding positions:
      for (std::size_t i = 0; i < site2string_[site].size(); i++) {
        int string_nr = site2string_[site][i][0];
        string_tochange[string_nr].push_back(site2string_[site][i][1]);
        string_newconf[string_nr].push_back(confindex_[newconf[j]]);
      }
    }
    results.push_back(string_tochange);
    results.push_back(string_newconf);
    return results;
  };

  // Derivative with full calculation
  VectorType DerLog(const Eigen::VectorXd &v) override {
    VectorType der = VectorType::Zero(npar_);
    int seg_init = 0;

    for (int i = 0; i < M_; i++) {
      der.segment(seg_init, strings_[i]->Npar()) =
          strings_[i]->DerLog(extract(v, i));
      seg_init += strings_[i]->Npar();
    }

    return der;
  };

  const Hilbert &GetHilbert() const { return hilbert_; };

  // Json functions
  void to_json(json &j) const override {
    j["Machine"]["Name"] = "SBS";
    j["Machine"]["Nsites"] = N_;
    j["Machine"]["PhysDim"] = d_;
    j["Machine"]["Strings"] = {};
    for (int i = 0; i < M_; i++) {
      strings_[i]->to_json_strings(j, string2site_[i]);
    }
  };

  void from_json(const json &pars) override {
    int Mdiag = 0;
    std::vector<int> empty_vector;

    if (pars.at("Machine").at("Name") != "SBS") {
      throw InvalidInputError("Error while constructing SBS from Json input");
    }

    if (FieldExists(pars["Machine"], "Nsites")) {
      N_ = pars["Machine"]["Nsites"];
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

    if (FieldExists(pars["Machine"], "Strings")) {
      Lstr_cumsum_.push_back(0);
      npar_ = 0;
      M_ = pars["Machine"]["Strings"].size();

      for (int string = 0; string < M_; string++) {
        json stringpars = pars["Machine"]["Strings"][string];
        bool diagonal_flag = false;
        int symperiod;

        // Assign sites to each string (string2site)
        if (FieldExists(stringpars, "SiteNumbers")) {
          string2site_.push_back(empty_vector);
          int temp_string_length = 0;
          for (auto const &j : stringpars["SiteNumbers"]) {
            string2site_.back().push_back(j);
            temp_string_length++;
          }
          Lstr_.push_back(temp_string_length);
          Lstr_cumsum_.push_back(Lstr_cumsum_.back() + 2 * temp_string_length);
        } else {
          // Default is strings that cover the whole configuration
          if (FieldExists(stringpars, "Length")) {
            if (stringpars["Length"] != N_) {
              throw InvalidInputError("Unspecified sites of the string");
            }
          }
          string2site_.push_back(empty_vector);
          for (int j = 0; j < N_; j++) {
            string2site_.back().push_back(j);
          }
          Lstr_.push_back(N_);
          Lstr_cumsum_.push_back(Lstr_cumsum_.back() + 2 * N_);
        }

        // Create string objects
        if (FieldExists(stringpars, "BondDim")) {
          Dstr_.push_back(stringpars["BondDim"]);
        } else {
          throw InvalidInputError("Unspecified bond dimension");
        }

        if (FieldExists(stringpars, "Diagonal")) {
          diagonal_flag = stringpars["Diagonal"];
        }
        if (FieldExists(stringpars, "SymmetryPeriod")) {
          symperiod = stringpars["SymmetryPeriod"];
          if (Lstr_.back() % symperiod != 0) {
            throw InvalidInputError(
                "Symmetry period is not a divisor of string length");
          }
        } else {
          symperiod = Lstr_.back();
        }

        // Initialize MPS objects and npar_ with the correct dimensions and
        // calculate npar_ String properties are defined in from_json function
        if (diagonal_flag) {
          strings_.push_back(Ptype(new MPSDiagonal<T>(
              hilbert_, Lstr_.back(), Dstr_.back(), symperiod)));
          Mdiag++;
        } else {
          strings_.push_back(Ptype(new MPSPeriodic<T>(
              hilbert_, Lstr_.back(), Dstr_.back(), symperiod)));
        }
        npar_ += strings_.back()->Npar();

        // Loading parameters, if defined in the input
        strings_.back()->from_jsonWeights(stringpars);
      }
    } else {
      throw InvalidInputError("Insufficient information to create strings");
    }

    Init(Mdiag);
  };
};

}  // namespace netket

#endif
