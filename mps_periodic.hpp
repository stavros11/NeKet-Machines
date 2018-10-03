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

#include <Eigen/Dense>
#include <iostream>
#include <vector>
#include "Lookup/lookup.hpp"
#include "Utils/all_utils.hpp"

#ifndef NETKET_MPS_TRANSLATION_HPP
#define NETKET_MPS_TRANSLATION_HPP

namespace netket {

	template <typename T>
	class MPSPeriodic : public AbstractMachine<T> {
		using VectorType = typename AbstractMachine<T>::VectorType;
		using MatrixType = typename AbstractMachine<T>::MatrixType;
		using Ptype = std::unique_ptr<AbstractMachine<T>>;

		// Number of sites
		int N_;
		// Physical dimension
		int d_;
		// Bond dimension
		int D_;
		// Period of translational symmetry (has to be a divisor of N)
		int symperiod_;
		// Use diagonal matrices
		bool diagonal_;

		Ptype m_;

		const Hilbert &hilbert_;

	public:
		using StateType = T;
		using LookupType = Lookup<T>;

		// constructor as a machine
		explicit MPSPeriodic(const Hilbert &hilbert, const json &pars)
			: N_(hilbert.Size()), hilbert_(hilbert), d_(hilbert.LocalSize()) {
			from_json(pars);
		};

		// Auxiliary function that defines the matrices
		void Init() {
			if (diagonal_) {
				m_ = new MPSDiagonal<T>(hilbert_, N_, D_, symperiod_);
				InfoMessage() << "Periodic diagonal MPS machine with " << N_ << " sites created" << std::endl;
			}
			else {
				m_ = new MPSTranslation<T>(hilbert_, N_, D_, symperiod_);
				InfoMessage() << "Periodic MPS machine with " << N_ << " sites created" << std::endl;
			}
			int npar_ = m_->Npar();

			InfoMessage() << "Physical dimension d = " << d_ << " and bond dimension D = " << D_ << std::endl;
			if (symperiod_ < N_) {
				InfoMessage() << "Translation invariance is used. Number of variational parameters is " << npar_ << " instead of " << npar_ * N_ / symperiod_ << std::endl;
			}
			else {
				InfoMessage() << "Number of variational parameters is " << npar_ << std::endl;
			}
		};

		int Npar() const override { 
			return m_->Npar(); 
		};

		VectorType GetParameters() override {
			return m_->
		};

		void SetParameters(const VectorType &pars) override {
			
		};

		// Auxiliary function used for setting initial random parameters and adding identities in every matrix
		inline void SetParametersIdentity(const VectorType &pars) override {
			
		};
		
		void InitRandomPars(int seed, double sigma) override {
			VectorType pars(npar_);

			netket::RandomGaussian(pars, seed, sigma);
			SetParametersIdentity(pars);
		};

		int Nvisible() const override {  };

		void InitLookup(const Eigen::VectorXd &v, LookupType &lt) override {
			
		};

		// Auxiliary function
		inline void _InitLookup_check(LookupType &lt, int i) {
			
		};

		// Check lookups later
		void UpdateLookup(const Eigen::VectorXd &v,
			const std::vector<int> &tochange,
			const std::vector<double> &newconf,
			LookupType &lt) override {
				
			
		};

		T LogVal(const Eigen::VectorXd &v) override {

		};

		T LogVal(const Eigen::VectorXd &v, const LookupType &lt) override {

		};

		VectorType LogValDiff(
			const Eigen::VectorXd &v, const std::vector<std::vector<int>> &tochange,
			const std::vector<std::vector<double>> &newconf) override {

		};

		T LogValDiff(const Eigen::VectorXd &v, const std::vector<int> &toflip,
			const std::vector<double> &newconf,
			const LookupType &lt) override {

		};

		// Derivative with full calculation
		VectorType DerLog(const Eigen::VectorXd &v) override {

		};

		// Json functions
		const Hilbert &GetHilbert() const { return hilbert_; };

		void to_json(json &j) const override {
		  j["Machine"]["Name"] = "MPSperiodic";
		  j["Machine"]["Nspins"] = N_;
		  j["Machine"]["BondDim"] = D_;
		  j["Machine"]["PhysDim"] = d_;
		  j["Machine"]["SymmetryPeriod"] = symperiod_;
		  //j["Machine"]["W"] = W_;
		}; 

		void from_json(const json &pars) override {
		  if (pars.at("Machine").at("Name") != "MPStranslation") {
			throw InvalidInputError("Error while constructing MPS from Json input");
		  }

		  if (FieldExists(pars["Machine"], "Nspins")) {
			N_ = pars["Machine"]["Nspins"];
		  }
		  if (N_ != hilbert_.Size()) {
			throw InvalidInputError("Number of spins is incompatible with given Hilbert space");
		  }

		  if (FieldExists(pars["Machine"], "PhysDim")) {
			  d_ = pars["Machine"]["PhysDim"];
		  }
		  if (d_ != hilbert_.LocalSize()) {
			  throw InvalidInputError("Number of spins is incompatible with given Hilbert space");
		  }

		  if (FieldExists(pars["Machine"], "BondDim")) {
			D_ = pars["Machine"]["BondDim"];
		  }
		  else {
			  throw InvalidInputError("Unspecified bond dimension");
		  }

		  if (FieldExists(pars["Machine"], "SymmetryPeriod")) {
			  symperiod_ = pars["Machine"]["SymmetryPeriod"];
		  }
		  else {
			  // Default is symperiod = N, resp. no translational symmetry - normal periodic MPS
			  symperiod_ = N_;
		  }

		  if (FieldExists(pars["Machine"], "Diagonal")) {
			  diagonal_ = pars["Machine"]["Diagonal"];
		  }
		  else {
			  diagonal_ = false;
		  }

		  Init();

		  // Loading parameters, if defined in the input
		  if (FieldExists(pars["Machine"], "W")) {
			  for (int i = 0; i < symperiod_; i++) {
				  for (int j = 0; j < d_; j++) {
					  W_[i][j] = pars["Machine"]["W"][d_ * i + j];
				  }
			  }
		  }
		};
	};

} // namespace netket

#endif
