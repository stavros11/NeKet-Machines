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
	class MPSTranslation : public AbstractMachine<T> {
		using VectorType = typename AbstractMachine<T>::VectorType;
		using MatrixType = typename AbstractMachine<T>::MatrixType;

		// Number of sites
		int N_;
		// Physical dimension
		int d_;
		// Bond dimension
		int D_;
		// Number of variational parameters
		int npar_;
		// Period of translational symmetry (has to be a divisor of N)
		int symperiod_;

		// MPS Matrices (stored as [symperiod * d, D, D]
		std::vector<MatrixType> W_;
		
		// Map from Hilbert states to MPS indices
		std::map<double, int> confindex_;

		const Hilbert &hilbert_;

	public:
		using StateType = T;
		using LookupType = Lookup<T>;

		// constructor
		explicit MPSTranslation(const Hilbert &hilbert, const json &pars)
			: N_(hilbert.Size()), hilbert_(hilbert), d_(hilbert.LocalSize()) {
			from_json(pars);
		}

		// Auxiliary function that defines the matrices
		void Init() {
			MatrixType x = MatrixType::Zero(D_, D_);
			npar_ = symperiod_ * d_ * D_ * D_;

			for (int i = 0; i < symperiod_*d_; i++) {
				W_.push_back(x);
			}

			// Machine creation messages
			InfoMessage() << "Periodic MPS machine with " << N_ << " sites created" << std::endl;
			InfoMessage() << "Physical dimension d = " << d_ << " and bond dimension D = " << D_ << std::endl;
			InfoMessage() << "Translation invariance is used. Number of variational parameters is " << npar_ << " instead of " << npar_ * N_ / symperiod_ << std::endl;

			// Initialize map from Hilbert space states to MPS indices
			auto localstates = hilbert_.LocalStates();
			for (int i = 0; i < d_; i++) {
				confindex_[localstates[i]] = i;
			}
		};

		// Auxiliary function that computes local vector vtilde_
		//inline void ComputeVtilde(const Eigen::VectorXd &v, Eigen::VectorXd &vtilde) {
		//	vtilde = (v + Eigen::VectorXd::Ones(N_)) / 2;
		//}

		// Auxiliary function that transforms newconf to {0,1}
		//inline int ComputeNewConftilde(const double x) {
		//	return (int)((x + 1) / 2);
		//}

		int Npar() const override { return npar_; };

		VectorType GetParameters() override {
			int k = 0;
			VectorType pars(npar_);

			for (int p = 0; p < symperiod_*d_; p++) {
				for (int i = 0; i < D_; i++) {
					for (int j = 0; j < D_; j++) {
						pars(k) = W_[p](i, j);
						k++;
					}
				}
			}
			return pars;
		};

		void SetParameters(const VectorType &pars) override {
			int k = 0;

			for (int p = 0; p < symperiod_*d_; p++) {
				for (int i = 0; i < D_; i++) {
					for (int j = 0; j < D_; j++) {
						W_[p](i, j) = pars(k);
						k++;
					}
				}
			}
		};

		// Auxiliary function used for setting initial random parameters and adding identities in every matrix
		inline void SetParametersIdentity(const VectorType &pars) {
			int k = 0;

			for (int p = 0; p < symperiod_*d_; p++) {
				for (int i = 0; i < D_; i++) {
					for (int j = 0; j < D_; j++) {
						W_[p](i, j) = pars(k);
						if (i == j) {
							W_[p](i, j) += T(1, 0);
						}
						k++;
					}
				}
			}
		};
		
		void InitRandomPars(int seed, double sigma) override {
			VectorType pars(npar_);

			netket::RandomGaussian(pars, seed, sigma);
			SetParametersIdentity(pars);
		};

		int Nvisible() const override { return N_; };

		// Check lookups later
		void InitLookup(const Eigen::VectorXd &v, LookupType &lt) override {
			int site;

			// First (left) site
			_InitLookup_check(lt, 0);
			lt.M(0) = W_[confindex_[v(0)]];

			// Last (right) site
			_InitLookup_check(lt, 1);
			lt.M(1) = W_[d_ * ((N_ - 1) % symperiod_) + confindex_[v(N_ - 1)]];

			// Rest sites
			for (int i = 2; i < 2 * N_; i += 2) {
				_InitLookup_check(lt, i);
				site = i / 2;
				lt.M(i) = lt.M(i - 2) * W_[d_ * (site % symperiod_) + confindex_[v(site)]];

				_InitLookup_check(lt, i + 1);
				site = N_ - 1 - site;
				lt.M(i + 1) = W_[d_ * (site % symperiod_) + confindex_[v(site)]] * lt.M(i - 1);
			}
		};

		// Auxiliary function
		inline void _InitLookup_check(LookupType &lt, int i) {
			if (lt.MatrixSize() == i) {
				lt.AddMatrix(D_, D_);
			}
			else {
				lt.M(i).resize(D_, D_);
			}
		};

		// Auxiliary function for sorting indeces 
		// (copied from stackexchange - original answer by Lukasz Wiklendt)
		inline std::vector<std::size_t> sort_indeces(const std::vector<int> &v) {
			// initialize original index locations
			std::vector<std::size_t> idx(v.size());
			std::iota(idx.begin(), idx.end(), 0);
			// sort indexes based on comparing values in v
			std::sort(idx.begin(), idx.end(), [&v](std::size_t i1, std::size_t i2) {return v[i1] < v[i2]; });
			return idx;
		};

		// Check lookups later
		void UpdateLookup(const Eigen::VectorXd &v,
			const std::vector<int> &tochange,
			const std::vector<double> &newconf,
			LookupType &lt) override {
			
			std::size_t nchange = tochange.size();
			if (nchange <= 0) {
				return;
			}
			std::vector<std::size_t> sorted_ind = sort_indeces(tochange);
			int site = tochange[sorted_ind[0]];

			//InfoMessage() << "Lookup update called" << std::endl;
			//for (std::size_t k = 0; k < nchange; k++) {
			//	InfoMessage() << tochange[sorted_ind[k]] << std::endl;
			//}

			// Update left (site++)
			if (site == 0) {
				lt.M(0) = W_[confindex_[newconf[sorted_ind[0]]]];
			}
			else {
				lt.M(2 * site) = lt.M(2 * (site - 1)) * W_[d_ * (site % symperiod_) + confindex_[newconf[sorted_ind[0]]]];
			}
			for (std::size_t k = 1; k < nchange; k++) {
				for (site = tochange[sorted_ind[k - 1]] + 1; site < tochange[sorted_ind[k]]; site++) {
					lt.M(2 * site) = lt.M(2 * (site - 1)) * W_[d_ * (site % symperiod_) + confindex_[v(site)]];
				}
				site = tochange[sorted_ind[k]];
				lt.M(2 * site) = lt.M(2 * (site - 1)) * W_[d_ * (site % symperiod_) + confindex_[newconf[sorted_ind[k]]]];
			}
			for (int site = tochange[sorted_ind[nchange - 1]] + 1; site < N_; site++) {
				lt.M(2 * site) = lt.M(2 * (site - 1)) * W_[d_ * (site % symperiod_) + confindex_[v(site)]];
			}

			//InfoMessage() << "Lookup update left completed" << std::endl;

			// Update right (site--)
			site = tochange[sorted_ind[nchange - 1]];
			if (site == N_ - 1) {
				lt.M(1) = W_[d_ * ((N_ - 1) % symperiod_) + confindex_[newconf[sorted_ind[nchange - 1]]]];
			}
			else {
				lt.M(2 * (N_ - site) - 1) = W_[d_ * (site % symperiod_) + confindex_[newconf[sorted_ind[nchange - 1]]]] * lt.M(2 * (N_ - site) - 3);
			}

			//InfoMessage() << "First right assigned" << std::endl;

			for (int k = nchange - 2; k >= 0; k--) {
				for (site = tochange[sorted_ind[k + 1]] - 1; site > tochange[sorted_ind[k]]; site--) {
					lt.M(2 * (N_ - site) - 1) = W_[d_ * (site % symperiod_) + confindex_[v(site)]] * lt.M(2 * (N_ - site) - 3);
				}
				site = tochange[sorted_ind[k]];
				lt.M(2 * (N_ - site) - 1) = W_[d_ * (site % symperiod_) + confindex_[newconf[sorted_ind[k]]]] * lt.M(2 * (N_ - site) - 3);
			}

			//InfoMessage() << "Middle loops done" << std::endl;

			for (site = tochange[sorted_ind[0]] - 1; site >= 0; site--) {
				lt.M(2 * (N_ - site) - 1) = W_[d_ * (site % symperiod_) + confindex_[v(site)]] * lt.M(2 * (N_ - site) - 3);
			}

			//InfoMessage() << "Lookup update ended" << std::endl;
		};

		T LogVal(const Eigen::VectorXd &v) override {
			//InfoMessage() << "LogVal called" << std::endl;

			MatrixType p = W_[confindex_[v(0)]];
			for (int site = 1; site < N_; site++) {
				p *= W_[d_ * (site % symperiod_) + confindex_[v(site)]];
			}
			return std::log(p.trace());
		};

		// Ignore lookups for now
		T LogVal(const Eigen::VectorXd &v, const LookupType &lt) override {
			return LogVal(v);
		};

		VectorType LogValDiff(
			const Eigen::VectorXd &v, const std::vector<std::vector<int>> &tochange,
			const std::vector<std::vector<double>> &newconf) override {
			const std::size_t nconn = tochange.size();
			int site = 0;

			//InfoMessage() << "LogValDiff full called" << std::endl;

			std::size_t nchange;
			std::vector<std::size_t> sorted_ind;
			VectorType logvaldiffs=VectorType::Zero(nconn);
			StateType current_psi = mps_contraction(v, 0, N_).trace();
			MatrixType new_prods(D_, D_);

			// current_prod calculation only needs to be done once. Fix that
			for (std::size_t k = 0; k < nconn; k++) {
				nchange = tochange[k].size();

				//InfoMessage() << "k = " << k << " nchange = " << nchange << std::endl;

				if (nchange > 0) {
					sorted_ind = sort_indeces(tochange[k]);
					site = tochange[k][sorted_ind[0]];

					if (site == 0) {
						new_prods = W_[confindex_[newconf[k][sorted_ind[0]]]];
					}
					else {
						new_prods = mps_contraction(v, 0, site) * W_[d_ * (site % symperiod_) + confindex_[newconf[k][sorted_ind[0]]]];
					}

					for (std::size_t i = 1; i < nchange; i++) {
						site = tochange[k][sorted_ind[i]];
						new_prods *= mps_contraction(v, tochange[k][sorted_ind[i - 1]] + 1, site) * W_[d_ * (site % symperiod_) + confindex_[newconf[k][sorted_ind[i]]]];
					}
					site = tochange[k][sorted_ind[nchange - 1]];
					if (site < N_ - 1) {
						new_prods *= mps_contraction(v, site + 1, N_);
					}
					
					logvaldiffs(k) = std::log(new_prods.trace() / current_psi);
				}
			}
			
			//InfoMessage() << "LogValDiff full ended" << std::endl;
			
			return logvaldiffs;
		};

		//Auxiliary function that calculates contractions from site1 to site2
		inline MatrixType mps_contraction(const Eigen::VectorXd &v,
			const int &site1, const int &site2) {
			MatrixType c = MatrixType::Identity(D_, D_);
			for (int site = site1; site < site2; site++) {
				c *= W_[d_ * (site % symperiod_) + confindex_[v(site)]];
			}
			return c;
		};

		/** // Ignore lookups for now (copy the previous function)
		T LogValDiff(const Eigen::VectorXd &v, const std::vector<int> &toflip,
			const std::vector<double> &newconf,
			const LookupType &lt) override {
			const std::size_t nflip = toflip.size();
			if (nflip <= 0) {
				return T(0, 0);
			}

			//InfoMessage() << "LogValDiff lookup called" << std::endl;

			std::vector<std::size_t> sorted_ind = sort_indeces(toflip);
			StateType current_psi = mps_contraction(v, 0, N_).trace();
			MatrixType new_prods(D_, D_);

			if (toflip[0] == 0) {
				new_prods = W_[confindex_[newconf[sorted_ind[0]]]];
			}
			else {
				new_prods = mps_contraction(v, 0, toflip[sorted_ind[0]]) * W_[d_ * (toflip[sorted_ind[0]] % symperiod_) + confindex_[newconf[sorted_ind[0]]]];
			}
			for (std::size_t i = 1; i < nflip; i++) {
				new_prods *= mps_contraction(v, toflip[sorted_ind[i-1]] + 1, toflip[sorted_ind[i]]) * W_[d_ * (toflip[sorted_ind[i]] % symperiod_) + confindex_[newconf[sorted_ind[i]]]];
			}
			if (toflip[nflip - 1] < N_ - 1) {
				new_prods *= mps_contraction(v, toflip[sorted_ind[nflip - 1]] + 1, N_);
			}

			//InfoMessage() << "LogValDiff lookup called" << std::endl;

			return std::log(new_prods.trace() / current_psi);
		}; */

		T LogValDiff(const Eigen::VectorXd &v, const std::vector<int> &toflip,
			const std::vector<double> &newconf,
			const LookupType &lt) override {
			// Assumes that vector toflip is in ascending order
			const std::size_t nflip = toflip.size();
			if (nflip <= 0) {
				return T(0, 0);
			}
			MatrixType new_prod;
			std::vector<std::size_t> sorted_ind = sort_indeces(toflip);
			int site = toflip[sorted_ind[0]];
			
			//InfoMessage() << "LogValDiff lookup called" << std::endl;
			//for (std::size_t k = 0; k < nflip; k++) {
			//	InfoMessage() << toflip[k] << std::endl;
			//}

			if (site == 0) {
				new_prod = W_[confindex_[newconf[sorted_ind[0]]]];
			}
			else {
				new_prod = lt.M(2 * (site - 1)) * W_[d_ * (site % symperiod_) + confindex_[newconf[sorted_ind[0]]]];
			}

			for (std::size_t k = 1; k < nflip; k++) {
				site = toflip[sorted_ind[k]];
				new_prod *= mps_contraction(v, toflip[sorted_ind[k - 1]] + 1, site) * W_[d_ * (site % symperiod_) + confindex_[newconf[sorted_ind[k]]]];
			}

			//InfoMessage() << "LogValDiff lookup ended" << std::endl;

			site = toflip[sorted_ind[nflip - 1]];
			if (site == N_ - 1) {
				return std::log(new_prod.trace() / lt.M(2 * N_ - 2).trace()); // A log is needed here
			}
			else {
				return std::log((new_prod * lt.M(2 * (N_ - site) - 3)).trace()/ lt.M(2 * N_ - 2).trace()); // A log is needed here
			}
		};

		// Derivative with full calculation
		VectorType DerLog(const Eigen::VectorXd &v) override {
			const int Dsq = D_ * D_;
			//ComputeVtilde(v, vtilde_);
			MatrixType temp_product(D_, D_);
			std::vector<MatrixType> left_prods, right_prods;
			VectorType der = VectorType::Zero(npar_);

			//InfoMessage() << "Derivative called" << std::endl;
			
			// Calculate products
			left_prods.push_back(W_[confindex_[v(0)]]);
			right_prods.push_back(W_[d_ * ((N_ - 1) % symperiod_) + confindex_[v(N_ - 1)]]);
			for (int site = 1; site < N_- 1; site++) {
				left_prods.push_back(left_prods[site - 1] * W_[d_ * (site % symperiod_) + confindex_[v(site)]]);
				right_prods.push_back(W_[d_ * ((N_ - 1 - site) % symperiod_) + confindex_[v(N_ - 1 - site)]] * right_prods[site-1]);
			}
			left_prods.push_back(left_prods[N_ -2] * W_[d_ * ((N_-1) % symperiod_) + confindex_[v(N_ -1)]]);
			right_prods.push_back(W_[confindex_[v(0)]] * right_prods[N_ - 2]);

			der.segment(confindex_[v(0)] * Dsq, Dsq) += Eigen::Map<VectorType>((right_prods[N_ - 2]).transpose().data(), Dsq);
			for (int site = 1; site < N_ - 1; site++) {
				temp_product = right_prods[N_ - site - 2] * left_prods[site - 1];
				der.segment((d_ * (site % symperiod_) + confindex_[v(site)])* Dsq, Dsq) = Eigen::Map<VectorType>((temp_product).transpose().data(), Dsq);
			}
			der.segment((d_ * ((N_ - 1) % symperiod_) + confindex_[v(N_ - 1)])* Dsq, Dsq) = Eigen::Map<VectorType>((left_prods[N_ - 2]).transpose().data(), Dsq);

			//InfoMessage() << "Derivative ended" << std::endl;

			return der / left_prods[N_ - 1].trace();
		};

		// Json functions
		const Hilbert &GetHilbert() const { return hilbert_; };

		void to_json(json &j) const override {
		  j["Machine"]["Name"] = "MPSperiodic";
		  j["Machine"]["Nspins"] = N_;
		  j["Machine"]["BondDim"] = D_;
		  j["Machine"]["PhysDim"] = d_;
		  j["Machine"]["SymmetryPeriod"] = symperiod_;
		  j["Machine"]["W"] = W_;
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
			  //throw InvalidInputError("Unspecified period of symmetry");
			  // Default is symperiod = N, resp. no translational symmetry
			  symperiod_ = N_;
		  }

		  Init();

		  // Loading parameters, if defined in the input
		 // if (FieldExists(pars["Machine"], "W")) {
		//	W_ = pars["Machine"]["W"];
		  //}
		};
	};

} // namespace netket

#endif
