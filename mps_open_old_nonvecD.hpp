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

#ifndef NETKET_MPS_OPEN_HPP
#define NETKET_MPS_OPEN_HPP

namespace netket {

	template <typename T>
	class MPSOpen : public AbstractMachine<T> {
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

		// MPS Matrices (stored as [N * d, D, D]
		std::vector<MatrixType> W_;

		// Map from Hilbert states to MPS indices
		std::map<double, int> confindex_;
		
		// Local lookup matrices
		//std::vector<MatrixType> loc_lt;
		// Local vectors to transform {-1, 1} to {0, 1}
		//Eigen::VectorXd vtilde_;

		const Hilbert &hilbert_;

	public:
		using StateType = T;
		using LookupType = Lookup<T>;

		// constructor
		explicit MPSOpen(const Hilbert &hilbert, const json &pars)
			: N_(hilbert.Size()), hilbert_(hilbert), d_(hilbert.LocalSize()) {
			from_json(pars);
		}

		// Auxiliary function that defines the matrices
		void Init() {
			MatrixType left = MatrixType::Zero(1, D_), middle = MatrixType::Zero(D_, D_), right = MatrixType::Zero(D_, 1);
			npar_ = (N_ - 2) * d_ * D_ * D_ + 2 * d_ * D_;

			W_.push_back(left);
			for (int i = 0; i < N_*d_; i++) {
				W_.push_back(middle);
			}
			W_.push_back(right);

			// Machine creation messages
			InfoMessage() << "Open MPS machine with " << N_ << " sites created" << std::endl;
			InfoMessage() << "Physical dimension d = " << d_ << " and bond dimension D = " << D_ << std::endl;
			InfoMessage() << "The number of variational parameters is " << npar_ << std::endl;

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

		int Npar() const override { return npar_; };

		VectorType GetParameters() override {
			int k = 0;
			VectorType pars(npar_);

			// Left boundary
			for (int p = 0; p < d_; p++) {
				for (int j = 0; j < D_; j++) {
					pars(k) = W_[p](0, j);
					k++;
				}
			}

			// Middle
			for (int p = d_; p < (N_ - 1)*d_; p++) {
				for (int i = 0; i < D_; i++) {
					for (int j = 0; j < D_; j++) {
						pars(k) = W_[p](i, j);
						k++;
					}
				}
			}

			// Right boundary
			for (int p = (N_ - 1)*d_; p < N_ * d_; p++) {
				for (int i = 0; i < D_; i++) {
					pars(k) = W_[p](i, 0);
					k++;
				}
			}

			return pars;
		};

		void SetParameters(const VectorType &pars) override {
			int k = 0;

			// Left boundary
			for (int p = 0; p < d_; p++) {
				for (int j = 0; j < D_; j++) {
					W_[p](0, j) = pars(k);
					k++;
				}
			}

			// Middle
			for (int p = d_; p < (N_ - 1)*d_; p++) {
				for (int i = 0; i < D_; i++) {
					for (int j = 0; j < D_; j++) {
						W_[p](i, j) = pars(k);
						k++;
					}
				}
			}

			// Right boundary
			for (int p = (N_ - 1)*d_; p < N_ * d_; p++) {
				for (int i = 0; i < D_; i++) {
					W_[p](i, 0) = pars(k);
					k++;
				}
			}
		};

		// Auxiliary function used for setting initial random parameters and adding identities in every matrix
		inline void SetParametersIdentity(const VectorType &pars) {
			int k = 0;

			// Left boundary
			for (int p = 0; p < d_; p++) {
				for (int j = 0; j < D_; j++) {
					W_[p](0, j) = pars(k) + T(1, 0);
					k++;
				}
			}

			// Middle
			for (int p = d_; p < (N_ - 1)*d_; p++) {
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

			// Right boundary
			for (int p = (N_ - 1)*d_; p < N_ * d_; p++) {
				for (int i = 0; i < D_; i++) {
					W_[p](i, 0) = pars(k) + T(1, 0);
					k++;
				}
			}
		};

		void InitRandomPars(int seed, double sigma) override {
			VectorType pars(npar_);

			netket::RandomGaussian(pars, seed, sigma);
			SetParametersIdentity(pars);
		};

		int Nvisible() const override { return N_; };

		void InitLookup(const Eigen::VectorXd &v, LookupType &lt) override {
			int site;
			//ComputeVtilde(v, vtilde_);
			// Initializes local lookup too! (commented for now)
			//std::vector<MatrixType> loc_lt;

			// First (left) site
			_InitLookupLeft_check(lt, 0);
			lt.M(0) = W_[confindex_[v(0)]];
			//loc_lt.push_back(W_[(int)vtilde_(0)]);

			// Last (right) site
			_InitLookupRight_check(lt, 1);
			lt.M(1) = W_[d_ * (N_ - 1) + confindex_[v(N_ - 1)]];
			//loc_lt.push_back(W_[d_ * (N_ - 1) + (int)vtilde_(N_ - 1)]);

			// Rest sites
			for (int i = 2; i < 2 * (N_ - 1); i += 2) {
				_InitLookupLeft_check(lt, i);
				site = i / 2;
				lt.M(i) = lt.M(i - 2) * W_[d_ * site + confindex_[v(site)]];
				//loc_lt.push_back(lt.M(i));

				_InitLookupRight_check(lt, i + 1);
				site = N_ - 1 - site;
				lt.M(i + 1) = W_[d_ * site + confindex_[v(site)]] * lt.M(i - 1);
				//loc_lt.push_back(lt.M(i + 1));
			}

			// Last lookups which are just numbers
			_InitLookupBoundary_check(lt, 2 * N_ - 2);
			lt.M(2 * N_ - 2) = lt.M(2 * N_ - 4) * W_[d_ * (N_ - 1) + confindex_[v(N_ - 1)]];
			_InitLookupBoundary_check(lt, 2 * N_ - 1);
			lt.M(2 * N_ - 1) = W_[confindex_[v(0)]] * lt.M(2 * N_ - 3);

			//InfoMessage() << "InitLookup update ended" << std::endl;
		};

		// Auxiliary function
		inline void _InitLookupLeft_check(LookupType &lt, int i) {
			if (lt.MatrixSize() == i) {
				lt.AddMatrix(1, D_);
			}
			else {
				lt.M(i).resize(1, D_);
			}
		};

		// Auxiliary function
		inline void _InitLookupRight_check(LookupType &lt, int i) {
			if (lt.MatrixSize() == i) {
				lt.AddMatrix(D_, 1);
			}
			else {
				lt.M(i).resize(D_, 1);
			}
		};

		// Auxiliary function (the last lookups are just numbers - shape=(1,1))
		inline void _InitLookupBoundary_check(LookupType &lt, int i) {
			if (lt.MatrixSize() == i) {
				lt.AddMatrix(1, 1);
			}
			else {
				lt.M(i).resize(1, 1);
			}
		};

		//Auxiliary function that calculates all left and right contractions (not ready yet)
		/**
		inline std::vector<MatrixType> initialize_left_right(const Eigen::VectorXd &vtilde) {
			int site;
			std::vector<MatrixType> prods;

			prods.push_back(W_[(int)vtilde(0)]);
			prods.push_back(W_[(int)vtilde(N_ - 1)]);

			for (int i = 2; i < 2 * N_; i += 2) {
				site = i / 2;
				prods.push_back(prods[i - 2] * W_[d_ * site + (int)vtilde(site)]);
				prods.push_back(W_[(int)vtilde(site) * prods[i - 2]]);
			}
			return prods;
		}; */

		// Auxiliary function for sorting indeces 
		// (copied from stackexchange - original answer by Lukasz Wiklendt)
		inline std::vector<std::size_t> sort_indeces(const std::vector<int> &v) {
			// initialize original index locations
			std::vector<std::size_t> idx(v.size());
			std::iota(idx.begin(), idx.end(), 0);
			// sort indexes based on comparing values in v
			std::sort(idx.begin(), idx.end(), [&v](std::size_t i1, std::size_t i2) {return v[i1] < v[i2];});
			return idx;
		};

		void UpdateLookup(const Eigen::VectorXd &v,
			const std::vector<int> &tochange,
			const std::vector<double> &newconf,
			LookupType &lt) override {
			//ComputeVtilde(v, vtilde_);
			// Updates local lookup too! (commented for now)
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
				//loc_lt[0] = lt.M(0);
			}
			else {
				lt.M(2 * site) = lt.M(2 * (site - 1)) * W_[d_ * site + confindex_[newconf[sorted_ind[0]]]];
				//loc_lt[2 * site] = lt.M(2 * site);
			}
			for (std::size_t k = 1; k < nchange; k++) {
				for (site = tochange[sorted_ind[k - 1]] + 1; site < tochange[sorted_ind[k]]; site++) {
					lt.M(2 * site) = lt.M(2 * (site - 1)) * W_[d_ * site + confindex_[v(site)]];
					//loc_lt[2 * site] = lt.M(2 * site);
				}
				site = tochange[sorted_ind[k]];
				lt.M(2 * site) = lt.M(2 * (site - 1)) * W_[d_ * site + confindex_[newconf[sorted_ind[k]]]];
				//loc_lt[2 * site] = lt.M(2 * site);
			}
			for (int site = tochange[sorted_ind[nchange - 1]] + 1; site < N_; site++) {
				lt.M(2 * site) = lt.M(2 * (site - 1)) * W_[d_ * site + confindex_[v(site)]];
				//loc_lt[2 * site] = lt.M(2 * site);
			}

			
			//InfoMessage() << "Lookup update left completed" << std::endl;

			// Update right (site--)
			site = tochange[sorted_ind[nchange - 1]];
			if (site == N_ - 1) {
				lt.M(1) = W_[d_ * (N_ - 1) + confindex_[newconf[sorted_ind[nchange - 1]]]];
			}
			else {
				lt.M(2 * (N_ - site) - 1) = W_[d_ * site + confindex_[newconf[sorted_ind[nchange - 1]]]] * lt.M(2 * (N_ - site) - 3);
			}

			//InfoMessage() << "First right assigned" << std::endl;

			for (int k = nchange - 2; k >= 0; k--) {
				for (site = tochange[sorted_ind[k + 1]] - 1; site > tochange[sorted_ind[k]]; site--) {
					lt.M(2 * (N_ - site) - 1) = W_[d_ * site + confindex_[v(site)]] * lt.M(2 * (N_ - site) - 3);
				}
				site = tochange[sorted_ind[k]];
				lt.M(2 * (N_ - site) - 1) = W_[d_ * site + confindex_[newconf[sorted_ind[k]]]] * lt.M(2 * (N_ - site) - 3);
			}

			//InfoMessage() << "Middle loops done" << std::endl;

			for (site = tochange[sorted_ind[0]] - 1; site >= 0; site--) {
				lt.M(2 * (N_ - site) - 1) = W_[d_ * site + confindex_[v(site)]] * lt.M(2 * (N_ - site) - 3);
			}

			//InfoMessage() << "Lookup update ended" << std::endl;
		};

		T LogVal(const Eigen::VectorXd &v) override {
			//ComputeVtilde(v, vtilde_);
			MatrixType p = W_[confindex_[v(0)]];
			for (int site = 1; site < N_ - 1; site++) { 
				p *= W_[d_ * site + confindex_[v(site)]];
			}
			return std::log((p * W_[d_ * (N_ - 1) + confindex_[v(N_ - 1)]]).trace());
		};

		T LogVal(const Eigen::VectorXd &v, const LookupType &lt) override {
			return std::log(lt.M(2 * N_ - 2).trace());
		};

		VectorType LogValDiff(
			const Eigen::VectorXd &v, const std::vector<std::vector<int>> &tochange,
			const std::vector<std::vector<double>> &newconf) override {
			const std::size_t nconn = tochange.size();
			int site = 0;
			//ComputeVtilde(v, vtilde_);
			std::size_t nchange;
			std::vector<std::size_t> sorted_ind;
			VectorType logvaldiffs=VectorType::Zero(nconn);
			StateType current_psi = mps_contractionLeft(v, N_).trace();
			MatrixType new_prods(1, D_);

			//InfoMessage() << "LogValDiff full called" << std::endl;

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
						new_prods = mps_contractionLeft(v, site) * W_[d_ * site + confindex_[newconf[k][sorted_ind[0]]]];
					}
					for (std::size_t i = 1; i < nchange; i++) {
						site = tochange[k][sorted_ind[i]];
						new_prods *= mps_contraction(v, tochange[k][sorted_ind[i - 1]] + 1, site) * W_[d_ * site + confindex_[newconf[k][sorted_ind[i]]]];
					}
					site = tochange[k][sorted_ind[nchange - 1]];
					if (site < N_ - 1) {
						logvaldiffs(k) = std::log((new_prods * mps_contractionRight(v, site + 1)).trace() / current_psi);
					}
					else {
						logvaldiffs(k) = std::log(new_prods.trace() / current_psi);
					}
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
				c *= W_[d_ * site + confindex_[v(site)]];
			}
			return c;
		};

		inline MatrixType mps_contractionLeft(const Eigen::VectorXd &v, const int &site) {
			MatrixType c = W_[confindex_[v(0)]];
			for (int i = 1; i < site; i++) {
				c *= W_[d_ * i + confindex_[v(i)]];
			}
			return c;
		};

		inline MatrixType mps_contractionRight(const Eigen::VectorXd &v, const int &site) {
			MatrixType c = W_[(N_ - 1) * d_ + confindex_[v(N_ - 1)]];
			for (int i = N_ - 2; i >= site; i--) {
				c = W_[d_ * i + confindex_[v(i)]] * c;
			}
			return c;
		};


		// Ignore lookups for now (copy the previous function)
		/**
		T LogValDiff(const Eigen::VectorXd &v, const std::vector<int> &toflip,
			const std::vector<double> &newconf,
			const LookupType &lt) override {
			const std::size_t nflip = toflip.size();
			ComputeVtilde(v, vtilde_);
			if (nflip <= 0) {
				return T(0, 0);
			}
			StateType current_psi = mps_contraction(vtilde_, 0, N_).trace();
			MatrixType new_prods(D_, D_);

			if (toflip[0] == 0) {
				new_prods = W_[confindex_[newconf[0]]];
			}
			else {
				new_prods = mps_contraction(vtilde_, 0, toflip[0]) * W_[d_ * toflip[0] + confindex_[newconf[0]]];
			}
			for (std::size_t i = 1; i < nflip; i++) {
				//InfoMessage() << "toflip = " << toflip[i] << std::endl;
				new_prods *= mps_contraction(vtilde_, toflip[i - 1] + 1, toflip[i]) * W_[d_ * toflip[i] + confindex_[newconf[i]]];
			}
			if (toflip[nflip - 1] < N_ - 1) {
				new_prods *= mps_contraction(vtilde_, toflip[nflip - 1] + 1, N_);
			}
			//InfoMessage() << "LogValDiff lookup ended " << std::log(new_prods.trace() / current_psi) << std::endl;
			return std::log(new_prods.trace() / current_psi);
		};*/

		
		T LogValDiff(const Eigen::VectorXd &v, const std::vector<int> &toflip,
			const std::vector<double> &newconf,
			const LookupType &lt) override {
			// Assumes that vector toflip is in ascending order
			//ComputeVtilde(v, vtilde_);
			const std::size_t nflip = toflip.size();
			if (nflip <= 0) {
				return T(0, 0);
			}
			const std::vector<std::size_t> sorted_ind = sort_indeces(toflip);
			int site = toflip[sorted_ind[0]];
			MatrixType new_prod;
			
			//InfoMessage() << "LogValDiff lookup called" << std::endl;
			//for (std::size_t k = 0; k < nflip; k++) {
			//	InfoMessage() << toflip[k] << std::endl;
			//}

			if (site == 0) {
				new_prod = W_[confindex_[newconf[sorted_ind[0]]]];
			}
			else {
				new_prod = lt.M(2 * (site - 1)) * W_[d_ * site + confindex_[newconf[sorted_ind[0]]]];
			}

			for (std::size_t k = 1; k < nflip; k++) {
				site = toflip[sorted_ind[k]];
				new_prod *= mps_contraction(v, toflip[sorted_ind[k - 1]] + 1, site) * W_[d_ * site + confindex_[newconf[sorted_ind[k]]]];
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

		// Derivative that uses local lookups (deleted)

		// Derivative with full calculation
		VectorType DerLog(const Eigen::VectorXd &v) override {
			const int Dsq = D_ * D_;
			//ComputeVtilde(v, vtilde_);
			std::vector<MatrixType> left_prods, right_prods;
			VectorType der = VectorType::Zero(npar_);

			//InfoMessage() << "Derivative called" << std::endl;

			// Calculate products
			left_prods.push_back(W_[confindex_[v(0)]]);
			right_prods.push_back(W_[d_ * (N_ - 1) + confindex_[v(N_ - 1)]]);
			for (int site = 1; site < N_- 1; site++) {

				//InfoMessage() << "Right shape = " << right_prods[site-1] << std::endl;
				//InfoMessage() << "Left shape = " << left_prods[site - 1] << std::endl;

				left_prods.push_back(left_prods[site - 1] * W_[d_ * site + confindex_[v(site)]]);
				right_prods.push_back(W_[d_ * (N_ - 1 - site) + confindex_[v(N_ - 1 - site)]] * right_prods[site-1]);
			}
			left_prods.push_back(left_prods[N_ -2] * W_[d_ * (N_-1) + confindex_[v(N_-1)]]);
			right_prods.push_back(W_[confindex_[v(0)]] * right_prods[N_ - 2]);

			//InfoMessage() << "Products calculated" << std::endl;

			der.segment(confindex_[v(0)] * D_, D_) = Eigen::Map<VectorType>((right_prods[N_ - 2]).data(), D_);
			//InfoMessage() << "Left derivative assigned" << std::endl;
			for (int site = 1; site < N_ - 1; site++) {
				der.segment(d_ * D_ + ((site-1) * d_ + confindex_[v(site)])* Dsq, Dsq) = middle_tensor_product(left_prods[site - 1], right_prods[N_ - site - 2]);

				//InfoMessage() << "site = " << site << std::endl;
			}
			der.segment(d_ * D_ + (N_ - 2) * d_ * Dsq + confindex_[v(N_ - 1)]* D_, D_) = Eigen::Map<VectorType>((left_prods[N_ - 2]).data(), D_);

			//InfoMessage() << "Derivative ended" << std::endl;
			//der = der / left_prods[N_ - 1].trace();

			return der / left_prods[N_ - 1].trace();
		};

		inline VectorType middle_tensor_product(const MatrixType left, const MatrixType right) {
			VectorType der_seg = VectorType::Zero(D_ * D_);
			int k = 0;
			for (int i = 0; i < D_; i++) {
				for (int j = 0; j < D_; j++) {
					der_seg(k) = left(0, i) * right(j, 0);
					k++;
				}
			}
			return der_seg;
		};

		// Json functions
		const Hilbert &GetHilbert() const { return hilbert_; };

		void to_json(json &j) const override {
		  j["Machine"]["Name"] = "MPSopen";
		  j["Machine"]["Nspins"] = N_;
		  j["Machine"]["BondDim"] = D_;
		  j["Machine"]["PhysDim"] = d_;
		  j["Machine"]["W"] = W_;
		}; 

		void from_json(const json &pars) override {
		  if (pars.at("Machine").at("Name") != "MPSopen") {
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

		  Init();

		  // Loading parameters, if defined in the input
		 // if (FieldExists(pars["Machine"], "W")) {
		//	W_ = pars["Machine"]["W"];
		  //}
		};

		// Still to do: 
		// Do vectorization for spins more than 1/2. (Done in a sketchy way)
		// Don't forget the logarithms where needed. (Done)
		// Upgrade look ups to be more efficient (tree calculation?).
		// Look ups for LogDer?.
		// Slight modifications for OBC.
	};

} // namespace netket

#endif
