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

#ifndef NETKET_MPS_CALCULATOR_HPP
#define NETKET_MPS_CALCULATOR_HPP
// A class that assists MPS and SBS calculations

namespace netket {

	template <typename T>
	class MPSCalculator {
		using VectorType = typename AbstractMachine<T>::VectorType;
		using MatrixType = typename AbstractMachine<T>::MatrixType;

		// Number of sites in the MPS
		int N_;
		// Physical dimension
		int d_;
		// Bond dimension (currently keep constant along the MPS)
		int D_;
		// Number of variational parameters in this MPS
		int npar_;
		// Period of translational symmetry (has to be a divisor of N)
		int symperiod_;

		// Map from Hilbert states to MPS indices
		//std::map<double, int> confindex_;

		// MPS Matrices (stored as [symperiod, d, D, D]
		std::vector<std::vector<MatrixType>> W_;

	public:
		using StateType = T;
		using LookupType = Lookup<T>;

		// constructor
		explicit MPSCalculator(const int &N, const int &d, const int &D, const int &symperiod)
			: N_(N), d_(d), D_(D), symperiod_(symperiod) {
			Init();
		};

		void Init() {
			// Initialize parameters
			std::vector<MatrixType> pushback_vec;
			MatrixType init_mat = MatrixType::Zero(D_, D_);
			npar_ = symperiod_ * d_ * D_ * D_;

			for (int site = 0; site < symperiod_; site++) {
				W_.push_back(pushback_vec);
				for (int spin = 0; spin < d_; spin++) {
					W_[site].push_back(init_mat);
				}
			}
		};

		int Npar() {
			return npar_;
		};

		VectorType GetParameters() {
			int k = 0;
			VectorType pars(npar_);

			for (int site = 0; site < symperiod_; site++) {
				for (int spin = 0; spin < d_; spin++) {
					for (int i = 0; i < D_; i++) {
						for (int j = 0; j < D_; j++) {
							pars(k) = W_[site][spin](i, j);
							k++;
						}
					}
				}
			}
			return pars;
		};

		void SetParameters(const VectorType &pars) {
			int k = 0;

			for (int site = 0; site < symperiod_; site++) {
				for (int spin = 0; spin < d_; spin++) {
					for (int i = 0; i < D_; i++) {
						for (int j = 0; j < D_; j++) {
							W_[site][spin](i, j) = pars(k);
							k++;
						}
					}
				}
			}
		};

		// Auxiliary function used for setting initial random parameters and adding identities in every matrix
		inline void SetParametersIdentity(const VectorType &pars) {
			int k = 0;

			for (int site = 0; site < symperiod_; site++) {
				for (int spin = 0; spin < d_; spin++) {
					for (int i = 0; i < D_; i++) {
						for (int j = 0; j < D_; j++) {
							W_[site][spin](i, j) = pars(k);
							if (i == j) {
								W_[site][spin](i, j) += T(1, 0);
							}
							k++;
						}
					}
				}
			}
		};

		int Nsites() {
			return N_;
		};

		void InitLookup(const std::vector<int> &v, LookupType &lt, const int &start_ind) {
			int site;

			// We need 2 * L matrix lookups for each string, where L is the string's length
			// other than that lookups work exactly as in the MPS case

			// First (left) site
			_InitLookup_check(lt, start_ind);
			lt.M(start_ind) = W_[0][v[0]];

			// Last (right) site
			_InitLookup_check(lt, start_ind+1);
			lt.M(start_ind + 1) = W_[(N_ - 1) % symperiod_][v[N_ - 1]];

			// Rest sites
			for (int i = 2; i < 2 * N_; i += 2) {
				_InitLookup_check(lt, start_ind + i);
				site = i / 2;
				lt.M(start_ind + i) = lt.M(start_ind + i - 2) * W_[(site % symperiod_)][v[site]];

				_InitLookup_check(lt, start_ind + i + 1);
				site = N_ - 1 - site;
				lt.M(start_ind + i + 1) = W_[site % symperiod_][v[site]] * lt.M(start_ind + i - 1);
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

		void UpdateLookup(const std::vector<int> &v,
			const std::vector<int> &tochange,
			const std::vector<int> &newconf,
			LookupType &lt,
			const int &start_ind) {

			std::size_t nchange = tochange.size();
			if (nchange <= 0) {
				return;
			}
			std::vector<std::size_t> sorted_ind = sort_indeces(tochange);
			int site = tochange[sorted_ind[0]];

			//InfoMessage() << "Lookup update called" << std::endl;
			//for (std::size_t k = 0; k < nchange; k++) {
			//	InfoMessage() << tochange[sorted_ind[k]] << " , " << newconf[sorted_ind[k]] << std::endl;
			//}

			// Update left (site++)
			if (site == 0) {
				lt.M(start_ind) = W_[0][newconf[sorted_ind[0]]];
			}
			else {
				lt.M(start_ind + 2 * site) = lt.M(start_ind + 2 * (site - 1)) * W_[site % symperiod_][newconf[sorted_ind[0]]];
			}

			//InfoMessage() << "Lookup check1" << std::endl;

			for (std::size_t k = 1; k < nchange; k++) {
				for (site = tochange[sorted_ind[k - 1]] + 1; site < tochange[sorted_ind[k]]; site++) {
					lt.M(start_ind + 2 * site) = lt.M(start_ind + 2 * (site - 1)) * W_[site % symperiod_][v[site]];
				}
				site = tochange[sorted_ind[k]];
				lt.M(start_ind + 2 * site) = lt.M(start_ind + 2 * (site - 1)) * W_[site % symperiod_][newconf[sorted_ind[k]]];
			}

			//InfoMessage() << "Lookup check2" << std::endl;

			for (site = tochange[sorted_ind[nchange - 1]] + 1; site < N_; site++) {
				lt.M(start_ind + 2 * site) = lt.M(start_ind + 2 * (site - 1)) * W_[site % symperiod_][v[site]];
			}

			//InfoMessage() << "Lookup update left completed" << std::endl;

			// Update right (site--)
			site = tochange[sorted_ind[nchange - 1]];
			if (site == N_ - 1) {
				lt.M(start_ind + 1) = W_[(N_ - 1) % symperiod_][newconf[sorted_ind[nchange - 1]]];
			}
			else {
				lt.M(start_ind + 2 * (N_ - site) - 1) = W_[site % symperiod_][newconf[sorted_ind[nchange - 1]]] * lt.M(start_ind + 2 * (N_ - site) - 3);
			}

			//InfoMessage() << "First right assigned" << std::endl;

			for (int k = nchange - 2; k >= 0; k--) {
				for (site = tochange[sorted_ind[k + 1]] - 1; site > tochange[sorted_ind[k]]; site--) {
					lt.M(start_ind + 2 * (N_ - site) - 1) = W_[site % symperiod_][v[site]] * lt.M(start_ind + 2 * (N_ - site) - 3);
				}
				site = tochange[sorted_ind[k]];
				lt.M(start_ind + 2 * (N_ - site) - 1) = W_[site % symperiod_][newconf[sorted_ind[k]]] * lt.M(start_ind + 2 * (N_ - site) - 3);
			}

			//InfoMessage() << "Middle loops done" << std::endl;

			for (site = tochange[sorted_ind[0]] - 1; site >= 0; site--) {
				lt.M(start_ind + 2 * (N_ - site) - 1) = W_[site % symperiod_][v[site]] * lt.M(start_ind + 2 * (N_ - site) - 3);
			}

			//InfoMessage() << "Lookup update ended" << std::endl;
		};

		T FullProduct(const std::vector<int> &v) {
			MatrixType p = W_[0][v[0]];
			for (int site = 1; site < N_; site++) {
				p *= W_[site % symperiod_][v[site]];
			};
			return p.trace();
		};

		T LookupProduct(const LookupType &lt, const int &start_ind) {
			return lt.M(start_ind + 2 * N_ - 2).trace();
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

		//Auxiliary function that calculates contractions from site1 to site2
		inline MatrixType mps_contraction(const std::vector<int> &v,
			const int &site1, const int &site2) {
			MatrixType c = MatrixType::Identity(D_, D_);
			for (int site = site1; site < site2; site++) {
				c *= W_[site % symperiod_][v[site]];
			}
			return c;
		};

		// No k-dependent version
		T LogValDiff(const std::vector<int> &v, const std::vector<int> &toflip,
			const std::vector<int> &newconf) {

			const std::size_t nflip = toflip.size();
			if (nflip <= 0) {
				return T(0, 0);
			}

			std::vector<std::size_t> sorted_ind = sort_indeces(toflip);
			StateType current_psi = mps_contraction(v, 0, N_).trace();
			MatrixType new_prods(D_, D_);

			if (toflip[sorted_ind[0]] == 0) {
				new_prods = W_[0][newconf[sorted_ind[0]]];
			}
			else {
				new_prods = mps_contraction(v, 0, toflip[sorted_ind[0]]) * W_[toflip[sorted_ind[0]] % symperiod_][newconf[sorted_ind[0]]];
			}
			for (std::size_t i = 1; i < nflip; i++) {
				//InfoMessage() << "toflip = " << toflip[i] << std::endl;
				new_prods *= mps_contraction(v, toflip[sorted_ind[i - 1]] + 1, toflip[sorted_ind[i]]) * W_[toflip[sorted_ind[i]] % symperiod_][newconf[sorted_ind[i]]];
			}
			if (toflip[sorted_ind[nflip - 1]] < N_ - 1) {
				new_prods *= mps_contraction(v, toflip[sorted_ind[nflip - 1]] + 1, N_);
			}

			//InfoMessage() << "LogValDiff lookup ended " << std::log(new_prods.trace() / current_psi) << std::endl;

			return std::log(new_prods.trace() / current_psi);
		};

		T LogValDiff(const std::vector<int> &v, const std::vector<int> &toflip,
			const std::vector<int> &newconf,
			const LookupType &lt, const int &startind) {

			const std::size_t nflip = toflip.size();
			if (nflip <= 0) {
				return T(0, 0);
			}

			std::vector<std::size_t> sorted_ind = sort_indeces(toflip);
			MatrixType new_prods(D_, D_);
			int site = toflip[sorted_ind[0]];

			if (site == 0) {
				new_prods = W_[0][newconf[sorted_ind[0]]];
			}
			else {
				new_prods = lt.M(startind + 2 * (site - 1)) * W_[site % symperiod_][newconf[sorted_ind[0]]];
			}

			for (std::size_t i = 1; i < nflip; i++) {
				site = toflip[sorted_ind[i]];
				new_prods *= mps_contraction(v, toflip[sorted_ind[i - 1]] + 1, site) * W_[site % symperiod_][newconf[sorted_ind[i]]];
			}

			site = toflip[sorted_ind[nflip - 1]];
			if (site < N_ - 1) {
				//new_prods *= mps_contraction(v, toflip[sorted_ind[nflip - 1]] + 1, N_);
				new_prods *= lt.M(startind + 2 * (N_ - site) - 3);
			}

			//InfoMessage() << "LogValDiff lookup ended " << std::log(new_prods.trace() / current_psi) << std::endl;

			return std::log(new_prods.trace() / lt.M(startind + 2 * N_ - 2).trace());

		};

		// Lookup LogValDiff that doesn't require v
		// works only if nflip = 1
		T FastLogValDiff(const std::vector<int> &toflip,
			const std::vector<int> &newconf,
			const LookupType &lt, const int &startind) {

			T new_trace;

			//InfoMessage() << "FastLogValDiff lookup called" << std::endl;
			//for (std::size_t k = 0; k < nflip; k++) {
			//	InfoMessage() << toflip[k] << std::endl;
			//}

			if (toflip[0] == 0) {
				new_trace = (W_[0][newconf[0]] * lt.M(startind + 2 * N_ - 3)).trace();
			}
			else if (toflip[0] == N_ - 1) {
				new_trace = (lt.M(startind + 2 * (N_ - 2)) * W_[(N_ - 1) % symperiod_][newconf[0]]).trace();
			}
			else {
				new_trace = (lt.M(startind + 2 * (toflip[0] - 1)) * W_[toflip[0] % symperiod_][newconf[0]] * lt.M(startind + 2 * (N_ - toflip[0]) - 3)).trace();
			}

			return std::log(new_trace / lt.M(startind + 2 * N_ - 2).trace());
		};

		// Derivative with full calculation
		VectorType DerLog(const std::vector<int> &v) {
			const int Dsq = D_ * D_;
			//ComputeVtilde(v, vtilde_);
			MatrixType temp_product(D_, D_);
			std::vector<MatrixType> left_prods, right_prods;
			VectorType der = VectorType::Zero(npar_);

			//InfoMessage() << "Derivative called" << std::endl;
			// Calculate products
			left_prods.push_back(W_[0][v[0]]);
			right_prods.push_back(W_[(N_ - 1) % symperiod_][v[N_ - 1]]);
			for (int site = 1; site < N_ - 1; site++) {
				left_prods.push_back(left_prods[site - 1] * W_[site][v[site]]);
				right_prods.push_back(W_[(N_ - 1 - site) % symperiod_][v[N_ - 1 - site]] * right_prods[site - 1]);
			}
			left_prods.push_back(left_prods[N_ - 2] * W_[(N_ - 1) % symperiod_][v[N_ - 1]]);
			right_prods.push_back(W_[0][v[0]] * right_prods[N_ - 2]);

			der.segment(v[0] * Dsq, Dsq) += Eigen::Map<VectorType>((right_prods[N_ - 2]).transpose().data(), Dsq);
			for (int site = 1; site < N_ - 1; site++) {
				temp_product = right_prods[N_ - site - 2] * left_prods[site - 1];
				der.segment((d_ * (site % symperiod_) + v[site])*Dsq, Dsq) += Eigen::Map<VectorType>((temp_product).transpose().data(), Dsq);
			}
			der.segment((d_ * ((N_ - 1) % symperiod_) + v[N_ - 1])*Dsq, Dsq) += Eigen::Map<VectorType>((left_prods[N_ - 2]).transpose().data(), Dsq);

			//InfoMessage() << "Derivative ended" << std::endl;

			return der / left_prods[N_ - 1].trace();
		};
	};
} // namespace netket

#endif
