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
			MatrixType init_mat = MatrixType::Zeros(D, D);
			npar_ = symperiod_ * d_ * D_ * D_;

			for (int site = 0; site < symperiod_; site++) {
				W_.push_back(std::vector<MatrixType> x);
				for (int spin = 0; spin < d_; spin++) {
					W_[site].push_back(init_mat);
				}
			}
		};

		int Npar() {
			return npar_;
		};

		T FullProduct(const std::vector<int> &v) {
			MatrixType p = W_[0][v[0]];
			for (int site = 0; site < N_; site++) {
				p *= W[site % symperiod_][v[site]]
			};
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

		VectorType LogValDiff(
			const std::vector<int> &v, const std::vector<std::vector<int>> &tochange,
			const std::vector<std::vector<int>> &newconf) {

			const std::size_t nconn = tochange.size();
			VectorType logvaldiffs = VectorType::Zero(nconn);
			if (nconn <= 0) {
				return VectorType::Zero(nconn);
			}

			int site = 0;
			std::size_t nchange;
			std::vector<std::size_t> sorted_ind;
			StateType current_psi = mps_contraction(v, 0, N_).trace();
			MatrixType new_prods(D_, D_);

			for (std::size_t k = 0; k < nconn; k++) {
				nchange = tochange[k].size();

				//InfoMessage() << "k = " << k << " nchange = " << nchange << std::endl;
				if (nchange > 0) {
					sorted_ind = sort_indeces(tochange[k]);
					site = tochange[k][sorted_ind[0]];

					if (site == 0) {
						new_prods = W_[0][newconf[k][sorted_ind[0]]];
					}
					else {
						new_prods = mps_contraction(v, 0, site) * W_[site % symperiod_][newconf[k][sorted_ind[0]]];
					}

					for (std::size_t i = 1; i < nchange; i++) {
						site = tochange[k][sorted_ind[i]];
						new_prods *= mps_contraction(v, tochange[k][sorted_ind[i - 1]] + 1, site) * W_[site % symperiod_][newconf[k][sorted_ind[i]]];
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

		// Ignore lookups for now
		T LogValDiff(const std::vector<int> &v, const std::vector<int> &toflip,
			const std::vector<int> &newconf) {

			const std::size_t nflip = toflip.size();
			if (nflip <= 0) {
				return T(0, 0);
			}

			std::vector<std::size_t> sorted_ind = sort_indeces(toflip);
			StateType current_psi = mps_contraction(v, 0, N_).trace();
			MatrixType new_prods(D_, D_);

			if (toflip[0] == 0) {
				new_prods = W_[0][newconf[0]];
			}
			else {
				new_prods = mps_contraction(v, 0, toflip[sorted_ind[0]]) * W_[toflip[0] % symperiod_][newconf[sorted_ind[0]]];
			}
			for (std::size_t i = 1; i < nflip; i++) {
				//InfoMessage() << "toflip = " << toflip[i] << std::endl;
				new_prods *= mps_contraction(v, toflip[sorted_ind[i-1]] + 1, toflip[sorted_ind[i]]) * W_[toflip[sorted_ind[i]] % symperiod_][newconf[sorted_ind[i]]];
			}
			if (toflip[nflip - 1] < N_ - 1) {
				new_prods *= mps_contraction(v, toflip[sorted_ind[nflip - 1]] + 1, N_);
			}
			//InfoMessage() << "LogValDiff lookup ended " << std::log(new_prods.trace() / current_psi) << std::endl;
			return std::log(new_prods.trace() / current_psi);
		};

		// Derivative with full calculation
		VectorType DerLog(const std::vector<int> &v) override {
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
} // namespace netket

#endif
