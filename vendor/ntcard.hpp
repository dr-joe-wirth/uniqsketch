/*
 * Copyright Hamid Mohamadi.
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License. See the LICENSE accompanying this file
 * for the specific language governing permissions and limitations under
 * the License.
 */

// Adapted from ntCard https://github.com/bcgsc/ntCard

#ifndef NTCARD_HPP_
#define NTCARD_HPP_

#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include "ntHashIterator.hpp"
#include "kseq_gz.h"

#ifdef _OPENMP
#include <omp.h>
#endif

namespace opt {
size_t rBuck;
unsigned rBits(27);
unsigned sBits(11);
unsigned sMask;
unsigned covMax(65536);
size_t nSamp(2);
size_t nK(1);
bool samH(true);
}

/**
 * Classify a hash value into a sampling bucket and increment the counter.
 *
 * @param hVal  Hash value to classify.
 * @param t_Counter  Counter array for sampling buckets.
 */
inline void ntComp(uint64_t hVal, uint16_t* t_Counter) {
    uint64_t indBit = opt::nSamp;
    if (hVal >> (64 - opt::sBits) == 1) {
        indBit = 0;
    }
    if (hVal >> (64 - opt::sBits) == opt::sMask) {
        indBit = 1;
    }
    if (indBit < opt::nSamp) {
        size_t shVal = hVal & (opt::rBuck - 1);
        #pragma omp atomic
        ++t_Counter[indBit * opt::rBuck + shVal];
    }
}

/**
 * Hash all k-mers in a sequence and update sampling counters.
 *
 * @param seq  DNA sequence.
 * @param kList  List of k-mer sizes.
 * @param t_Counter  Counter array for sampling buckets.
 * @param totKmer  Per-k total k-mer counts (output).
 */
inline void ntRead(const std::string& seq, const std::vector<unsigned>& kList,
                   uint16_t* t_Counter, size_t totKmer[]) {
    for (unsigned k = 0; k < kList.size(); k++) {
        ntHashIterator itr(seq, 1, kList[k]);
        while (itr != itr.end()) {
            ntComp((*itr)[0], t_Counter + k * opt::nSamp * opt::rBuck);
            ++itr;
            ++totKmer[k];
        }
    }
}

/**
 * Compute cardinality estimates from sampling counters.
 *
 * @param t_Counter  Counter array from sampling.
 * @param F0Mean  Estimated number of distinct k-mers (output).
 * @param fMean  Estimated k-mer frequency histogram (output).
 */
void compEst(const uint16_t* t_Counter, double& F0Mean, std::vector<double>& fMean) {
    std::vector<std::vector<unsigned>> p(opt::nSamp, std::vector<unsigned>(opt::covMax, 0));
    for (size_t i = 0; i < opt::nSamp; i++) {
        for (size_t j = 0; j < opt::rBuck; j++) {
            ++p[i][t_Counter[i * opt::rBuck + j]];
        }
    }

    std::vector<double> pMean(opt::covMax, 0.0);
    for (size_t i = 0; i < opt::covMax; i++) {
        for (size_t j = 0; j < opt::nSamp; j++) {
            pMean[i] += p[j][i];
        }
        pMean[i] /= static_cast<double>(opt::nSamp);
    }

    F0Mean = static_cast<double>(
        (opt::rBits * std::log(2) - std::log(pMean[0])) *
        (static_cast<size_t>(1) << (opt::sBits + opt::rBits)));

    fMean[1] = -1.0 * pMean[1] / (pMean[0] * (std::log(pMean[0]) - opt::rBits * std::log(2)));
    for (size_t i = 2; i < opt::covMax; i++) {
        double sum = 0.0;
        for (size_t j = 1; j < i; j++) {
            sum += j * pMean[i - j] * fMean[j];
        }
        fMean[i] = -1.0 * pMean[i] / (pMean[0] * (std::log(pMean[0]) - opt::rBits * std::log(2)))
                    - sum / (i * pMean[0]);
    }
    for (size_t i = 1; i < opt::covMax; i++) {
        fMean[i] = std::abs(static_cast<ssize_t>(fMean[i] * F0Mean));
    }
}

/**
 * Estimate k-mer cardinality from input sequence files.
 *
 * @param kmerNum  Total number of k-mers processed (output).
 * @param dbsize  Estimated number of distinct k-mers (output).
 * @param sbsize  Estimated number of solid k-mers (output).
 * @param kmer_len  K-mer length.
 * @param num_thread  Number of threads.
 * @param inFiles  Input sequence file paths.
 * @return Always returns false.
 */
bool getCardinality(size_t& kmerNum, size_t& dbsize, size_t& sbsize,
                    unsigned kmer_len, unsigned num_thread,
                    const std::vector<std::string>& inFiles) {
    std::vector<unsigned> kList = {kmer_len};
    std::vector<size_t> totalKmers(kList.size(), 0);
    opt::rBuck = static_cast<size_t>(1) << opt::rBits;
    opt::sMask = (static_cast<size_t>(1) << (opt::sBits - 1)) - 1;
    uint16_t* t_Counter = new uint16_t[opt::nK * opt::nSamp * opt::rBuck]();

#ifdef _OPENMP
    omp_set_num_threads(num_thread);
#endif

    #pragma omp parallel for schedule(dynamic)
    for (unsigned file_i = 0; file_i < inFiles.size(); ++file_i) {
        size_t totKmer[kList.size()];
        for (unsigned k = 0; k < kList.size(); k++) {
            totKmer[k] = 0;
        }
        gzFile fp = gzopen(inFiles[file_i].c_str(), "r");
        kseq_t* record = kseq_init(fp);
        while (kseq_read(record) >= 0) {
            std::string seq(record->seq.s, record->seq.l);
            ntRead(seq, kList, t_Counter, totKmer);
        }
        kseq_destroy(record);
        gzclose(fp);
        for (unsigned k = 0; k < kList.size(); k++) {
            #pragma omp atomic
            totalKmers[k] += totKmer[k];
        }
    }

    double F0Mean = 0.0;
    std::vector<double> fMean(opt::covMax, 0.0);
    compEst(t_Counter, F0Mean, fMean);

    kmerNum = totalKmers[0];
    dbsize = static_cast<size_t>(F0Mean);
    sbsize = dbsize - static_cast<size_t>(fMean[1]);

    delete[] t_Counter;
    return false;
}

#endif // NTCARD_HPP_