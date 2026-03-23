/*
 * Copyright Hamid Mohamadi.
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License. See the LICENSE accompanying this file
 * for the specific language governing permissions and limitations under
 * the License.
 */

// Adapted from ntHash https://github.com/bcgsc/ntHash

#ifndef BLOOMFILTER_HPP_
#define BLOOMFILTER_HPP_

#include <climits>
#include <cstdint>
#include <cstring>
#include <fstream>

#include "nthash.hpp"

/**
 * Count the number of set bits in a byte.
 *
 * @param x  Input byte.
 * @return Number of bits set to 1.
 */
inline unsigned popCount(uint8_t x) {
    return ((0x876543210 >>
             (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2))
           & 0xf;
}

class BloomFilter {
public:
    /**
     * Load a Bloom filter from an existing file.
     *
     * @param filterSize  Bit size of the Bloom filter.
     * @param hashNum  Number of hash functions.
     * @param kmerSize  K-mer size used to build the filter.
     * @param fPath  Path to the stored Bloom filter file.
     */
    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize, const char* fPath)
        : m_filter(nullptr), m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        size_t bytes = (m_size + CHAR_BIT - 1) / CHAR_BIT;
        m_filter = new unsigned char[bytes];
        std::ifstream file(fPath, std::ios::in | std::ios::binary);
        file.seekg(0, std::ios::beg);
        file.read(reinterpret_cast<char*>(m_filter), bytes);
        file.close();
    }

    /**
     * Construct an empty Bloom filter.
     *
     * @param filterSize  Bit size of the Bloom filter.
     * @param hashNum  Number of hash functions.
     * @param kmerSize  K-mer size used to build the filter.
     */
    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize)
        : m_filter(nullptr), m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize) {
        size_t bytes = (m_size + CHAR_BIT - 1) / CHAR_BIT;
        m_filter = new unsigned char[bytes]();
    }

    // Non-copyable
    BloomFilter(const BloomFilter&) = delete;
    BloomFilter& operator=(const BloomFilter&) = delete;

    ~BloomFilter() {
        delete[] m_filter;
    }

    /**
     * Insert a k-mer using precomputed hash values and report if it was new.
     *
     * @param hVal  Array of hash values for the k-mer.
     * @return True if the k-mer was inserted for the first time.
     */
    bool insert_make_change(const uint64_t* hVal) {
        bool change = false;
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            unsigned char mask = 1 << (CHAR_BIT - 1 - hLoc % CHAR_BIT);
            unsigned char oldByte = __sync_fetch_and_or(&m_filter[hLoc / CHAR_BIT], mask);
            if ((oldByte & mask) == 0) {
                change = true;
            }
        }
        return change;
    }

    /**
     * Insert a k-mer using precomputed hash values.
     *
     * @param hVal  Array of hash values for the k-mer.
     */
    void insert(const uint64_t* hVal) {
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            __sync_or_and_fetch(&m_filter[hLoc / CHAR_BIT],
                                1 << (CHAR_BIT - 1 - hLoc % CHAR_BIT));
        }
    }

    /**
     * Insert a k-mer by its sequence string.
     *
     * @param kmer  The k-mer sequence.
     */
    void insert(const char* kmer) {
        uint64_t hVal[m_hashNum];
        NTMC64(kmer, m_kmerSize, m_hashNum, hVal);
        insert(hVal);
    }

    /**
     * Query a k-mer using precomputed hash values.
     *
     * @param hVal  Array of hash values for the k-mer.
     * @return True if all hash positions are set.
     */
    bool contains(const uint64_t* hVal) const {
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t hLoc = hVal[i] % m_size;
            if ((m_filter[hLoc / CHAR_BIT] & (1 << (CHAR_BIT - 1 - hLoc % CHAR_BIT))) == 0) {
                return false;
            }
        }
        return true;
    }

    /**
     * Query a k-mer by its sequence string.
     *
     * @param kmer  The k-mer sequence.
     * @return True if the k-mer is in the Bloom filter.
     */
    bool contains(const char* kmer) const {
        uint64_t hVal[m_hashNum];
        NTMC64(kmer, m_kmerSize, m_hashNum, hVal);
        return contains(hVal);
    }

    /**
     * Store the Bloom filter to a binary file.
     *
     * @param fPath  Output file path.
     */
    void storeFilter(const char* fPath) const {
        size_t bytes = (m_size + CHAR_BIT - 1) / CHAR_BIT;
        std::ofstream file(fPath, std::ios::out | std::ios::binary);
        file.write(reinterpret_cast<char*>(m_filter), bytes);
        file.close();
    }

    /**
     * Compute the population count (number of set bits) in the filter.
     *
     * @return Number of set bits.
     */
    size_t get_pop() const {
        size_t bytes = (m_size + CHAR_BIT - 1) / CHAR_BIT;
        size_t popBF = 0;
        #pragma omp parallel for reduction(+:popBF)
        for (size_t i = 0; i < bytes; i++) {
            popBF += popCount(m_filter[i]);
        }
        return popBF;
    }

    unsigned getHashNum() const { return m_hashNum; }
    unsigned getKmerSize() const { return m_kmerSize; }
    unsigned char* getFilter() { return m_filter; }

private:
    unsigned char* m_filter;
    size_t m_size;
    unsigned m_hashNum;
    unsigned m_kmerSize;
};

#endif // BLOOMFILTER_HPP_