/*
 * Copyright Hamid Mohamadi.
 * SPDX-License-Identifier: MIT
 *
 * Licensed under the MIT License. See the LICENSE accompanying this file
 * for the specific language governing permissions and limitations under
 * the License.
 */

// Adapted from ntHash https://github.com/bcgsc/ntHash

#ifndef NTHASH_ITERATOR_HPP_
#define NTHASH_ITERATOR_HPP_

#include <cstdint>
#include <limits>
#include <string>

#include "nthash.hpp"

/**
 * Iterator over canonical ntHash values for k-mers in a DNA sequence.
 * Efficiently computes hash values for successive k-mers using rolling hash.
 */
class ntHashIterator {
public:
    /** Default constructor — creates an end-of-range sentinel. */
    ntHashIterator()
        : m_hVec(nullptr)
        , m_pos(std::numeric_limits<size_t>::max()) {}

    /**
     * Construct an iterator over k-mers in a sequence.
     *
     * @param seq  DNA sequence to hash.
     * @param h  Number of hash functions.
     * @param k  K-mer size.
     */
    ntHashIterator(const std::string& seq, unsigned h, unsigned k)
        : m_seq(seq), m_h(h), m_k(k)
        , m_hVec(new uint64_t[h])
        , m_pos(0), m_fhVal(0), m_rhVal(0) {
        init();
    }

    /** Copy constructor. */
    ntHashIterator(const ntHashIterator& other)
        : m_seq(other.m_seq), m_h(other.m_h), m_k(other.m_k)
        , m_hVec(new uint64_t[other.m_h])
        , m_pos(other.m_pos), m_fhVal(other.m_fhVal), m_rhVal(other.m_rhVal) {
        for (unsigned i = 0; i < m_h; i++) {
            m_hVec[i] = other.m_hVec[i];
        }
    }

    // Non-assignable (would need proper swap semantics)
    ntHashIterator& operator=(const ntHashIterator&) = delete;

    ~ntHashIterator() {
        delete[] m_hVec;
    }

    /** Get the current position in the sequence. */
    size_t get_pos() const { return m_pos; }

    /** Dereference — get pointer to hash values for current k-mer. */
    const uint64_t* operator*() const { return m_hVec; }

    bool operator==(const ntHashIterator& it) const { return m_pos == it.m_pos; }
    bool operator!=(const ntHashIterator& it) const { return m_pos != it.m_pos; }

    /** Pre-increment — advance to the next valid k-mer. */
    ntHashIterator& operator++() {
        next();
        return *this;
    }

    /** Sentinel iterator marking end of range. */
    static ntHashIterator end() { return ntHashIterator(); }

private:
    /** Initialize to the first valid k-mer position. */
    void init() {
        if (m_k > m_seq.length()) {
            m_pos = std::numeric_limits<size_t>::max();
            return;
        }
        unsigned locN = 0;
        while (m_pos < m_seq.length() - m_k + 1 &&
               !NTMC64(m_seq.data() + m_pos, m_k, m_h, m_fhVal, m_rhVal, locN, m_hVec)) {
            m_pos += locN + 1;
        }
        if (m_pos >= m_seq.length() - m_k + 1) {
            m_pos = std::numeric_limits<size_t>::max();
        }
    }

    /** Advance to the next valid k-mer. */
    void next() {
        ++m_pos;
        if (m_pos >= m_seq.length() - m_k + 1) {
            m_pos = std::numeric_limits<size_t>::max();
            return;
        }
        if (seedTab[static_cast<unsigned char>(m_seq.at(m_pos + m_k - 1))] == seedN) {
            m_pos += m_k;
            init();
        } else {
            NTMC64(m_seq.at(m_pos - 1), m_seq.at(m_pos - 1 + m_k),
                   m_k, m_h, m_fhVal, m_rhVal, m_hVec);
        }
    }

    std::string m_seq;      // DNA sequence
    unsigned m_h;           // number of hash functions
    unsigned m_k;           // k-mer size
    uint64_t* m_hVec;       // hash values array
    size_t m_pos;           // current k-mer position
    uint64_t m_fhVal;       // forward-strand hash value
    uint64_t m_rhVal;       // reverse-complement hash value
};

#endif // NTHASH_ITERATOR_HPP_