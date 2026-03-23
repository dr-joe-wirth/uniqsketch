#ifndef CLUSTERUTIL_HPP_
#define CLUSTERUTIL_HPP_

#include <climits>
#include <cstring>
#include <functional>
#include <iostream>
#include <fstream>
#include <numeric>
#include <string>
#include <unordered_map>
#include <vector>

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "kseq_gz.h"

/**
 * Load all k-mer hashes from a FASTA/FASTA.gz reference into a flat vector.
 *
 * @param fPath  Path to a reference file (plain or gzipped).
 * @param nhash  Number of hash functions.
 * @param kmerLen  K-mer length.
 * @return Flattened vector of hash values (nhash values per k-mer).
 */
inline std::vector<uint64_t> clusterLoadKmerHashes(const std::string& fPath,
                                                   unsigned nhash, unsigned kmerLen) {
    std::vector<uint64_t> hashes;
    gzFile fp = gzopen(fPath.c_str(), "r");
    kseq_t* ks = kseq_init(fp);
    while (kseq_read(ks) >= 0) {
        std::string seq(ks->seq.s, ks->seq.l);
        ntHashIterator itr(seq, nhash, kmerLen);
        while (itr != itr.end()) {
            const uint64_t* hVec = *itr;
            for (unsigned h = 0; h < nhash; h++) {
                hashes.push_back(hVec[h]);
            }
            ++itr;
        }
    }
    kseq_destroy(ks);
    gzclose(fp);
    return hashes;
}

/**
 * Cluster references by pairwise unique k-mer count using single-linkage clustering.
 * References with fewer than `threshold` unique k-mers between them are merged.
 *
 * @param refFiles  All reference file paths.
 * @param threshold  Minimum unique k-mers to consider two references as distinct.
 * @param kmerLen  K-mer length.
 * @param nhash  Number of hash functions for Bloom filter.
 * @param bits  Bits per element in Bloom filter.
 * @param dbfSize  Approximate genome size for Bloom filter sizing.
 * @param clusterFile  Output TSV path mapping each reference to its cluster representative.
 * @return Vector of representative file paths (one per cluster).
 */
inline std::vector<std::string> clusterReferences(
        const std::vector<std::string>& refFiles,
        unsigned threshold,
        unsigned kmerLen,
        unsigned nhash,
        unsigned bits,
        size_t dbfSize,
        const std::string& clusterFile) {
    size_t n = refFiles.size();
    if (n <= 1) return refFiles;

    std::cerr << "Clustering " << n << " references (threshold=" << threshold
              << " unique k-mers, k=" << kmerLen << ")...\n";

    // Pre-compute k-mer hashes for all references
    std::vector<std::vector<uint64_t>> hashSets(n);
    std::vector<std::string> names(n);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < n; i++) {
        hashSets[i] = clusterLoadKmerHashes(refFiles[i], nhash, kmerLen);
        names[i] = getBaseId(refFiles[i]);
    }

    // Union-Find for single-linkage clustering
    std::vector<unsigned> parent(n);
    std::iota(parent.begin(), parent.end(), 0);

    std::function<unsigned(unsigned)> find = [&](unsigned x) -> unsigned {
        while (parent[x] != x) {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        return x;
    };

    size_t bfBits = bits * dbfSize;

    #pragma omp parallel
    {
        BloomFilter dbFilter(bfBits, nhash, kmerLen);
        unsigned char* filterPtr = dbFilter.getFilter();
        size_t filterBytes = (bfBits + CHAR_BIT - 1) / CHAR_BIT;

        #pragma omp for schedule(dynamic)
        for (unsigned i = 0; i < n; i++) {
            std::memset(filterPtr, 0, filterBytes);
            for (size_t h = 0; h < hashSets[i].size(); h += nhash) {
                dbFilter.insert(&hashSets[i][h]);
            }

            for (unsigned j = i + 1; j < n; j++) {
                // Count how many k-mers in j are NOT in i's BF
                size_t countUnique = 0;
                size_t countAll = hashSets[j].size() / nhash;
                for (size_t h = 0; h < hashSets[j].size(); h += nhash) {
                    if (!dbFilter.contains(&hashSets[j][h])) {
                        ++countUnique;
                    }
                }

                if (countUnique < threshold) {
                    #pragma omp critical(union_op)
                    {
                        unsigned ri = find(i), rj = find(j);
                        if (ri != rj) parent[ri] = rj;
                    }
                }
            }
        }
    }

    // Build clusters
    std::unordered_map<unsigned, std::vector<unsigned>> clusters;
    for (unsigned i = 0; i < n; i++) {
        clusters[find(i)].push_back(i);
    }

    // Pick representative per cluster: member with the most k-mers
    std::vector<std::string> representatives;
    std::ofstream clusterOut(clusterFile);
    clusterOut << "cluster_id\trepresentative\tmember\tnum_kmers\n";

    unsigned clusterId = 0;
    for (auto& [root, members] : clusters) {
        unsigned repIdx = members[0];
        size_t maxKmers = hashSets[repIdx].size() / nhash;
        for (unsigned m : members) {
            size_t nk = hashSets[m].size() / nhash;
            if (nk > maxKmers) {
                maxKmers = nk;
                repIdx = m;
            }
        }
        representatives.push_back(refFiles[repIdx]);

        for (unsigned m : members) {
            clusterOut << clusterId << "\t" << names[repIdx] << "\t" << names[m]
                       << "\t" << hashSets[m].size() / nhash << "\n";
        }
        ++clusterId;
    }
    clusterOut.close();

    std::cerr << "Clustering complete: " << n << " references -> "
              << representatives.size() << " clusters\n";

    return representatives;
}

#endif // CLUSTERUTIL_HPP_