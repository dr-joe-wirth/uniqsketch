#ifndef UNIQSKETCHUTIL_HPP_
#define UNIQSKETCHUTIL_HPP_

#include <algorithm>
#include <cmath>
#include <cstring>
#include <iomanip>
#include <numeric>
#include <random>
#include <sstream>
#include <unordered_map>

#include "SequenceUtil.hpp"
#include "BloomFilter.hpp"
#include "ntHashIterator.hpp"

namespace opt {
unsigned numRef;                                        // number of references in the universe
unsigned threads(1);                                    // number of threads
unsigned bits(128);                                     // bits per element in Bloom filter
unsigned nhash1(5);                                     // number of hashes for distinct BF
unsigned nhash2(5);                                     // number of hashes for solid BF
unsigned kmerLen(81);                                   // length of k-mer
int sketchnum(100);                                     // unique k-mers to represent a reference
size_t m1;                                              // bit size of distinct BF
size_t m2;                                              // bit size of solid BF
size_t dbfSize(0);                                      // distinct BF size (number of elements)
size_t sbfSize(0);                                      // solid BF size (number of elements)
std::string outdir("outdir_uniqsketch");                // dir for universe unique k-mers
std::string refstat("db_uniq_count.tsv");               // unique k-mers stat filename
std::string outfile("sketch_uniq.tsv");                 // output sketch file name
unsigned kRange(5);                                     // k-mer spectrum range for entropy filter
double maxEntropy(12.0);                                // initial entropy for k={1,2,3}
double entropyThreshold(0.65);                          // entropy score rate threshold
}

using SketchHash = std::unordered_map<std::string, unsigned>;

/**
 * Insert k-mers from a FASTA reference into distinct and solid Bloom filters.
 *
 * @param fPath  Path to a reference file.
 * @param dbFilter  Distinct Bloom filter for tracking distinct k-mers.
 * @param sbFilter  Solid Bloom filter for tracking repeat k-mers.
 */
void loadBFfa(const std::string& fPath, BloomFilter& dbFilter, BloomFilter& sbFilter) {
    std::ifstream refFile(fPath);
    std::string seq;
    bool good = static_cast<bool>(std::getline(refFile, seq));
    while (good) {
        std::string faSeq;
        good = static_cast<bool>(std::getline(refFile, seq));
        while (good && seq[0] != '>') {
            faSeq += seq;
            good = static_cast<bool>(std::getline(refFile, seq));
        }
        ntHashIterator itr(faSeq, std::max(opt::nhash1, opt::nhash2), opt::kmerLen);
        while (itr != itr.end()) {
            if (!dbFilter.insert_make_change(*itr)) {
                sbFilter.insert(*itr);
            }
            ++itr;
        }
    }
}

/**
 * Scan all k-mers in a reference and record unique k-mers.
 *
 * @param fPath  Path to a reference file.
 * @param refId  Integer id assigned to the reference.
 * @param sbFilter  Solid Bloom filter for tracking repeat k-mers.
 * @param refCount  Vector of unique k-mer counts per reference.
 */
void checkRef(const std::string& fPath, int refId,
              const BloomFilter& sbFilter, std::vector<int>& refCount) {
    std::ifstream refFile(fPath);
    std::string sampleId = getBaseId(fPath);
    std::ofstream uniqFile(opt::outdir + "/" + sampleId + ".tsv");
    std::string seq;
    size_t countUnique = 0;

    bool good = static_cast<bool>(std::getline(refFile, seq));
    std::string headerSeq = seq.substr(1);
    std::replace(headerSeq.begin(), headerSeq.end(), ' ', '_');

    while (good) {
        std::string faSeq;
        good = static_cast<bool>(std::getline(refFile, seq));
        while (good && seq[0] != '>') {
            faSeq += seq;
            good = static_cast<bool>(std::getline(refFile, seq));
        }
        ntHashIterator itr(faSeq, opt::nhash2, opt::kmerLen);
        while (itr != itr.end()) {
            if (!sbFilter.contains(*itr)) {
                ++countUnique;
                size_t seqLoc = itr.get_pos();
                uniqFile << get_canonical(faSeq.substr(seqLoc, opt::kmerLen))
                         << "\t" << seqLoc << "\t" << headerSeq << "\n";
            }
            ++itr;
        }
        if (seq.size() > 1) {
            headerSeq = seq.substr(1);
            std::replace(headerSeq.begin(), headerSeq.end(), ' ', '_');
        }
    }
    uniqFile.close();
    refCount[refId] = countUnique;
}

/**
 * Identify all unique k-mers in the reference universe.
 *
 * @param refFiles  Vector of all reference file paths.
 */
void identifyUniqKmers(const std::vector<std::string>& refFiles) {
    std::cout << "k-mer length: " << opt::kmerLen << "\n"
              << "Number of distinct k-mers: " << opt::dbfSize << "\n"
              << "Number of repeat k-mers: " << opt::sbfSize << "\n"
              << "Number of unique k-mers: " << opt::dbfSize - opt::sbfSize << "\n";

    BloomFilter dbFilter(opt::m1, opt::nhash1, opt::kmerLen);
    BloomFilter sbFilter(opt::m2, opt::nhash2, opt::kmerLen);

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < refFiles.size(); ++i) {
        loadBFfa(refFiles[i], dbFilter, sbFilter);
    }

    std::vector<int> refCount(opt::numRef);
    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < refFiles.size(); ++i) {
        checkRef(refFiles[i], i, sbFilter, refCount);
    }

    std::ofstream tsvOut(opt::refstat);
    for (size_t i = 0; i < refCount.size(); i++) {
        tsvOut << opt::outdir << "/" << getBaseId(refFiles[i]) << ".tsv\t"
               << refCount[i] << "\t" << i << "\n";
    }
    tsvOut.close();

    std::cout << "Distinct BF actual fpr: " << std::setprecision(4) << std::fixed
              << std::pow(static_cast<double>(dbFilter.get_pop()) / opt::m1, opt::nhash1) << "\n"
              << "Solid BF actual fpr: " << std::setprecision(4) << std::fixed
              << std::pow(static_cast<double>(sbFilter.get_pop()) / opt::m2, opt::nhash2) << "\n";
}

/**
 * Digitize a DNA sequence into an integer representation.
 *
 * @param seq  DNA sequence to digitize.
 * @return Integer representation of the DNA sequence.
 */
unsigned digitize(const std::string& seq) {
    if (seq.length() > opt::kRange) {
        return 0;
    }
    uint64_t hVal = 0;
    for (unsigned i = 0; i < seq.length(); i++) {
        hVal = (hVal << 2) | b2f[static_cast<unsigned char>(seq[i])];
    }
    return static_cast<unsigned>(hVal);
}

/**
 * Check a signature for low-complexity content using aggregate Shannon entropy.
 *
 * @param signature  Signature sequence to check.
 * @return True if the signature is low-complexity, false otherwise.
 */
bool lowComplexity(const std::string& signature) {
    if (opt::kmerLen <= opt::kRange) {
        return true;
    }

    // Max spectrum size for kRange=5 is 4^5 = 1024
    static constexpr unsigned MAX_SPEC = 1024;
    unsigned kmerSpec[MAX_SPEC];

    double aggEntropy = 0.0;
    for (unsigned k = 1; k <= opt::kRange; k++) {
        unsigned specSize = 1u << (2 * k);
        unsigned mask = specSize - 1;
        std::memset(kmerSpec, 0, specSize * sizeof(unsigned));

        // Seed the rolling hash with the first k-mer
        unsigned hVal = 0;
        for (unsigned j = 0; j < k; j++) {
            hVal = (hVal << 2) | b2f[static_cast<unsigned char>(signature[j])];
        }
        kmerSpec[hVal]++;

        // Roll through remaining positions: shift out old base, shift in new base
        for (unsigned i = 1; i <= signature.length() - k; i++) {
            hVal = ((hVal << 2) | b2f[static_cast<unsigned char>(signature[i + k - 1])]) & mask;
            kmerSpec[hVal]++;
        }

        double entropy = 0.0;
        double denom = signature.size() - k + 1;
        for (unsigned i = 0; i < specSize; i++) {
            if (kmerSpec[i]) {
                double p = static_cast<double>(kmerSpec[i]) / denom;
                entropy -= p * std::log2(p);
            }
        }
        aggEntropy += entropy;
    }
    return (aggEntropy / opt::maxEntropy <= opt::entropyThreshold);
}

/**
 * Extract the unique k-mer set for a given reference.
 *
 * @param fPath  Path to a reference's candidate unique k-mers.
 * @param refId  Integer id assigned to the reference.
 * @param sketchFilter  Bloom filter for the unique sketch.
 * @param out  Output file stream for the uniqsketch.
 */
void getUniqSet(const std::string& fPath, unsigned refId,
                BloomFilter& sketchFilter, std::ofstream& out) {
    // Read all candidate lines
    std::vector<std::string> lines;
    std::ifstream sketchFile(fPath);
    std::string line;
    while (std::getline(sketchFile, line)) {
        lines.push_back(std::move(line));
    }
    sketchFile.close();

    if (lines.empty()) return;

    // Build a lightweight index: extract (contig, pos) for sorting without extra string copies
    size_t total = lines.size();
    std::vector<unsigned> indices(total);
    std::iota(indices.begin(), indices.end(), 0);

    // Parse position from each line for sorting (format: seq \t pos \t contig)
    // Use a fast scan: skip first field, read pos, read contig
    struct PosKey {
        unsigned pos;
        const char* contigStart;
        size_t contigLen;
    };
    std::vector<PosKey> keys(total);
    for (size_t i = 0; i < total; i++) {
        const char* s = lines[i].c_str();
        // skip seq field
        while (*s && *s != '\t') ++s;
        if (*s) ++s;
        // parse pos
        keys[i].pos = 0;
        while (*s >= '0' && *s <= '9') {
            keys[i].pos = keys[i].pos * 10 + (*s - '0');
            ++s;
        }
        if (*s) ++s;
        // contig starts here
        keys[i].contigStart = s;
        keys[i].contigLen = lines[i].c_str() + lines[i].size() - s;
    }

    std::sort(indices.begin(), indices.end(),
              [&keys](unsigned a, unsigned b) {
                  int cmp = std::strncmp(keys[a].contigStart, keys[b].contigStart,
                                         std::min(keys[a].contigLen, keys[b].contigLen));
                  if (cmp != 0) return cmp < 0;
                  if (keys[a].contigLen != keys[b].contigLen) return keys[a].contigLen < keys[b].contigLen;
                  return keys[a].pos < keys[b].pos;
              });

    // Select evenly spaced candidates, checking low-complexity lazily
    int count = 0;
    size_t needed = std::min(static_cast<size_t>(opt::sketchnum), total);
    double step = static_cast<double>(total) / needed;

    for (size_t s = 0; s < needed && count < opt::sketchnum; s++) {
        size_t slotStart = static_cast<size_t>(s * step);
        size_t slotEnd = (s + 1 < needed) ? static_cast<size_t>((s + 1) * step) : total;

        for (size_t i = slotStart; i < slotEnd; i++) {
            unsigned idx = indices[i];
            std::istringstream seqstm(lines[idx]);
            std::string useq, contig;
            unsigned pos;
            seqstm >> useq >> pos >> contig;

            if (lowComplexity(useq)) continue;

            ntHashIterator itr(useq, opt::nhash1, opt::kmerLen);
            bool found = false;
            while (itr != itr.end()) {
                if (sketchFilter.insert_make_change(*itr)) {
                    ++count;
                    found = true;
                    #pragma omp critical(out)
                    out << useq << "\t" << refId << "\t" << getBaseId(fPath)
                        << "\t" << opt::numRef << "\t" << contig << "\t" << pos << "\n";
                }
                ++itr;
            }
            if (found) break;
        }
    }
}

/**
 * Get the unique k-mer count from a stat line for sorting.
 *
 * @param sketchLine  Tab-separated line: path, count, refId.
 * @return Number of unique k-mers.
 */
inline int sketchRank(const std::string& sketchLine) {
    std::istringstream ss(sketchLine);
    std::string sketchPath;
    int sketchCount;
    ss >> sketchPath >> sketchCount;
    return sketchCount;
}

/**
 * Construct the unique sketch for the reference universe.
 *
 * @param statFile  Path to the unique k-mer stat file for all references.
 */
void buildSketch(const std::string& statFile) {
    std::ifstream tsvIn(statFile);
    std::vector<std::string> uniqStats;
    std::string line;
    while (std::getline(tsvIn, line)) {
        uniqStats.push_back(line);
    }

    std::sort(uniqStats.begin(), uniqStats.end(),
              [](const std::string& s1, const std::string& s2) {
                  return sketchRank(s1) < sketchRank(s2);
              });

    std::ofstream out(opt::outfile);
    BloomFilter sketchFilter(opt::bits * opt::sketchnum * uniqStats.size(),
                             opt::nhash1, opt::kmerLen);

    for (unsigned i = 4; i <= opt::kRange; i++) {
        opt::maxEntropy += std::log2(opt::kmerLen - i + 1);
    }

    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < uniqStats.size(); i++) {
        std::istringstream uniqstm(uniqStats[i]);
        std::string sketchPath;
        size_t sketchCount, sketchRefId;
        uniqstm >> sketchPath >> sketchCount >> sketchRefId;
        getUniqSet(sketchPath, sketchRefId, sketchFilter, out);
    }
    out.close();
}

#endif // UNIQSKETCHUTIL_HPP_