#include <iostream>
#include <cstdio>
#include <string>
#include <cassert>

#include "BloomFilter.hpp"

namespace opt {
size_t size;
unsigned hash;
unsigned k;
}

void initialize() {
    opt::size = 65536;
    opt::hash = 3;
    opt::k = 32;
}

void testInsertContain() {
    std::cerr << "START: Test insert, contain... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    uint64_t hVec1[] = {1, 3000, 50000000};
    uint64_t hVec2[] = {2000000, 700000000, 4000000000};

    bFilter.insert(hVec1);
    assert(bFilter.contains(hVec1) == true);
    assert(bFilter.contains(hVec2) == false);

    std::string kmer1("GCCTAGCTAGCTTTTAGCTGGGATTTTTT");
    std::string kmer2("ATTTAGCTAGCGCTAAGCTGGGACCCTGC");

    bFilter.insert(kmer1.c_str());
    assert(bFilter.contains(kmer1.c_str()) == true);
    assert(bFilter.contains(kmer2.c_str()) == false);

    std::cerr << "PASSED: Test insert, contain\n";
}

void testInsertMakeChange() {
    std::cerr << "START: Test insert_make_change... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    uint64_t hVec1[] = {34987645632, 89, 37692833};
    uint64_t hVec2[] = {34359738368, 4398046511238, 678};

    bFilter.insert(hVec1);
    assert(bFilter.insert_make_change(hVec1) == false);
    assert(bFilter.insert_make_change(hVec2) == true);

    std::cerr << "PASSED: Test insert_make_change\n";
}

// Test that an empty Bloom filter contains nothing
void testEmptyFilter() {
    std::cerr << "START: Test empty filter... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    uint64_t hVec[] = {42, 100, 999};
    assert(bFilter.contains(hVec) == false);
    assert(bFilter.get_pop() == 0);

    std::cerr << "PASSED: Test empty filter\n";
}

// Test get_pop returns correct population count
void testGetPop() {
    std::cerr << "START: Test get_pop... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    assert(bFilter.get_pop() == 0);

    uint64_t hVec1[] = {1, 3000, 50000000};
    bFilter.insert(hVec1);
    size_t pop1 = bFilter.get_pop();
    assert(pop1 > 0);
    // With 3 hash functions, at most 3 bits set (could be fewer if collisions)
    assert(pop1 <= opt::hash);

    // Inserting the same element again should not change population
    bFilter.insert(hVec1);
    assert(bFilter.get_pop() == pop1);

    // Inserting a new element should increase population
    uint64_t hVec2[] = {2000000, 700000000, 4000000000};
    bFilter.insert(hVec2);
    assert(bFilter.get_pop() >= pop1);

    std::cerr << "PASSED: Test get_pop\n";
}

// Test store and load roundtrip
void testStoreLoad() {
    std::cerr << "START: Test store/load... \t";

    BloomFilter bFilter(opt::size, opt::hash, opt::k);
    uint64_t hVec1[] = {1, 3000, 50000000};
    uint64_t hVec2[] = {2000000, 700000000, 4000000000};

    bFilter.insert(hVec1);
    bFilter.insert(hVec2);

    const char* tmpPath = "test_bf_roundtrip.bin";
    bFilter.storeFilter(tmpPath);

    // Load into a new filter and verify contents match
    BloomFilter loaded(opt::size, opt::hash, opt::k, tmpPath);
    assert(loaded.contains(hVec1) == true);
    assert(loaded.contains(hVec2) == true);

    uint64_t hVec3[] = {999999, 888888, 777777};
    assert(loaded.contains(hVec3) == false);

    assert(loaded.get_pop() == bFilter.get_pop());

    std::remove(tmpPath);

    std::cerr << "PASSED: Test store/load\n";
}

int main() {
    initialize();
    testInsertContain();
    testInsertMakeChange();
    testEmptyFilter();
    testGetPop();
    testStoreLoad();

    std::cerr << "BloomFilter: All tests PASSED!\n\n";
    return 0;
}
