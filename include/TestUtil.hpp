#ifndef TESTUTIL_HPP_
#define TESTUTIL_HPP_

#include <algorithm>
#include <fstream>
#include <iterator>
#include <string>

/**
 * Compare two files byte-by-byte for equality.
 *
 * @param path1  Path to the first file.
 * @param path2  Path to the second file.
 * @return True if files are identical, false otherwise.
 */
bool compareFiles(const std::string& path1, const std::string& path2) {
    std::ifstream file1(path1);
    std::ifstream file2(path2);
    return std::equal(std::istreambuf_iterator<char>(file1),
                      std::istreambuf_iterator<char>(),
                      std::istreambuf_iterator<char>(file2));
}

#endif // TESTUTIL_HPP_
