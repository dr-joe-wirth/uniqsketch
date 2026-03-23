/*
 * Shared kseq.h instantiation for gzFile.
 * Include this header wherever you need kseq_t / kseq_init / kseq_read with gzFile.
 */

#ifndef KSEQ_GZ_H_
#define KSEQ_GZ_H_

#include <zlib.h>
#include "kseq.h"

KSEQ_INIT(gzFile, gzread)

#endif // KSEQ_GZ_H_
