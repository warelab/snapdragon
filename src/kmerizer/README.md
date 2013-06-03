Kmerizer
========
Kmerizer is a library that implements functions for counting k-mers. It counts k-mers using a hybrid approach that initially counts k-mers by sorting them, and afterward by consulting a read-only lookup table of the most common k-mers. Rare k-mers that are not found in the lookup table are counted by sorting them. Each batch of distinct kmers and their counts are stored using compressed bitvectors. Without this compression, the largest component of the run time would be the IO required for output files. Therefore, reducing the size of the output files will have a huge impact on the overall performance of the software.

Overview:
K-mers from an input sequence are packed (2 bits per nucleotide) into an 8 bit bin and up to 8 additional 64 bit words. Each k-mer is optionally canonicalized by comparing it to its reverse complement and selecting the minimum. Memory is allocated in advance for holding the bit packed k-mers. For parallelization, the k-mers are hashed into one of 256 bins based on the first 8 bits. If one of the bins fills before the input sequences have been processed, the k-mers in each of the 256 buffers are checked against the lookup table, and if at enough rare k-mers remain, they are counted, compressed and written to disk. Then, memory is reused for loading more k-mers. After all the sequences have been read in, another round of serialization happens for the last batch of k-mers. Then, if necessary, batches are merged with the in memory lookup table of common kmers.

Counting:
K-mers are counted by sorting them in place and then iterating over them to identify and count the distinct k-mers. When k>36 multiple words are required to hold the packed sequence, and a custom multi-word comparison function is used. We have arbitrarily coded 8 such functions thereby limiting the maximum kmer to 260 nucleotides.

Compression:
The sorted distinct k-mers each occupy 64*ceil((k-4)/32) bits. The first 8 bits are implicitly stored in the bin number. We reduce the overall run time significantly by spending some CPU cycles to create a bit-sliced bitmap index of compressed bitvectors (one per bit position.) The high-order bits compress extremely well, while low order bits only require a small amount of overhead. Similarly, the k-mer counts are converted to a range encoded bitmap index with one compressed bitvector for each k-mer frequency. In this index, bitvector f marks the distinct k-mers that occur <= f times. The bitvector caches the number of set bits, so it is easy to calculate a histogram of k-mer frequencies.


