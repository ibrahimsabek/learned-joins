# Exotic Hashing

A header only cmake/c++ library that exposes various state-of-the-art
exotic hash functions, i.e., hash functions that possess at least one
of the following attributes:

* *perfect*: produce **no** collisions, i.e., is injective
* *minimal perfect*: perfect and output range is [0, N] for N-element keyset.
* *order preserving*: keys retain a certain order after being hashed
* *monotone*: k1 <= k2 -->  h(k1) <= h(k2)
