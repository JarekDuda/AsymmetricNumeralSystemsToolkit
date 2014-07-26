AsymmetricNumeralSystemsToolkit
===============================

Testing various methods for choosing tANS entropy coding automata

ANS is new approach to entropy coding, which adds fractional bits into consideration into Huffman-like decoder, combining its speed with accuracy of arithmetic coding, like in implementation of Yann Collet: https://github.com/Cyan4973/FiniteStateEntropy
The choice of such finite L state entropy coding automaton consists of:
- quantization of symbol probability distribution as p[s] ~ q[s]/L fractions (q is a natural number)
- speading these symbols in range [0, L-1], such that symbol s appears q[s] times
This toolkit contains various choices of these functions, allows to test obtained compression rates, compare with Huffman. Currently it allows to choose betwen:

1) Symobol probability distributions: 
- init_binary(p, n): n binary variables (2^n alphabet size)
- init_power(rho,m): p[i] is proportional to rho^i
- init_rand_unif(m): simple random distribution;

2) quantizer (L=2^R):
- quantize_fast(R): first q[s] = round(p[s]*L), then modify q for the most probable symbol to fulfill nomalization (linear complexity),
- quantize_prec(R): find the quantization which minimizes sum_s (p[s] - q[s]/L)^2 / p[s] approximation of Kullback_Leibler distance (n log n).

3) symbol spread:
- spread_fast(): basic Yann Collet's random spread - very fast,
- spread_prec():  very good using only q - wants to distributie symbols in 1/q[s] distance (still linear, but slower),
- spread_tuned(): uses both q and p - wants to get close to 1/x stationary probability of states (still linear, best compression rate).

Some sources:

arxiv: http://arxiv.org/abs/1311.2540

slides: https://dl.dropboxusercontent.com/u/12405967/ANSsem.pdf

Feel free to add new probability distributions, better quantizers and spreads.

Jarek Duda, July 2014
