AsymmetricNumeralSystemsToolkit
===============================

Testing various methods for choosing tANS entropy coding automata

ANS is new approach to entropy coding, which adds fractional bits into consideration into Huffman-like decoder, combining its speed with accuracy of arithmetic coding, like in implementation of Yann Collet ( https://github.com/Cyan4973/FiniteStateEntropy ). Another advantage in comparison to Huffman coding is that we choose the size of coding tables (L) here, corresponding to 2^depth of Huffman tree, and that there is no need to sort symbol probabilities (linear initialization).

The choice of such finite L state entropy coding automaton consists of:
- quantization of symbol probability distribution as p[s] ~ q[s]/L fractions (q is a natural number)
- spreading these symbols in range [0, L-1], such that symbol s appears q[s] times

This toolkit contains various choices of these functions, allows to test obtained compression rates, compare with Huffman. Currently it allows to choose betwen:

1) Symbol probability distributions: 
- init_binary(p, n): n binary variables (2^n alphabet size)
- init_power(rho,m): p[i] is proportional to rho^i
- init_rand_unif(m): simple random distribution;

2) quantizer (L=2^R):
- quantize_fast(R): first q[s] = round(p[s]*L), then modify q for the most probable symbol to fulfill nomalization (linear complexity),
- quantize_prec(R): find the quantization which minimizes sum_s (p[s] - q[s]/L)^2 / p[s] approximation of Kullback_Leibler distance (n log n complexity).

3) symbol spread:
- spread_fast(): basic Yann Collet's random spread - very fast,
- spread_prec():  very good using only q - wants to distributie symbols in 1/q[s] distance (still linear, but slower),
- spread_tuned(): uses both q and p - wants to get close to 1/x stationary probability of states (still linear, best compression rate).

For example for 4 symbol and L=16 states:

p: 0.04 0.16 0.16 0.64

q/L: 0.0625 0.1875 0.125 0.625

this quantization itself has dH/H ~ 0.011 rate penalty

spread_fast() gives 0233233133133133 and dH/H ~ 0.020

spread_prec() gives 3313233103332133 and dH/H ~ 0.015 - generally it is close to quantization dH/H

while spread_tuned(): 3233321333313310 and dH/H ~ 0.0046 - better than quantization dH/H due to using also p

Tuning shifts symbols right when q[s]/L > p[s] and left otherwise, getting better agreement and so compression rate. 

Some sources: [article](http://arxiv.org/abs/1311.2540), [slides](https://dl.dropboxusercontent.com/u/12405967/ANSsem.pdf), [discussion](http://encode.ru/threads/2013-Asymmetric-numeral-system-toolkit-and-fast-tuned-symbol-spread).

Feel free to add new probability distributions, better quantizers and spreads.

Jarek Duda, July 2014
