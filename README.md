AsymmetricNumeralSystemsToolkit
===============================

Testing various methods for choosing tANS entropy coding automata

ANS is new approach to entropy coding, which adds fractional bits into consideration into Huffman-like decoder, combining its speed with accuracy of arithmetic coding, like in [implementation of Yann Collet](https://github.com/Cyan4973/FiniteStateEntropy). Another advantage in comparison to Huffman coding is that we choose the size of coding tables (L) here, corresponding to 2^depth of Huffman tree, and that there is no need to sort symbol probabilities (linear initialization).

The choice of such finite L state entropy coding automaton consists of 2-3 parts (which generally can be merged):
- quantization of symbol probability distribution as p[s] ~ q[s]/L fractions (q is a natural number)
- spreading these symbols in range [0, L-1], such that symbol s appears q[s] times
- eventual scrambler to simultaneously encrypt encoded message: disturb the found spread in a pseudoranom way accordingly to cryptographic key.

This toolkit contains various choices of these functions, allows to test obtained compression rates, compare with Huffman. Currently it allows to choose betwen:

1) Symbol probability distributions: 
- init_binary(p, n): n binary variables (2^n alphabet size)
- init_power(rho,m): p[i] is proportional to rho^i
- init_rand_unif(m): simple random distribution;

2) quantizer (L=2^R):
- quantize_fast(R): first q[s] = round(p[s]*L), then modify q for the most probable symbol to fulfill nomalization (linear complexity),
- quantize_prec(R): find the quantization which minimizes sum_s (p[s] - q[s]/L)^2 / p[s] approximation of Kullback_Leibler distance (n log n complexity).

3) symbol spread:
- spread_range_i(), spread_range_d(), : [Huffman coding can be seen as a special case of tANS](http://fastcompression.blogspot.fr/2014/01/huffman-comparison-with-fse.html): if spreading in size 2^k ranges - these spreads take it to general frequencies (in increasing or decreasing order),
- spread_fast(): basic Yann Collet's random spread - very fast,
- spread_prec():  very good using only q - wants to distributie symbols in 1/q[s] distance (still linear, but slower),
- spread_tuned(): uses both q and p - wants to get close to 1/x stationary probability of states (still linear, best compression rate),
- spread_tuned_s(): uses sort instead of buckets - can be slighlty better, but slower,
- spread_tuned_p(): uses 1/i approximation of 1/(p*ln(1+1/i)) formula for tuned spread,
- spread_uABS(): available arithmetic coding/decoding formulas for binary alphabet, e.g. used in [Matt Mahoney's fpaqc](http://www.mattmahoney.net/dc/).
 
4) scramblers:
- scrambler0(): for each 2i-1 and 2i positions, with some probability switch symbols - accordingly to PRNG initialized with cryptographic key (the current version assumes at most 256 size alphabet),
- scrambler1(): takes blocks of four symbols and randomly rotate them cyclically (also 256 size alphabet limit) - faster and stronger disturbtion.
 
For example for 4 symbol and L=16 states: p=(0.04 0.16 0.16 0.64), q/L=(0.0625 0.1875 0.125 0.625). 
<table>
  <tr>
    <th>method</th><th>symbol spread</th><th>dH/H rate loss</th><th>comment</th>
  </tr>
  <tr>
    <th> - </th><th> - </th><th>~0.011</th><th>penalty of quantizer itself</th>
  </tr>
  <tr>
    <th> Huffman </th><th> 0011222233333333 </th><th>~0.080</th><th>would give Huffman decoder</th>
  </tr>
  <tr>
    <th>spread_range_i()</th><th>0111223333333333</th><th>~0.059</th><th> analogous to Huffman </th>
  </tr>
  <tr>
    <th>spread_range_d()</th><th>3333333333221110</th><th>~0.022</th><th> decreasing order </th>
  </tr>
  <tr>
    <th>spread_fast()</th><th>0233233133133133</th><th>~0.020</th><th> fast </th>
  </tr>
  <tr>
    <th>spread_prec()</th><th>3313233103332133</th><th>~0.015</th><th>generally close to quantization dH/H</th>
  </tr>
    <tr>
    <th>spread_tuned()</th><th>3233321333313310</th><th>~0.0046</th><th>better than quantization dH/H due to using also p</th>
  </tr>
  <tr>
    <th>spread_tuned_s()</th><th>2333312333313310</th><th>~0.0040</th><th>L log L complexity (sort)</th>
  </tr>
  <tr>
    <th>spread_tuned_p()</th><th>2331233330133331</th><th>~0.0058</th><th>testing 1/(p ln(1+1/i)) ~ i/p approximation</th>
  </tr>
</table>

Tuning shifts symbols right when q[s]/L > p[s] and left otherwise, getting better agreement and so compression rate. 

Some sources: [article](http://arxiv.org/abs/1311.2540), [slides](https://dl.dropboxusercontent.com/u/12405967/ANSsem.pdf), [discussion](http://encode.ru/threads/2013-Asymmetric-numeral-system-toolkit-and-fast-tuned-symbol-spread), [list of  implementations of ANS](http://encode.ru/threads/2078-List-of-Asymmetric-Numeral-Systems-implementations).

Feel free to add new probability distributions, better quantizers and spreads.

Jarek Duda, July 2014
