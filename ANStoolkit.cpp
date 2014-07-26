// ANS toolkit     Jarek Duda, July 2014
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include <algorithm>
#include <vector>
#include <utility>
#include <random>
using namespace std;

typedef uint16_t    avar;   // change to 32 if alphabet larger than 65k
typedef uint16_t    tvar;   // change to 32 if table larger than 65k
typedef double      prec;   // change to float to reduce memory need 
const prec acc=1E-10;       // accuracy for hANS

struct ANS {
  avar m;               // alphabet size 
  tvar L;               // table size L - assumed to be a power of 2
  prec *p=NULL;         // probability table
  tvar *q=NULL;          // probabilities quantized to 1/L
  avar *s=NULL;         // symbol spread
  
  prec h, hq, hhuff, hANS;     // entropy of p, of q, for Huffman, and ANS
  int huff_depth;         // depth of Huffman tree
  prec *sp=NULL;          // stationary state probability distribution (after x -> floor(x/2) reductions)
  ~ANS(){delete[] p; delete[] s; delete[] q; delete[] sp;}
 void printp(ostream &out=cout) 
 {out<<"p: "; for(int i=0;i<m;i++) out<<p[i]<<" "; out<<endl;}
 void printq(ostream &out=cout) {
     out<<"q/L: ";for(int i=0;i<m;i++) out<<(float)q[i]/L<<" "; out<<endl;
     out<<"q: ";for(int i=0;i<m;i++) out<<q[i]<<" "; out<<endl;
 }
 void prints(ostream &out=cout) 
 {out<<"s: "; for(int i=0;i<L;i++) {out<<s[i];if(m>10)out<<" ";}; out<<endl;}
 void normalize(){ 
     prec sum=0; for(int i=0;i<m;i++) sum+=p[i]; for(int i=0;i<m;i++) p[i]/=sum;
 }
  void calc_h(){                             // calculate entropies
    h=0; for(int i=0;i<m;i++) h-=p[i]*log(p[i]); h/=log(2);
    hq=log(L); for(int i=0;i<m;i++) hq-=p[i]*log(q[i]); hq/=log(2);
    hhuff=0; vector<pair<prec,int>> v;               // calculate hhuff and huff_depth
    for(int i=0;i<m;i++)v.push_back({-p[i],0});      // to search for minimum not maximum
    make_heap(v.begin(),v.end());           
    for(int i=1;i<m;i++){
        auto t=v.front();pop_heap(v.begin(),v.end());v.pop_back();
        auto tt=v.front();pop_heap(v.begin(),v.end()); v.pop_back();
        hhuff -= t.first + tt.first;  
        v.push_back({t.first + tt.first, max(t.second,tt.second)+1}); 
        push_heap (v.begin(),v.end());
    }
    huff_depth=v.front().second;
 }
 
 // ----------------- QUANTIZERS ------------------------
 void quantize_fast(tvar LL=10){   // fast quantization - correction applied to largest probability symbol
     L = 1 << LL;                   // sometimes does not work!
     delete[] q; q = new tvar[m];
     tvar used=0, maxp;
     prec maxv=0;
     for(int i=0;i<m;i++){ 
         q[i]=round(L*p[i]); if(!q[i])q[i]++; used+=q[i];
         if(p[i]>maxv){maxv = p[i]; maxp = i;}
     }
     q[maxp]+=L-used;
 }
 void quantize_prec(tvar LL=10){   // precise quantization - minimizing dH~sum (p_i-q_i)^2/p_i
     L = 1 << LL;                      // can be speeded up by rewriting heap
     delete[] q; q = new tvar[m];                
     prec *pL = new prec[m], *ip = new prec[m];
     prec *cc = new prec[m];        // current approximation cost for symbols
     tvar used = 0;                   // currently used 
     for(int i=0;i<m;i++){
         pL[i]=p[i]*L; ip[i]=1/p[i];               // tabled for speedup
         q[i]=round(pL[i]); if(!q[i])q[i]++; used+=q[i];
         cc[i] = pow(pL[i]-q[i],2)*ip[i];                //L^2 (p_i-q_i)^2/p_i 
     }
     if(used!=L){                             //correction is needed
         int sgn=1; if(used > L) sgn=-1;     // direction of correction
         vector<pair<prec,avar>> v(m);      // heap of correction results of different symbols
         v.clear();                             // to prevent reallocation
         for(int i=0;i<m;i++){               // 
             if(q[i]+sgn) v.push_back({cc[i] - pow(pL[i]-(q[i]+sgn),2)*ip[i],i});
         }
         make_heap(v.begin(),v.end());  
         for( ;used!=L; used+=sgn){
             auto par = v.front(); pop_heap(v.begin(),v.end());v.pop_back();
             cc[par.second]= -par.first;
             if((q[par.second]+=sgn)+sgn){
                v.push_back({- par.first - pow(pL[par.second]-(q[par.second]+sgn),2)*ip[par.second], par.second}); 
                push_heap (v.begin(),v.end());    
             }
         }
     }
     delete[] cc; delete[] pL; delete[] ip;
 }
 
// -----------------  SYMBOL SPREADS ---------------------------- 
 void spread_fast(){                // the original Yann Collet's O(L) symbol spread
     tvar pos=0, step = (L >> 1) + (L >> 3) + 3;           // STARTING POSITION AND STEP
     delete[] s; s = new avar[L];
     tvar mask = L-1;                                   // assume L is a power of 2 !
     for(avar sym=0;sym<m;sym++){
        for(int i=0;i<q[sym];i++) {s[pos]=sym; pos = (pos+step) & mask;}    
     }       
 }
 void spread_prec(){     // O(L) precise spread - when we know only q
     delete[] s; s = new avar[L];
     tvar sym[L], first[L], next[L];
     int cp, pos=0; 
     for(int i=0; i<L; i++) first[i]=next[i]=L;       // empty all lists (L is the guardian)
     for(int sm=0; sm<m; sm++){                       // for each position, build lists of preferred symbols
         prec step = L/((prec)q[sm]), cur = step/2;          // STEP AND STARTING POSITION
         for(int i=0;i<q[sm];i++){
             sym[pos]=sm;
             tvar ins=round(cur); cur+=step;          // PREFERRED POSITION FOR THIS SYMBOL OCCURANCE
             next[pos]=first[ins]; first[ins]=pos++;         // insert in this position
         }
     }    
     pos=first[cp=0];
     for(int i=0;i<L;i++){                          // use succeeding lists for symbol spread
         while(pos==L) pos=first[++cp];             // is empty, look for another nonempty list
         s[i]=sym[pos];                             // choose the symbol
         pos=next[pos];                             // next symbol for this position
     }   
 }
 void spread_tuned(){        // O(L) tuned spread - uses both q and p
     delete[] s; s = new avar[L];
     tvar sym[L], first[L], next[L];         
     int cp, pos=0; 
     for(int i=0; i<=L; i++) first[i]=next[i]=L;       // empty all lists (L is the guardian)
     for(int sm=0; sm<m; sm++){                       // for each position, build lists of preferred symbols         
         for(int i=q[sm];i<2*q[sm];i++){
             sym[pos]=sm;
             int ins=round(1/(p[sm]*log(1+1/(prec)i))-L);    // PREFERRED POSITION FOR THIS SYMBOL OCCURANCE
             ins=min(max(ins,0),L-1);                        // low probable symbols would go above
             next[pos]=first[ins]; first[ins]=pos++;         // insert in this position
         }
     }    cout<<endl;
     pos=first[cp=0];
     for(int i=0;i<L;i++){                          // use succeeding lists for symbol spread
         while(pos==L) pos=first[++cp];             // is empty, look for another nonempty list
         s[i]=sym[pos];                             // choose the symbol
         pos=next[pos];                             // next symbol for this position
     }      
 }    
 // ---------- finding stationary probability and hANS  ----------------
  void heapify(){
      for(int x=2*L-2;x;x-=2) sp[x>>1]=sp[x]+sp[x+1];                     // build heap of sums
      if(abs(sp[1]-1)>acc) 
      {prec c=1/sp[1]; 
      for(int x=1;x<2*L;x++)sp[x]*=c;}    //normalize to 1
  }
 void make_step(){                       // propagate stationary probability by one step
     prec * nsp = new prec[2*L];         // new sp
     tvar cp[m];                         // current positions of symbols
     heapify();
     for(int i=0;i<m;i++) cp[i]=q[i];
     for(int i=0;i<L;i++) nsp[i+L] = p[s[i]] * sp[cp[s[i]]++];  // find new sp for position i        
     delete[] sp; sp=nsp;
 }
  void calc_hANS(){      
      hANS=0;
      prec nb[m];                         // number of bits used for symbols
      tvar cp[m];                         // current positions of symbols
      heapify();
      for(int i=0;i<m;i++){cp[i]=q[i]; nb[i]=0; int j=q[i]; while(j<L){j<<=1; nb[i]++;}};
      for(int i=0;i<L;i++){
          int sym=s[i];          
          hANS += nb[sym] * p[sym] * sp[cp[sym]++];
          if(!(cp[sym] &(cp[sym]-1))) nb[sym]--;        // if a power of 2
      }
  }
 void find_sp(){      // finds sp: stationary probability distribution of states 
    sp = new prec[2*L];
    for(int x=L;x < 2*L; x++) sp[x]=1/(1+(prec)x); 
    hANS=1; prec ht=0; 
    heapify();
    while(abs(ht-hANS)>acc){        
        for(int i=0;i<10;i++) make_step();            
        ht=hANS;calc_hANS();
    }
 }
};
// ----------- probability distribution initializers --------------
ANS init_binary(prec pr, int n){   // n binary distributions (m=2^n)
    ANS tmp; tmp.p = new prec[tmp.m = 1<<n];
    for(int i=0;i<tmp.m;i++){
        int j=i,s=0; while(j){s+=j&1;j>>=1;}; tmp.p[i]=pow(pr,s)*pow(1-pr,n-s);}
    sort(tmp.p,tmp.p+tmp.m);                                // not needed!
    return tmp;
}
ANS init_power(prec rho, int n){   //p_i propto rho^i distribution
    ANS tmp; tmp.p = new prec[tmp.m = n];
    for(int i=0; i<n; i++) tmp.p[i]=pow(rho,i);
    tmp.normalize();
    return tmp;
}
ANS init_rand_unif(int n){   // simple random distribution
    ANS tmp; tmp.p = new prec[tmp.m = n];
    srand(time(0)); for(int i=0; i<n; i++) tmp.p[i]=rand();
    tmp.normalize();
    sort(tmp.p,tmp.p+tmp.m);                                // not needed!
    return tmp;
}

int main()    // currently: single test
{// choose probability distribution
ANS test=init_binary(0.2,2);         // n binary  variables,  m=2^n             
// ANS test=init_power(0.95,256);    // p_i ~ rho^i,  m=256
//ANS test=init_rand_unif(256);      // m = 256
// test.printp();
test.quantize_prec(4);                        // choose quantization
// test.quantize_fast(12);                      // L = 2^value
// test.printq();
// int sum=0; for(int i=0;i<test.m;i++)sum+=test.q[i]; cout<<"sum:"<<sum<<endl;
test.calc_h();                              // find entropies
//test.spread_fast();                       // choose symbol spread
//test.spread_prec();
test.spread_tuned();
// test.prints();
test.find_sp();                             // find stationary probability and hANS

cout<<"Size m="<<test.m<<" alphabet with ENTROPY H="<<test.h<<endl;
if(test.m < 20) test.printp();
cout<<"HUFFMAN uses " << test.hhuff <<" ~ (1 + "<<(test.hhuff-test.h)/test.h<<
        ")H bits/symbol for depth "<< test.huff_depth<<" tree ("<<pow(2,test.huff_depth)<<" states)" <<endl; 
cout<<"We first QUANTIZE it to fractions with denominator L="<<test.L<<endl;
if(test.m < 20) test.printq();
cout<<"Entropy for QUANTIZATION grows to " <<test.hq<< " ~ (1 + "<<(test.hq-test.h)/test.h<<")H bits/symbol"<< endl;
cout<<"Then perform symbol spread"<<endl;
if(test.m < 11) test.prints();
cout<<"Entropy for its tANS is " <<test.hANS<< " ~ (1 + "<<(test.hANS-test.h)/test.h<<")H bits/symbol"<<endl;
return 0;
}
