#include "flint/flint.h"
#include "flint/arb.h"
#include "flint/fmpz.h"
#include "flint/arb_hypgeom.h"
#include"../lib/fftw-3.3.10/api/fftw3.h"
#include<random>
#include<string>
#include<cstdlib>
#include<cmath>
#include"randgen.h"
#include<chrono>
#include <arb_hypgeom.h>
#include<iostream>
#include "../lib/matplotlibcpp.h"
#include<set>
namespace plt = matplotlibcpp;

double igamc(double a, double x);
size_t count_substrings(std::string const &needle, std::string const &haystack);
std::vector<std::vector<int>> getCombinations(int len, std::set<int> alphabet);
std::pair<double, double> serial_test(std::string bitstream, int m);
double fourier_test(std::string bitstream);
double monobit_test(std::string bitstream);
double block_test(std::string bitstream, int m);
double runs_test(std::string bitstream);
void randomness_test(RandGen* generator, int streamlen, int iterations);

int main()
{   
    /*
    RandGen *randor = new ORpsr(50);
    printf("OR Pseudorandom Number Generator Initialized\n");

    randomness_test(randor, 1000, 1);
    
    */

   
    RandGen *rand = new RSApsr(512);
    printf("RSA Pseudorandom Number Generator Initialized\n");

    randomness_test(rand, 1000, 1000);
    

    
}

void randomness_test(RandGen* generator, int streamlen, int iterations) {

    int failcount = 0;

    std::vector<double> pdata;
    

    for(int i = 0; i < iterations; i++) {
        std::string bitstream;
        for(int j = 0; j < streamlen; j++) {
            char bit = (*generator).randbit() + '0';

            bitstream.push_back(bit);
        }

        double p = monobit_test(bitstream);

        pdata.push_back(p);
        

        if(p < 0.01) {
            failcount++;
        }
    }

    printf("MONOBIT | Fails: %d, Failrate: %f\n", failcount, (double)failcount / iterations);

    failcount = 0;

    for(int i = 0; i < iterations; i++) {
        std::string bitstream;
        for(int j = 0; j < streamlen; j++) {
            char bit = (*generator).randbit() + '0';

            bitstream.push_back(bit);
        }

        double p = block_test(bitstream, floor(log2((double)streamlen)));

        if(p < 0.01) {
            failcount++;
        }
    }
        
    printf("BLOCK | Fails: %d, Failrate: %f\n", failcount, (double)failcount / iterations);
    
    failcount = 0;

    for(int i = 0; i < iterations; i++) {
        std::string bitstream;
        for(int j = 0; j < streamlen; j++) {
            char bit = (*generator).randbit() + '0';

            bitstream.push_back(bit);
        }

        double p = runs_test(bitstream);

        if(p < 0.01) {
            failcount++;
        }
    }
        
    printf("RUNS | Fails: %d, Failrate: %f\n", failcount, (double)failcount / iterations);

    failcount = 0;

    for(int i = 0; i < iterations; i++) {
        std::string bitstream;
        for(int j = 0; j < streamlen; j++) {
            char bit = (*generator).randbit() + '0';

            bitstream.push_back(bit);
        }   

        double p = fourier_test(bitstream);

        if(p < 0.01) {
            failcount++;
        }
    }
        
    printf("FOURIER | Fails: %d, Failrate: %f\n", failcount, (double)failcount / iterations);

    failcount = 0;

    for(int i = 0; i < iterations; i++) {
        std::string bitstream;
        for(int j = 0; j < streamlen; j++) {
            char bit = (*generator).randbit() + '0';

            bitstream.push_back(bit);
        }   

        std::pair<double, double> p = serial_test(bitstream, 3);

        if(p.first < 0.01 || p.second < 0.01) {
            failcount++;
        }
    }
        
    printf("SERIAL | Fails: %d, Failrate: %f\n", failcount, (double)failcount / iterations);


}

std::pair<double, double> serial_test(std::string bitstream, int m) {


    std::vector<std::string> augstreams;
    for(int j = 1; j <= 3; j++) {
        std::string aug = bitstream;
        for(int i = 0; i < m - j; i++) {
            aug.push_back(aug[i]);
        }
        augstreams.push_back(aug);
    }

    std::vector<std::vector<std::vector<int>>> vpattern;

    std::vector<std::vector<int>> combinations;
    std::vector<std::vector<int>> permutations;

    std::set<int> alphabet;
    alphabet.insert(0);
    alphabet.insert(1);

    combinations = getCombinations(m, alphabet);
    for(std::vector<int> combo : combinations) {
        do {
            permutations.push_back(combo);

        } while(next_permutation(combo.begin(), combo.end()));
    }
    vpattern.push_back(permutations);
    permutations.clear();

    combinations = getCombinations(m - 1, alphabet);
    for(std::vector<int> combo : combinations) {
        do {
            permutations.push_back(combo);
        } while(next_permutation(combo.begin(), combo.end()));
    }
    vpattern.push_back(permutations);
    permutations.clear();

    combinations = getCombinations(m - 2, alphabet);
    for(std::vector<int> combo : combinations) {
        do {
            permutations.push_back(combo);
        } while(next_permutation(combo.begin(), combo.end()));
    }
    vpattern.push_back(permutations);
    permutations.clear();


    std::vector<std::string> substrings;
    std::vector<std::vector<std::string>> bitpattern;


    for(unsigned long i = 0; i < vpattern.size(); i++) {
        for(std::vector<int> bits : vpattern[i]) {
            std::string bitstring;
            for(unsigned long i = 0; i < bits.size(); i++) {
                if(bits[i] == 1) {
                    bitstring.push_back('1');
                }
                if(bits[i] == 0) {
                    bitstring.push_back('0');
                }
            }
            substrings.push_back(bitstring);
        }

        bitpattern.push_back(substrings);
        substrings.clear();
    }

    std::vector<double> psivec;

    for(unsigned long j = 0; j < bitpattern.size(); j++) {

        std::vector<std::string> bitlength = bitpattern[j];

        double psi = 0;
        int sum = 0;

        for(unsigned long i = 0; i < bitlength.size(); i++) {
            size_t occurences = count_substrings(bitlength[i], augstreams[j]);

            sum = sum + ((occurences * occurences));
        }

        psi = (sum * (pow(2, m - j) / bitstream.length())) - bitstream.length();

        psivec.push_back(psi);
    }


    double psisquared = 0;
    double psisquaredtwo = 0;

    psisquared = psivec[0] - psivec[1];
    psisquaredtwo = psivec[0] - (2 * psivec[1]) + psivec[2];

    double p1 = igamc(pow(2, m - 2), psisquared / 2);
    double p2 = igamc(pow(2, m - 3), psisquaredtwo / 2);

    return std::make_pair(p1, p2);
}

double igamc(double a, double x) {
    arf_t s;
    arf_init(s);
    arf_set_d(s, a);
    arf_t z;
    arf_init(z);
    arf_set_d(z, x);

    arb_t sarb;
    arb_init(sarb);
    arb_t zarb;
    arb_init(zarb);

    arb_set_arf(sarb, s);
    arb_set_arf(zarb, z);

    arb_t res;
    arb_init(res);
    
    arb_hypgeom_gamma_upper(res, sarb, zarb, 1, 10000);

    arf_struct num = (*res).mid;

    return arf_get_d(&num, ARF_RND_DOWN);
}

size_t count_substrings(std::string const &needle, std::string const &haystack) {
    size_t count = 0;

    for (size_t pos =0; (pos=haystack.find(needle, pos)) != std::string::npos; ++pos, ++count)
        ;

   return count;
}

std::vector<std::vector<int>> getCombinations(int len, std::set<int> alphabet) {

    if(len < 1) {
        std::vector<std::vector<int>> empty;
        return empty;
    }

    std::vector<std::vector<int>> combinations;
    std::set<int> subalphabet = alphabet;

    if(len == 1) {
        for(int letter : alphabet) {
            std::vector<int> part;
            part.push_back(letter);
            combinations.push_back(part);
        }
        return combinations;
    }

    for(int letter : alphabet) {
        std::vector<std::vector<int>> part;

        std::vector<std::vector<int>> result = getCombinations(len - 1, subalphabet);
        part.insert(part.end(), result.begin(), result.end());

        for(std::vector<int>& combo : part) {
            combo.insert(combo.begin(), letter);
        }

        combinations.insert(combinations.end(), part.begin(), part.end());

        subalphabet.erase(letter);
    }

    return combinations;
} 

double fourier_test(std::string bitstream) {

    fftw_complex *in, *out;
    fftw_plan plan;

    in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * bitstream.length());
    out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * bitstream.length());
    plan = fftw_plan_dft_1d(bitstream.length(), in, out, FFTW_FORWARD, FFTW_ESTIMATE);

    for(unsigned long i = 0; i < bitstream.length(); i++) {
        if(bitstream[i] == '1') {
            in[i][0] = 1;
        } else {
            in[i][0] = -1;
        }
        in[i][1] = 0;
    }

    fftw_execute(plan);

    std::vector<double> s;

    for(unsigned long i = 0; i < bitstream.length() / 2; i++) {
        double mag = out[i][0] * out[i][0] + out[i][1] * out[i][1];

        mag = sqrt(mag);

        s.push_back(mag);
    }

    double t = log(1 / 0.05) * bitstream.length();
    t = sqrt(t);

    double n_zero = 0.95 * bitstream.length() / 2;

    int n_one = 0;

    for(unsigned long i = 0; i < s.size(); i++) {
        if(s[i] < t) {
            n_one++;
        }
    }

    double d = (n_one - n_zero) / sqrt((bitstream.length() * (0.95) * (0.05)) / 4);

    double p = erfc(abs(d) / sqrt(2));

    fftw_free(in);
    fftw_free(out);
    fftw_destroy_plan(plan);

    return p;
}

double block_runs_test(std::string bitstream) {

    int m = 4;

    if(bitstream.length() > 128) {
        m = 8;
    } else if(bitstream.length() > 6272) {
        m = 128;
    } else if(bitstream.length() > 750000) {
        m = 10000;
    }


    std::vector<std::string> streams;

    std::string part;
    for(unsigned long i = 0; i < bitstream.length(); i++) {
        part.push_back(bitstream[i]);

        if((int)part.length() == m) {
            streams.push_back(part);
            part.clear();
        }
    }

    return 0.5;


}

double runs_test(std::string bitstream) {
    int ones = 0;
    for(long unsigned i = 0; i < bitstream.length(); i++) {
        if(bitstream[i] == '1') {
            ones++;
        }
    } 

    double pi = (double)ones / bitstream.length();
    

    int Vobs = 1;

    for(unsigned long i = 0; i < bitstream.length() - 1; i++) {
        if(bitstream[i] == bitstream[i + 1]) {
            Vobs = Vobs;
        } else {
            Vobs = Vobs + 1;
        }
    }

    double numer = Vobs - (2 * bitstream.length() * pi * (1 - pi));
    numer = abs(numer);

    double denom = 2 * sqrt(2 * bitstream.length()) * pi * (1 - pi);

    double p = erfc(numer / denom);

    return p;

}

double monobit_test(std::string bitstream) {
    
    int ones = 0;

    //std::vector<double> rdata;

    for(long unsigned i = 0; i < bitstream.length(); i++) {
        if(bitstream[i] == '1') {
            ones++;
        }

        //double ratio = (double)ones / (i + 1);

        //rdata.push_back(ratio);
    } 

    //plt::xlabel("Bit stream length (bits)");
    //plt::ylabel("Percentage of 1's");

    //plt::plot(rdata);

    //plt::save("./basic.png");

    int sum = 2 * ones - bitstream.length();

    sum = abs(sum);
    
    double sumobs = sum / sqrt(bitstream.length());

    double p = erfc(sumobs / sqrt(2));

    return p;
}

double block_test(std::string bitstream, int m) {
    std::vector<std::string> streams;

    std::string part;
    for(unsigned long i = 0; i < bitstream.length(); i++) {
        part.push_back(bitstream[i]);

        if((int)part.length() == m) {
            streams.push_back(part);
            part.clear();
        }
    }

    double chisq = 0;

    for(std::string bstream : streams) {

        int ones = 0;
        for(char bit : bstream) {
            if(bit == '1') {
                ones++;
            }
        }

        double prop = ((double)ones / m);

        prop = (prop - 0.5) * (prop - 0.5);
        chisq += prop;
    }

    chisq = chisq * 4 * m;

    double p = igamc((double)streams.size() / 2, chisq / 2);

    return p;
}
