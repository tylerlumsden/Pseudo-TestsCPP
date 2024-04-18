#include "flint/flint.h"
#include "flint/arb.h"
#include "flint/fmpz.h"
#include "randgen.h"
#include<random>

fmpz true_rand(int bits) {

    fmpz_t num;
    fmpz_init(num);

    std::random_device ran;
    std::uniform_int_distribution<int> dist(0, 1);

    for(int i = 0; i < bits; i++) {
        ulong bit = dist(ran); 
        if(bit == 1) {
            fmpz_setbit(num, i);
        }
    }

    return *num;
}

fmpz true_randprime(int bits) {

    fmpz prime;

    do {
        prime = true_rand(bits);
    } while(!fmpz_is_probabprime(&prime) && !fmpz_is_prime(&prime));

    return prime;
}