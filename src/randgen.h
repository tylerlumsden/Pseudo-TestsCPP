#include "flint/flint.h"
#include "flint/arb.h"
#include "flint/fmpz.h"
#include<random>

fmpz true_rand(int bits);
fmpz true_randprime(int bits);

class RandGen {
    public:
        virtual int randbit() = 0;
};

class ORpsr : public RandGen {
    private:
        fmpz seed;
        fmpz shiftseed;
        fmpz randqueue;
        int queue;
        int max;
    public:
        ORpsr(int maxbits) {
            
            max = maxbits;
            fmpz_init(&seed);
            fmpz_init(&shiftseed);
            fmpz_init(&randqueue);

            seed = true_rand(maxbits);

            shiftseed = seed;

            //shift bits of shiftseed twice to the left
            fmpz_mul_ui(&shiftseed, &shiftseed, 32);

            fmpz_or(&randqueue, &seed, &shiftseed);

            shiftseed = randqueue;

            queue = maxbits;
        }

        int randbit() override {
            if(queue == 0) {
                fmpz_mul_ui(&shiftseed, &shiftseed, 32);
                fmpz_or(&randqueue, &seed, &shiftseed);
                shiftseed = randqueue;
                queue = max;
            }

            //retrieve bit at lowest index
            int bit = fmpz_fdiv_ui(&randqueue, 2);
            //shift randqueue to the right by one
            fmpz_fdiv_q_ui(&randqueue, &randqueue, 2);

            queue--;

            return bit;
        }
};

class RSApsr : public RandGen {
    private:
        fmpz primex;
        fmpz primey;
        fmpz n;
        fmpz x;
        ulong e;
    public: 
        RSApsr(int bits) {

            fmpz_init(&primex);
            fmpz_init(&primey);
            fmpz_init(&n);
            fmpz_init(&x);

            primex = true_randprime(bits);
            primey = true_randprime(bits);
            fmpz_mul(&n, &primex, &primey);
            e = 3;
            
            do {
                x = true_rand(fmpz_bits(&n));
            } while(!fmpz_cmp(&x, &n));
        }

        int randbit() override {

            fmpz z;
            fmpz_init(&z);

            fmpz_pow_ui(&x, &x, e);
            fmpz_mod(&x, &x, &n);

            fmpz_mod_ui(&z, &x, 2);

            return z;
        }
};
