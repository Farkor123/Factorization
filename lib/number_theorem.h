#pragma once
#include "utility.h"

namespace number_theorem {
  int least_divisor(const int& num) {
    for(int i = 3; i <= sqrt(num); i += 2) {
      if (num % i == 0) {
        return i;
      }
    }
    return num;
  }

  std::vector<mpz_class> lookup_table (mpz_class& num) {
    mpz_class arr[65537] = {2};
    std::vector<mpz_class> vec;
    for (int i = 2; i < 65537; i++) {
      if (i % 2 == 0) {
        arr[i] = 2;
      }
      else {
        arr[i] = least_divisor(i);
      }
    }
    while (num != 1) {
      vec.push_back(arr[num.get_ui()]);
      num /= arr[num.get_ui()];
      }
    std::sort(vec.begin(), vec.end());
    return vec;
  }

  bool miller_rabin_test (mpz_class d, const mpz_class& n) {
    mpz_class x = 2 + rand() % (n-4);
    x = utility::pow_mod(x, d, n);
    if (x == 1 || x == n-1) {
      return true;
    }
    while (d != n-1) {
      x = (x * x) % n;
      d *= 2;
      if (x == 1) {
        return false;
      }
      if (x == n-1) {
        return true;
      }
    }
    return false;
  }

  bool is_prime(const mpz_class &num, int k = 256) {
    if (num <= 1 || num == 4) {
      return false;
    }
    if (num == 2 || num == 3) {
      return true;
    }
    mpz_class d = num - 1;
    while (d % 2 == 0) {
      d /= 2;
    }
    for (int i = 0; i < k; i++) {
      if(!miller_rabin_test(d, num)) {
        return false;
      }
    }
    return true;
  }

  std::vector<mpz_class> pollard_rho (mpz_class& num) {
    /*
     *  Actually, this implementation comes directly from J. M. Pollard's paper "AMonte Carlo method for factorization"
     */
    std::vector<mpz_class> x;
    std::vector<mpz_class> ret;
    //simply dividing by 2, ez one
    while (num % 2 == 0) {
      ret.push_back(2);
      num /= 2;
    }
    //if num was 2**n, no point to continue
    if (num == 1) {
      return ret;
    }
    //innitialise algorithm
    x.push_back(mpz_class(2));
    x.push_back(mpz_class(3));
    x.push_back(mpz_class(8));
    mpz_class q = 1, d;
    for (int i = 1; num != 1; i++) {
      //if at any point num is prime, we are finished
      if(is_prime(num)) {
        ret.push_back(num);
        std::sort(ret.begin(), ret.end());
        return ret;
      }
      //if num < 2**16 we can use lookup_table
      if(num < 65537) {
        auto w = lookup_table(num);
        ret.insert(ret.end(), w.begin(), w.end());
        std::sort(ret.begin(), ret.end());
        return ret;
      }
      q *= (x[2*i] - x[i]) % num;
      d = utility::gcd(abs(q), num);
      if (d > 1 && d < num) {
        //if d is prime, divide num by it
        if (is_prime(d)) {
          ret.push_back(d);
          num /= d;
          for (auto n : x) {
            n %= num;
          }
          i = 0;
          q = 1;
        }
        //if d isn't prime and d < 65537 we can factor it by lookup_table
        else if (d < 65537) {
          //std::cout << "DUPA";
          num /= d;
          auto w = lookup_table(d);
          ret.insert(ret.end(), w.begin(), w.end());
          for (auto n : x) {
            n %= num;
          }
          i = 0;
          q = 1;
        }
        //else we can factor it by pollard_rho
        else {
          num /= d;
          auto w = pollard_rho(d);
          ret.insert(ret.end(), w.begin(), w.end());
          for (auto n : x) {
            n %= num;
          }
          i = 0;
          q = 1;
        }
      }
      //if d is equal to num, we have to start with different x0
      else if (d == num) {
        srand(time(NULL));
        x.clear();
        x.push_back(mpz_class(rand() % (num-2) + 2));
        x.push_back(utility::x_squared_minus_one(x.back()) % num);
        x.push_back(utility::x_squared_minus_one(x.back()) % num);
        i = 0;
        q = 1;
      }
      //finally, if none of the above applies, we increment i by 1
      else {
        //any non-factorable polynomial applies, except for x**2-2
        x.push_back(utility::x_squared_minus_one(x.back()) % num);
        x.push_back(utility::x_squared_minus_one(x.back()) % num);
      }
    }
    std::sort(x.begin(), x.end());
    return x;
  }

  std::vector<mpz_class> simpler_pollard_rho (mpz_class& num) {
    std::vector<mpz_class> ret;
    gmp_randclass r(gmp_randinit_default);
    while (num % 2 == 0) {
      ret.push_back(2);
      num /= 2;
    }
    if (num == 1) {
      std::sort(ret.begin(), ret.end());
      return ret;
    }
    if (num < 65537) {
      auto w = lookup_table(num);
      ret.insert(ret.end(), w.begin(), w.end());
      std::sort(ret.begin(), ret.end());
      return ret;
    }
    if (is_prime(num)) {
      ret.push_back(num);
      std::sort(ret.begin(), ret.end());
      return ret;
    }
    mpz_class x = 2, y = 2, d = 1;
    FLAG:
    while (d == 1) {
      x = utility::x_squared_minus_one(x) % num;
      y = utility::x_squared_minus_one(utility::x_squared_minus_one(y)) % num;
      d = gcd(abs(x - y), num);
    }
    if (d == num) {
      x = y = r.get_z_range(num-3)+3;
      d = 1;
      goto FLAG;
    }
    else {
      if (is_prime(d)) {
        num /= d;
        ret.push_back(d);
        if(is_prime(num)) {
          ret.push_back(num);
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        if (num < 65537) {
          auto w = lookup_table(num);
          ret.insert(ret.end(), w.begin(), w.end());
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        x = y = 2;
        d = 1;
        goto FLAG;
      }
      else if(d < 65537) {
        num /= d;
        auto w = lookup_table(d);
        ret.insert(ret.end(), w.begin(), w.end());
        if (is_prime(num)) {
          ret.push_back(num);
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        if (num < 65537) {
          auto w = lookup_table(num);
          ret.insert(ret.end(), w.begin(), w.end());
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        x = y = 2;
        d = 1;
        goto FLAG;
      }
      else {
        num /= d;
        auto w = simpler_pollard_rho(d);
        ret.insert(ret.end(), w.begin(), w.end());
        if (is_prime(num)) {
          ret.push_back(num);
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        if (num < 65537) {
          auto w = lookup_table(num);
          ret.insert(ret.end(), w.begin(), w.end());
          std::sort(ret.begin(), ret.end());
          return ret;
        }
        x = y = 2;
        d = 1;
        goto FLAG;
      }
    }
  }

  std::vector<mpz_class> brent_pollard_rho (mpz_class& num) {
    /*
     * An implementation from Brent's paper "An improved Monte Carlo factorization algorithm."
     */

  }
}
