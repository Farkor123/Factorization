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
    std::vector<mpz_class> x;
    std::vector<mpz_class> ret;
    while (num % 2 == 0) {
      ret.push_back(mpz_class("2", 10));
      num /= 2;
    }
    x.push_back(mpz_class("2", 10));
    x.push_back(mpz_class("3", 10));
    x.push_back(mpz_class("8", 10));
    mpz_class q = 1, d;
    for (int i = 1; num != 1; i++) {
      bool flag = true;
      q *= (x[2*i] - x[i]) % num;
      d = utility::gcd(q, num);
      if (d > 1 && d < num) {
        if (is_prime(d)) {
          ret.push_back(d);
          num /= d;
          for (auto n : x) {
            n %= num;
          }
          i = 1;
          q = 1;
          flag = false;
        }
        else if (d < 65537) {
          ret.insert(ret.end(), lookup_table(d).begin(), lookup_table(d).end());
          num /= d;
          for (auto n : x) {
            n %= num;
          }
          i = 1;
          q = 1;
          flag = false;
        }
      }
      if (flag) {
        x.push_back(mpz_class((x.back() * x.back()) % num));
        x.push_back(mpz_class((x.back() * x.back()) % num));
      }
    }
    std::sort(x.begin(), x.end());
    return x;
  }
}
