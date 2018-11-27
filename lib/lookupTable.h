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
}
