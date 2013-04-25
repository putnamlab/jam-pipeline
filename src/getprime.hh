// get_prime, Copyright 2004, 2005 by Paul Havlak, Baylor College of Medicine
//
// Function to return largest prime <= target.  Good when you want one
//     large prime near a certain value, for example while constructing
//     a hash table and hash function.
// Tuned for large primes:
//
//   -- Uses Sieve of Eratosthenes to identify primes up to sqrt(N)
//   -- Uses modified sieve to identify primes between (N-10*log(N)) and N
//   -- Returns largest prime <= N, else zero if failed to find a prime
//      (I've found in the literature that the average density of primes is 
//      related to log(N), and that the is_hi_prime range used always contains 
//      many primes, but haven't found any guarantees.)
//
// Paul Havlak, 17 May 2004

#include <cmath>
#include <cstdlib>
#include <climits>
#include <iostream>

unsigned long get_prime(unsigned long target) {
  if (target < 2) return (target);

  double d_targ = (double) target;

  // Search for primes in low range:
  //      [2, lo_range]
  //
  unsigned long lo_range = (unsigned long) (sqrt(d_targ));

  // Protect against truncation/rounding errors (even using ceil() on
  // the sqrt() above wouldn't protect against rounding down in long
  // to double conversion).
  for (; lo_range < (target-1) && lo_range * lo_range < target; lo_range++) {
  }

  // Search for primes in high range: 
  //      [hi_start, target] (a range of hi_range+1 values)
  //
  unsigned long hi_range = (unsigned long) ceil(10*log(d_targ));
  unsigned long hi_start = (target - hi_range);
  if (hi_range > target) {
    hi_start = 0;
    hi_range = target;
  }  

  bool is_lo_prime[lo_range+1]; // whether lo_idx is prime
  bool is_hi_prime[hi_range+1]; // whether hi_idx + hi_start is prime
  long lo_idx; // signed to avoid underflow to maxint
  long hi_idx; // signed to avoid underflow to maxint

  for (lo_idx = 2; lo_idx <= lo_range; lo_idx++) {
    is_lo_prime[lo_idx] = true;
  }
  for (hi_idx = 0; hi_idx <= hi_range; hi_idx++) {
    is_hi_prime[hi_idx] = true;
  }

  unsigned long prime = 2;
  while (prime <= lo_range) {

    // Sieve through low range
    //
    for (lo_idx = 2 * prime; lo_idx <= lo_range; lo_idx += prime) {
      is_lo_prime[lo_idx] = false;
    }

    // Sieve through high range -- downwards is simpler than upwards.
    // Doing this with hi_start subtracted from the numbers we're
    // really interested in.
    //
    for (hi_idx = (prime * (target/prime) - hi_start);
	 hi_idx >= 0; 
	 hi_idx -= prime) {
      is_hi_prime[hi_idx] = false;
    }

    // Advance to next prime in low range -- those are the only primes
    // we need to worry about dividing anything in 2..target
    //
    while (++prime <= lo_range && !is_lo_prime[prime]) {}
  }

  // for (j = hi_start + 1; j <= target; j++) {
  //   if (is_hi_prime[j - hi_start]) {
  //     std::cout << j << " is prime" << std::endl;
  //   }
  // }

  // Return largest prime less than or equal to target (= hi_range + hi_start)
  for (hi_idx = hi_range; hi_idx >= 0; hi_idx--) {
    if (is_hi_prime[hi_idx]) {
      return hi_idx + hi_start;
    }
  }

  // We should only fail for target <= 1, no lesser primes.
  // But it's possible I haven't made hi_range large enough, didn't find
  // conclusive theory guaranteeing a prime in that particular range.
  return 0;
}
