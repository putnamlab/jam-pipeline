#include <string>

using namespace std;

class Debugging {
private:
  unsigned long bitvector;
public:
  Debugging() : bitvector(0) { }
  Debugging(const string flags) { 
    record(flags); 
  }
  void record(const string flags) {
    int i;
    for (i = 0; i < flags.length(); i++) {
      if (flags[i] == '+') {
	bitvector = ~0;
      }
      else if (flags[i] >= '@') {
	bitvector |= (1 << (flags[i] - 64));
      }
    }
  }
  inline bool check(char flag) {
    if (flag >= '@')
      return !!(bitvector & (1 << (flag - 64)));
    else
      return true;
  }
};
