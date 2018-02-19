#ifndef PBGL2_EQN_BUILDER
#define PBGL2_EQN_BUILDER
class linear_equation {
private:
  double m;
  double c;

public:
  linear_equation(double _m, double _c) : m(_m),
					  c(_c){}

  uint64_t operator()(uint64_t x) {
    return static_cast<uint64_t>(m*x + c); 
  }
};
#endif
