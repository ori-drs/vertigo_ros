#include "AdaptiveLossFunction.h"

using namespace std;
namespace vertigo {

Adaptive::Adaptive(double c, double alpha, const ReweightScheme reweight) :
  Base(reweight), c_(c), alpha_(alpha), csquared_(c_*c_){
  if (c <= 0) {
    throw runtime_error("Adaptive mEstimator takes only positive double in constructor.");
  }
}

void Adaptive::print(const std::string &s="") const {
  cout << s << "apative (" << c_ << ")" << endl;
}

bool Adaptive::equals(const Base &expected, double tol) const {
  const Adaptive* p = dynamic_cast<const Adaptive*>(&expected);
  if (p == NULL) return false;
  return std::abs(c_ - p->c_) < tol;
}

Adaptive::shared_ptr Adaptive::Create(double c, double alpha, const ReweightScheme reweight) {
  return shared_ptr(new Adaptive(c, reweight));
}

} // vertigo namespace
