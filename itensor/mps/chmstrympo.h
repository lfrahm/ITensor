//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_CHMSTRYMPO_H
#define __ITENSOR_CHMSTRYMPO_H

#include "itensor/global.h"
#include "itensor/mps/mpo.h"

namespace itensor {

class ChmstryMPO;

//
// Given an ChmstryMPO representing a Hamiltonian H,
// returns an exact (IQ)MPO form of H.
//
template <typename Tensor>
MPOt<Tensor>
toMPO(ChmstryMPO const& a,
      Args const& args = Args::global());


// //
// // Given an AutoMPO representing a Hamiltonian H,
// // returns an IQMPO which approximates exp(-tau*H)
// //
// // Although the tau argument is of Complex type, passing a Real
// // tau (Real is auto convertible to Complex) will 
// // result in a real-valued MPO.
// //
// // Arguments recognized:
// // o "Approx":
// //   - (Default) "ZW1" - Zaletel et al. "W1" approximation
// //
// template <typename Tensor>
// MPOt<Tensor>
// toExpH(const AutoMPO& a,
//        Complex tau,
//        const Args& args = Args::global());



//Instantiations of templates to allow us to define them
//later in autompo.cc
template<> MPO toMPO<ITensor>(const ChmstryMPO& a, const Args& args);
template<> IQMPO toMPO<IQTensor>(const ChmstryMPO& a, const Args& args);
// template<> MPO toExpH<ITensor>(const AutoMPO& a, Complex tau, const Args& args);
// template<> IQMPO toExpH<IQTensor>(const AutoMPO& a, Complex tau, const Args& args);
 
struct OneBodyInt
    {
    int i;
    int j;

    Cplx coef;

    OneBodyInt(int i,
               int j,
               Cplx coef = 1);


    bool
    operator==(const OneBodyInt& other) const;

    bool
    operator!=(const OneBodyInt& other) const { return !operator==(other); }

    bool 
    operator<(const OneBodyInt& other) const { return (i < other.i || j < other.j); }

    bool
    proportialTo(const OneBodyInt& other) const;

    bool
    adjointOf(const OneBodyInt& other) const;

    };

class ChmstryMPO
    {
    SiteSet sites_;
    std::vector<OneBodyInt> terms_;

    public:

    ChmstryMPO(SiteSet const& sites) 
      : sites_(sites)
        { }

    SiteSet const&
    sites() const { return sites_; }

    std::vector<OneBodyInt> const&
    terms() const { return terms_; }

    operator MPO() const { return toMPO<ITensor>(*this); }

    operator IQMPO() const { return toMPO<IQTensor>(*this); }

    void
    add(OneBodyInt t) 
        { 
        if(abs(t.coef) != 0) 
            {
            terms_.push_back(t); 
            std::sort(terms_.begin(), terms_.end());
            }
        }

    void
    reset() { terms_.clear(); }

    };

std::ostream& 
operator<<(std::ostream& s, const OneBodyInt& a);

std::ostream& 
operator<<(std::ostream& s, const ChmstryMPO& a);

}

#endif
