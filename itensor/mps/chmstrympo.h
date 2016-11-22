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
    operator<(const OneBodyInt& other) const { return (i < other.i || (i == other.i && j < other.j)); }

    bool
    proportialTo(const OneBodyInt& other) const;

    bool
    adjointOf(const OneBodyInt& other) const;

    };

struct TwoBodyInt
    {
    int i;
    int j;
    int k;
    int l;

    Cplx coef;

    TwoBodyInt(int i,
               int j,
               int k,
               int l,
               Cplx coef = 1);


    bool
    operator==(const TwoBodyInt& other) const;

    bool
    operator!=(const TwoBodyInt& other) const { return !operator==(other); }

    bool
    operator<(const TwoBodyInt& other) const { return (i < other.i ||
                                                      (i == other.i && j < other.j) ||
                                                      (i == other.i && j == other.j && k < other.k)  ||
                                                      (i == other.i && j == other.j && k == other.k && l < other.l));
                                                 }

    bool
    proportialTo(const TwoBodyInt& other) const;

    bool
    adjointOf(const TwoBodyInt& other) const;

    bool
    similarTo(const TwoBodyInt& other) const;

    };

class ChmstryMPO
    {
    SiteSet sites_;
    std::vector<OneBodyInt> oneBodyTerms_;
    std::vector<TwoBodyInt> twoBodyTerms_;

    public:

    ChmstryMPO(SiteSet const& sites)
      : sites_(sites)
        { }

    SiteSet const&
    sites() const { return sites_; }

    std::vector<OneBodyInt> const&
    oneBodyTerms() const { return oneBodyTerms_; }

    std::vector<TwoBodyInt> const&
    twoBodyTerms() const { return twoBodyTerms_; }

    operator MPO() const { return toMPO<ITensor>(*this); }

    operator IQMPO() const { return toMPO<IQTensor>(*this); }

    void
    add(OneBodyInt t)
        {
        if(abs(t.coef) > 0)
            {
            oneBodyTerms_.push_back(t);
            std::sort(oneBodyTerms_.begin(), oneBodyTerms_.end());
            }
        }

    void
    add(TwoBodyInt t)
        {
        if(abs(t.coef) > 0)
            {
            twoBodyTerms_.push_back(t);
            std::sort(twoBodyTerms_.begin(), twoBodyTerms_.end());
            }
        }

    void
    reset() { oneBodyTerms_.clear(); twoBodyTerms_.clear(); }

    };

std::ostream&
operator<<(std::ostream& s, const OneBodyInt& a);

std::ostream&
operator<<(std::ostream& s, const TwoBodyInt& a);

std::ostream&
operator<<(std::ostream& s, const ChmstryMPO& a);

}

#endif
