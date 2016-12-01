//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_LOCALMPOSET
#define __ITENSOR_LOCALMPOSET
#include "itensor/mps/localmpo.h"

namespace itensor {

template <class Tensor>
class LocalMPOSet
    {
    std::vector<MPOt<Tensor>> const* Op_ = nullptr;
    std::vector<LocalMPO<Tensor>> lmpo_;
    public:

    LocalMPOSet() { }

    LocalMPOSet(std::vector<MPOt<Tensor> > const& Op,
                Args const& args = Args::global());

    void
    product(Tensor const& phi, 
            Tensor & phip) const;

    Real
    expect(Tensor const& phi) const;

    Tensor
    deltaRho(Tensor const& AA, 
             Tensor const& comb, 
             Direction dir) const;

    Tensor
    diag() const;

    template <class MPSType>
    void
    position(int b, 
             MPSType const& psi);

    int
    numCenter() const { return lmpo_.front().numCenter(); }
    void
    numCenter(int val);

    int
    size() const { return lmpo_.front().size(); }

    explicit
    operator bool() const { return bool(Op_); }

    bool
    doWrite() const { return false; }
    void
    doWrite(bool val) 
        { 
        if(val) Error("Write to disk not yet supported LocalMPOSet");
        }

    };

template <class Tensor>
inline LocalMPOSet<Tensor>::
LocalMPOSet(std::vector<MPOt<Tensor>> const& Op,
            Args const& args)
  : Op_(&Op),
    lmpo_(Op.size())
    { 
    using LocalMPOT = LocalMPO<Tensor>;
    for(auto n : range(lmpo_.size()))
        {
        lmpo_[n] = LocalMPOT(Op.at(n));
        }
    }

template <class Tensor>
void inline LocalMPOSet<Tensor>::
product(Tensor const& phi, 
        Tensor & phip) const
    {
    lmpo_.front().product(phi,phip);

    Tensor phi_n;
    for(auto n : range(1,lmpo_.size()))
        {
        lmpo_[n].product(phi,phi_n);
        phip += phi_n;
        }
    }

template <class Tensor>
Real inline LocalMPOSet<Tensor>::
expect(Tensor const& phi) const
    {
    Real ex_ = 0;
    for(size_t n = 0; n < lmpo_.size(); ++n)
    for(auto n : range(lmpo_.size()))
        {
        ex_ += lmpo_[n].expect(phi);
        }
    return ex_;
    }

template <class Tensor>
Tensor inline LocalMPOSet<Tensor>::
deltaRho(Tensor const& AA,
         Tensor const& comb, 
         Direction dir) const
    {
    Tensor delta = lmpo_.front().deltaRho(AA,comb,dir);
    for(auto n : range(1,lmpo_.size()))
        {
        delta += lmpo_[n].deltaRho(AA,comb,dir);
        }
    return delta;
    }

template <class Tensor>
Tensor inline LocalMPOSet<Tensor>::
diag() const
    {
    Tensor D = lmpo_.front().diag();
    for(auto n : range(1,lmpo_.size()))
        {
        D += lmpo_[n].diag();
        }
    return D;
    }

template <class Tensor>
template <class MPSType> 
void inline LocalMPOSet<Tensor>::
position(int b, 
         MPSType const& psi)
    {
    for(auto n : range(lmpo_.size()))
        {
        lmpo_[n].position(b,psi);
        }
    }

template <class Tensor>
void inline LocalMPOSet<Tensor>::
numCenter(int val)
    {
    for(auto n : range(lmpo_.size()))
        {
        lmpo_[n].numCenter(val);
        }
    }

} //namespace itensor

#endif
