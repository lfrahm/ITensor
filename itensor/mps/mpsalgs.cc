//
// Distributed under the ITensor Library License, Version 1.2
//    (See accompanying LICENSE file.)
//
#include <map>
#include "itensor/mps/mps.h"
#include "itensor/mps/mpo.h"
#include "itensor/mps/localop.h"
#include "itensor/util/print_macro.h"

namespace itensor {

using std::map;
using std::istream;
using std::ostream;
using std::cout;
using std::endl;
using std::vector;
using std::find;
using std::pair;
using std::make_pair;
using std::string;

void
plussers(Index const& l1,
         Index const& l2,
         Index      & sumind,
         ITensor    & first,
         ITensor    & second)
    {
    auto m = l1.m()+l2.m();
    if(m <= 0) m = 1;
    sumind = Index(sumind.rawname(),m);

    first = delta(l1,sumind);
    auto S = Matrix(l2.m(),sumind.m());
    for(auto i : range(l2.m()))
        {
        S(i,l1.m()+i) = 1;
        }
    second = matrixTensor(std::move(S),l2,sumind);
    }

// IQTensor
// constructWithBlock(IQIndexSet all, Index a, Index b)
// {
//   IQTensor result(all);
//   assert(all.r() == 2);
//   auto A = all.front();
//   auto B = all.back();
//
//   auto offsetA = offset(A, a);
//   auto offsetB = offset(B, b);
//
//   for (int i = 1; i<=a.m() || i<=b.m(); i++) {
//     result.set(A(offsetA + i), B(offsetB + i), 1.0);
//   }
//
//   return result;
// }
//
// void
// plussers(IQIndex const& l1, IQIndex const& l2,
//          IQIndex& sumind,
//          IQTensor& first, IQTensor& second)
//     {
//     map<Index,Index> l1map, l2map;
//     vector<IndexQN> iq;
//     for(IndexQN const& x : l1)
//         {
//         Index jj(x.index.rawname(),x.m(),x.type());
//         l1map[x.index] = jj;
//         iq.push_back(IndexQN(jj,x.qn));
//         }
//     for(IndexQN const& x : l2)
//         {
//         Index jj(x.index.rawname(),x.m(),x.type());
//         l2map[x.index] = jj;
//         iq.push_back(IndexQN(jj,x.qn));
//         }
//     sumind = IQIndex(sumind.rawname(),std::move(iq),sumind.dir(),sumind.primeLevel());
//     first = IQTensor(dag(l1),sumind);
//     for(IndexQN const& il1 : l1)
//         {
//         Index& s1 = l1map[il1.index];
//         auto t = delta(il1.index,s1);
//         first += t;
//         }
//     second = IQTensor(dag(l2),sumind);
//     for(IndexQN const& il2 : l2)
//         {
//         Index& s2 = l2map[il2.index];
//         auto t = delta(il2.index,s2);
//         second += t;
//         }
//     }

// void
// plussers(IQIndex const& l1, IQIndex const& l2,
//          IQIndex& sumind,
//          IQTensor& first, IQTensor& second)
//     {
//     PrintData(l1);
//     PrintData(l2);
//     vector<int> l2ints;
//     vector<IndexQN> iq;
//     for (IndexQN const& x : l1)
//         {
//         Index jj(x.index.rawname(),x.m(),x.type());
//         for(IndexQN const& y: l2)
//             {
//             if(x.qn == y.qn)
//                 {
//                 for (int i = jj.m() + 1; i <= jj.m() + y.m(); i++)
//                     {
//                     l2ints.push_back(i);
//                     }
//                 jj = Index(jj.rawname() + y.index.rawname(), jj.m() + y.m(), jj.type());
//                 }
//             }
//         iq.push_back(IndexQN(jj, x.qn));
//         }
//     for (IndexQN const& y : l2)
//         {
//         bool add = true;
//         for(IndexQN const& x: l1) { if(x.qn == y.qn) { add=false; } }
//         if (add)
//             {
//             iq.push_back(IndexQN(Index(y.index.rawname(),y.m(),y.type()), y.qn));
//             for (int i = l2ints.back(); i <= l2ints.back() + y.m(); i++)
//                 {
//                 l2ints.push_back(i);
//                 }
//             }
//         }
//     sumind = IQIndex(sumind.rawname(),std::move(iq),sumind.dir(),sumind.primeLevel());
//     PrintData(sumind);
//     // PrintData(sumind);
//     first = IQTensor(dag(l1), sumind);
//     for (int i = 1; i <= l1.m(); i++) {
//         if (std::find(l2ints.begin(), l2ints.end(), i) == l2ints.end())
//             first.set(l1(i), sumind(i), 1.0);
//     }
//
//     second = IQTensor(dag(l2), sumind);
//     for (int i = 1; i <= l2.m(); i++) {
//         for(int j : l2ints)
//             {
//             second.set(l2(i), sumind(j), 1.0);
//             }
//     }
//
//     PrintData(first);
//     PrintData(second);
//     }

void
reArrow(MPO& a, MPO const& b)
    {

    }

void
reArrow(IQMPS& a, IQMPS const& b)
    {

    }

void
reArrow(MPS& a, MPS const& b)
    {

    }

void
reArrow(IQMPO& a,
        IQMPO const& b)
// void
// reArrow(MPOt<IQTensor>& a, MPOt<IQTensor> const& b)
    {
    // int N = a.N();
    // assert(b.N() == N);

    // auto ra = rightLinkInd(a, 1);
    // auto rb = rightLinkInd(b, 1);
    //
    // if (ra.dir() != rb.dir())
    //     {
    //     a.Anc(1).dag();
    //     }

    // for (int i = 2; i < N; i++)
    //     {
    //     auto ra = rightLinkInd(a, i);
    //     auto la = leftLinkInd(a, i);
    //     auto rb = rightLinkInd(b, i);
    //     auto lb = leftLinkInd(b, i);
    //     if (la.dir() != lb.dir())
    //         {
    //         PrintData(a.A(i));
    //         PrintData(b.A(i));
    //         std::vector<IndexQN> iq;
    //         for (IndexQN const& x: la) { iq.push_back(x); }
    //         auto laN = IQIndex(la.rawname(), std::move(iq), lb.dir(), 0);
    //         PrintData(delta(dag(la), laN));
    //         PrintData(a.A(i) * delta(dag(la), laN));
    //         a.Anc(i) *= delta(dag(la), laN);
    //         PrintData(a.A(i));
    //         }
    //     if (ra.dir() != rb.dir())
    //         {
    //         PrintData(a.A(i));
    //         PrintData(b.A(i));
    //         PrintData(delta(dag(ra)));
    //         // a.Anc(i) *= delta(dag(ra), ra);
    //         PrintData(a.A(i));
    //         }
    //     }
    }

void
plussers(IQIndex const& l1, IQIndex const& l2,
         IQIndex& sumind,
         IQTensor& first, IQTensor& second)
    {

    vector<IndexQN> iq;
    vector<IQIndex> map;

    for(IndexQN const& x : l1)
        {
        if(hasQN(l2, x.qn))
            {
            Index y = findByQN(l2, x.qn);
            Index jj(x.index.rawname() + y.rawname(),x.m() + y.m(),x.type());
            iq.push_back(IndexQN(jj, x.qn));
            for (int i = 1; i <= x.m(); i++) { map.push_back(l1); }
            for (int i = 1; i <= y.m(); i++) { map.push_back(l2); }
            }
        else
            {
            Index jj(x.index.rawname(),x.m(),x.type());
            iq.push_back(IndexQN(jj, x.qn));
            for (int i = 1; i <= x.m(); i++) { map.push_back(l1); }
            }
        }
    for(IndexQN const& y : l2)
        {
        if(!hasQN(l1, y.qn))
            {
            Index jj(y.index.rawname(),y.m(),y.type());
            iq.push_back(IndexQN(jj, y.qn));
            for (int i = 1; i <= y.m(); i++) { map.push_back(l2); }
            }
        }

    sumind = IQIndex(sumind.rawname(),std::move(iq),sumind.dir(),sumind.primeLevel());

    assert(sumind.m() == l1.m() + l2.m());
    first = IQTensor(dag(l1), sumind);
    int count = 1;
    for (int i = 1; i <= sumind.m(); i++)
        {
        if (map[i-1] == l1)
            {
            first.set(l1(count), sumind(i), 1.0);
            count++;
            }
        }

    second = IQTensor(dag(l2), sumind);
    count = 1;
    for (int i = 1; i <= sumind.m(); i++)
        {
        if (map[i-1] == l2)
            {
            second.set(l2(count), sumind(i), 1.0);
            count++;
            }
        }
    if (l1.dir() != l2.dir())
        {
        PrintData(first);
        PrintData(second);
        }
    }


//
// Adds two MPSs but doesn't attempt to
// orthogonalize them first
//
template <class MPSType>
MPSType&
addAssumeOrth(MPSType      & L,
              MPSType const& R,
              Args const& args)
    {
    using Tensor = typename MPSType::TensorT;
    reArrow(L, R);
    // PrintData(L);
    // PrintData(R);

    auto N = L.N();
    if(R.N() != N) Error("Mismatched MPS sizes");

    L.primelinks(0,4);

    auto first = vector<Tensor>(N);
    auto second = vector<Tensor>(N);

    for(auto i : range1(N-1))
        {
        auto l1 = rightLinkInd(L,i);
        auto l2 = rightLinkInd(R,i);
        auto r = l1;
        plussers(l1,l2,r,first[i],second[i]);
        }

    // PrintData(L.A(1) * first[1]);
    // PrintData(R.A(1) * second[1]);

    L.Anc(1) = L.A(1) * first[1] + R.A(1) * second[1];
    for(auto i : range1(2,N-1))
        {
        L.Anc(i) = dag(first[i-1]) * L.A(i) * first[i]
                     + dag(second[i-1]) * R.A(i) * second[i];
        }
    L.Anc(N) = dag(first[N-1]) * L.A(N) + dag(second[N-1]) * R.A(N);

    L.noprimelink();

    L.orthogonalize(args);

    return L;
    }
template MPS& addAssumeOrth(MPS & L,MPS const& R, Args const& args);
template IQMPS& addAssumeOrth(IQMPS & L,IQMPS const& R, Args const& args);
template MPO& addAssumeOrth(MPO & L,MPO const& R, Args const& args);
template IQMPO& addAssumeOrth(IQMPO & L,IQMPO const& R, Args const& args);

template <class Tensor>
void
fitWF(const MPSt<Tensor>& psi_basis, MPSt<Tensor>& psi_to_fit)
    {
    if(!itensor::isOrtho(psi_basis))
        Error("psi_basis must be orthogonolized.");
    if(orthoCenter(psi_basis) != 1)
        Error("psi_basis must be orthogonolized to site 1.");

    auto N = psi_basis.N();
    if(psi_to_fit.N() != N)
        Error("Wavefunctions must have same number of sites.");

    auto A = psi_to_fit.A(N) * dag(prime(psi_basis.A(N),Link));
    for(int n = N-1; n > 1; --n)
        {
        A *= dag(prime(psi_basis.A(n),Link));
        A *= psi_to_fit.A(n);
        }
    A = psi_to_fit.A(1) * A;
    A.noprime();

    auto nrm = norm(A);
    if(nrm == 0) Error("Zero overlap of psi_to_fit and psi_basis");
    A /= nrm;

    psi_to_fit = psi_basis;
    psi_to_fit.Anc(1) = A;
    }
template void fitWF(const MPSt<ITensor>& psi_basis, MPSt<ITensor>& psi_to_fit);
template void fitWF(const MPSt<IQTensor>& psi_basis, MPSt<IQTensor>& psi_to_fit);

bool
checkQNs(const IQMPS& psi)
    {
    const int N = psi.N();

    QN Zero;

    int center = findCenter(psi);
    if(center == -1)
        {
        cout << "Did not find an ortho. center\n";
        return false;
        }

    //Check that all IQTensors have zero div
    //except possibly the ortho. center
    for(int i = 1; i <= N; ++i)
        {
        if(i == center) continue;
        if(!psi.A(i))
            {
            println("A(",i,") null, QNs not well defined");
            return false;
            }
        if(div(psi.A(i)) != Zero)
            {
            cout << "At i = " << i << "\n";
            Print(psi.A(i));
            cout << "IQTensor other than the ortho center had non-zero divergence\n";
            return false;
            }
        }

    //Check arrows from left edge
    for(int i = 1; i < center; ++i)
        {
        if(rightLinkInd(psi,i).dir() != In)
            {
            println("checkQNs: At site ",i," to the left of the OC, Right side Link not pointing In");
            return false;
            }
        if(i > 1)
            {
            if(leftLinkInd(psi,i).dir() != Out)
                {
                println("checkQNs: At site ",i," to the left of the OC, Left side Link not pointing Out");
                return false;
                }
            }
        }

    //Check arrows from right edge
    for(int i = N; i > center; --i)
        {
        if(i < N)
        if(rightLinkInd(psi,i).dir() != Out)
            {
            println("checkQNs: At site ",i," to the right of the OC, Right side Link not pointing Out");
            return false;
            }
        if(leftLinkInd(psi,i).dir() != In)
            {
            println("checkQNs: At site ",i," to the right of the OC, Left side Link not pointing In");
            return false;
            }
        }

    //Done checking arrows
    return true;
    }

QN
totalQN(const IQMPS& psi)
    {
    const int center = findCenter(psi);
    if(center == -1)
        Error("Could not find ortho. center");
    return div(psi.A(center));
    }

} //namespace itensor
