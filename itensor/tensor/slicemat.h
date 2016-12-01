//
// Distributed under the ITensor Library License, Version 1.2.
//    (See accompanying LICENSE file.)
//
#ifndef __ITENSOR_SLICEMAT_H_
#define __ITENSOR_SLICEMAT_H_

#include "itensor/tensor/mat.h"

namespace itensor {

template<typename Mat_>
auto
transpose(Mat_&& M)
    -> decltype(makeRef(std::forward<Mat_>(M),MatRange{}))
    {
    return makeRef(std::forward<Mat_>(M),transpose(M.range()));
    }

template<typename Mat_>
auto
subMatrix(Mat_&& M,
          size_t rstart,
          size_t rstop,
          size_t cstart,
          size_t cstop)
    -> decltype(makeRef(std::forward<Mat_>(M),MatRange{}))
    {
    static_assert(!std::is_same<Mat_&&,Matrix&&>::value,"Cannot pass temp/rvalue Matrix to subMatrix");
#ifdef DEBUG
    if(rstop > nrows(M) || rstart >= rstop) throw std::runtime_error("subMatrix invalid row start and stop");
    if(cstop > ncols(M) || cstart >= cstop) throw std::runtime_error("subMatrix invalid col start and stop");
#endif
    auto offset = rowStride(M)*rstart+colStride(M)*cstart;
    auto subrange = MatRange(rstop-rstart,rowStride(M),cstop-cstart,colStride(M));
    return makeRef(M.store()+offset,std::move(subrange));
    }

template<typename Mat_>
auto
rows(Mat_&& M,
     size_t rstart,
     size_t rstop)
    -> decltype(makeRef(std::forward<Mat_>(M),MatRange{}))
    {
    return subMatrix(std::forward<Mat_>(M),rstart,rstop,0,nrows(M));
    }

template<typename Mat_>
auto
columns(Mat_&& M,
        size_t cstart,
        size_t cstop)
    -> decltype(makeRef(std::forward<Mat_>(M),MatRange{}))
    {
    return subMatrix(std::forward<Mat_>(M),0,nrows(M),cstart,cstop);
    }

template<typename Mat_>
auto
diagonal(Mat_&& M)
    -> decltype(makeRef(std::forward<Mat_>(M).store(),VecRange{}))
    {
    static_assert(!std::is_same<Mat_&&,Matrix&&>::value,"Cannot pass temp/rvalue Matrix to diagonal(M)");
    auto drange = VecRange(std::min(nrows(M),ncols(M)),rowStride(M)+colStride(M));
    return makeRef(M.store(),std::move(drange));
    }

template<typename Mat_>
auto
row(Mat_&& M, size_t j)
    -> decltype(makeRef(std::forward<Mat_>(M).store(),VecRange{}))
    {
    static_assert(!std::is_same<Mat_&&,Matrix&&>::value,"Cannot pass temp/rvalue Matrix to row(M,n)");
#ifdef DEBUG
    if(j > nrows(M)) throw std::runtime_error("invalid row index");
#endif
    auto offset = j*rowStride(M);
    return makeRef(M.store()+offset,VecRange(ncols(M),colStride(M)));
    }

template<typename Mat_>
auto
column(Mat_&& M, size_t j)
    -> decltype(makeRef(std::forward<Mat_>(M).store(),VecRange{}))
    {
    static_assert(!std::is_same<Mat_&&,Matrix&&>::value,"Cannot pass temp/rvalue Matrix to column(M,n)");
#ifdef DEBUG
    if(j > ncols(M)) throw std::runtime_error("invalid column index");
#endif
    auto offset = j*colStride(M);
    return makeRef(M.store()+offset,VecRange(nrows(M),rowStride(M)));
    }


} //namespace itensor

#endif
