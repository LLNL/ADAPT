/*******************************************************************************
* Copyright (c) 2015, Lawrence Livermore National Security, LLC
* Produced at the Lawrence Livermore National Laboratory
* Written by Peer-Timo Bremer bremer5@llnl.gov
* LLNL-CODE-665196
* All rights reserved.
* 
* This file is part of ADAPT. For details, see
* https://github.com/scalability-llnl/ADAPT. Please also read the
* additional BSD notice below. Redistribution and use in source and
* binary forms, with or without modification, are permitted provided
* that the following conditions are met:
* 
* - Redistributions of source code must retain the above copyright
*    notice, this list of conditions and the disclaimer below.
* 
* - Redistributions in binary form must reproduce the above copyright
*    notice, this list of conditions and the disclaimer (as noted below) in
*    the documentation and/or other materials provided with the
*    distribution.
* 
* - Neither the name of the LLNS/LLNL nor the names of its contributors
*    may be used to endorse or promote products derived from this software
*    without specific prior written permission.
* 
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
* A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE
* LIVERMORE NATIONAL SECURITY, LLC, THE U.S. DEPARTMENT OF ENERGY OR
* CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
* EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
* PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
* PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
* LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
* NEGLIGENCE OR OTHERWISE) ARISING ￼ IN ANY WAY OUT OF THE USE OF THIS
* SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
* 
* Additional BSD Notice 
* 
* 1. This notice is required to be provided under our contract with the
* U.S. Department of Energy (DOE). This work was produced at Lawrence
* Livermore National Laboratory under Contract No. DE-AC52-07NA27344
* with the DOE. 
* 
* 2. Neither the United States Government nor Lawrence Livermore
* National Security, LLC nor any of their employees, makes any warranty,
* express or implied, or assumes any liability or responsibility for the
* accuracy, completeness, or usefulness of any information, apparatus,
* product, or process disclosed, or represents that its use would not
* infringe privately-owned rights. 
* 
* 3. Also, reference herein to any specific commercial products,
* process, or services by trade name, trademark, manufacturer or
* otherwise does not necessarily constitute or imply its endorsement,
* recommendation, or favoring by the United States Government or
* Lawrence Livermore National Security, LLC. The views and opinions of
* authors expressed herein do not necessarily state or reflect those of
* the United States Government or Lawrence Livermore National Security,
* LLC, and shall not be used for advertising or product endorsement
* purposes.
********************************************************************************/

#include "Neighborhood.h"


Neighborhood::iterator::iterator()
: mNeighbors(NULL), mOffsets(NULL), mOrigin(0),mCount(0),mCurrent(0)
{
  mDim[0] = 0;
  mDim[1] = 0;
  mDim[2] = 0;

  mCoords[0] = mCoords[1] = mCoords[2] = 0;
}


Neighborhood::iterator::iterator(GlobalIndexType v, GlobalIndexType dim[3],int8_t* neighbors,SignedGlobalIndexType* offsets, uint8_t count)
: mNeighbors(neighbors), mOffsets(offsets), mOrigin(v),mCount(count),mCurrent(0)
{
  mDim[0] = dim[0];
  mDim[1] = dim[1];
  mDim[2] = dim[2];

  mCoords[0] = v % mDim[0];
  mCoords[1] = (v / mDim[0]) % mDim[1];
  mCoords[2] = v / (mDim[0]*mDim[1]);

  while ((mCurrent < mCount) && (!inside(mCurrent)))
    mCurrent++;
}

Neighborhood::iterator::iterator(const iterator& it) : mNeighbors(it.mNeighbors),mOffsets(it.mOffsets),mOrigin(it.mOrigin),
mCount(it.mCount), mCurrent(it.mCurrent)
{
  mDim[0] = it.mDim[0];
  mDim[1] = it.mDim[1];
  mDim[2] = it.mDim[2];

  mCoords[0] = it.mCoords[0];
  mCoords[1] = it.mCoords[1];
  mCoords[2] = it.mCoords[2];
}

Neighborhood::iterator& Neighborhood::iterator::operator=(const iterator& it)
{
  mNeighbors = it.mNeighbors;
  mOffsets = it.mOffsets;
  mOrigin = it.mOrigin;
  mCount = it.mCount;
  mCurrent = it.mCurrent;

  mDim[0] = it.mDim[0];
  mDim[1] = it.mDim[1];
  mDim[2] = it.mDim[2];

  mCoords[0] = it.mCoords[0];
  mCoords[1] = it.mCoords[1];
  mCoords[2] = it.mCoords[2];

  return *this;
}


bool Neighborhood::iterator::operator!=(const iterator& it)
{
  if ((mCurrent != it.mCurrent) || (mOrigin != it.mOrigin))
    return true;

  return false;
}

Neighborhood::iterator& Neighborhood::iterator::operator++(int i)
{
  mCurrent++;
  while (!inside(mCurrent) && (mCurrent < mCount))
    mCurrent++;

  return *this;
}

bool Neighborhood::iterator::inside(uint8_t i)
{
  if (   (((int32_t)mCoords[0] + mNeighbors[3*i + 0] < 0) || ((int32_t)mCoords[0] + mNeighbors[3*i + 0] >= mDim[0]))
      || (((int32_t)mCoords[1] + mNeighbors[3*i + 1] < 0) || ((int32_t)mCoords[1] + mNeighbors[3*i + 1] >= mDim[1]))
      || (((int32_t)mCoords[2] + mNeighbors[3*i + 2] < 0) || ((int32_t)mCoords[2] + mNeighbors[3*i + 2] >= mDim[2]))) {


      return false;
  }

  return true;
}



Neighborhood::Neighborhood(GlobalIndexType dim[3]) : mNeighbors(NULL), mOffsets(NULL), mCount(0)
{
  mDim[0] = dim[0];
  mDim[1] = dim[1];
  mDim[2] = dim[2];
}


Neighborhood::iterator Neighborhood::begin(GlobalIndexType origin)
{
  iterator it(origin,mDim,mNeighbors,mOffsets,mCount);

  return it;
}

Neighborhood::iterator Neighborhood::end(GlobalIndexType origin)
{
  iterator it(origin,mDim,mNeighbors,mOffsets,mCount);

  it.mCurrent = mCount;
  return it;
}

void Neighborhood::computeOffsets()
{
  for (uint8_t i=0;i<mCount;i++)
    mOffsets[i] = mNeighbors[3*i+2]*mDim[0]*mDim[1] + mNeighbors[3*i+1]*mDim[0] + mNeighbors[3*i];
}


