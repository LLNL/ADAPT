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

#include <cassert>
#include <cstddef>
#include <stack>
#include "MergeTree.h"

extern FunctionType* gData;

extern uint32_t gDim[3];


std::vector<Node>* Node::sNodes = NULL;

Node::Node(GlobalIndexType id, LocalIndexType i) : mIndex(id), mDown(LNULL), mUp(LNULL), mNext(i), mRep(LNULL)
{
}

Node::Node(const Node& n) : mIndex(n.mIndex), mDown(n.mDown), mUp(n.mUp), mNext(n.mNext), mRep(n.mRep)
{
}


FunctionType Node::arcLength() const
{
  if (mDown == LNULL)
    return 0;

  return fabs(gData[mIndex] - gData[(*sNodes)[mDown].mIndex]);
}

MergeTree::MergeTree()
{
  // For now there can be only a asingle merge tree
  assert (Node::sNodes == NULL);

  Node::sNodes = &mNodes;
}


LocalIndexType MergeTree::addCriticalPoint(GlobalIndexType id)
{
  LocalIndexType i = (LocalIndexType)mNodes.size();

  mNodes.push_back(Node(id,i));
  mArcs.push_back(Arc(id));

  return i;
}

int MergeTree::addEdge(LocalIndexType up, LocalIndexType down)
{
  assert(up < mNodes.size());
  assert(down < mNodes.size());

  assert (mNodes[up].down() == LNULL);

  mNodes[up].down(down);

  if (mNodes[down].up() == LNULL) {
    mNodes[down].up(up);
  }
  else {
    mNodes[up].next(mNodes[mNodes[down].up()].next());
    mNodes[mNodes[down].up()].next(up);
  }

  return 1;
}

int MergeTree::removeEdge(LocalIndexType up, LocalIndexType down)
{
  // First we unlink our sibling if there are any
  if (mNodes[up].next() != up) {

    // Make sure that the up-pointer of down will remain valid
    mNodes[down].up(mNodes[up].next());

    LocalIndexType prev = up;
    while (mNodes[prev].next() != up)
      prev = mNodes[prev].next();

    mNodes[prev].next(mNodes[up].next());
    mNodes[up].next(up);
  }
  else {
    mNodes[down].up(LNULL);
  }

  // Now unlink
  mNodes[up].down(LNULL);

  return 1;
}

int MergeTree::addVertex(GlobalIndexType v, LocalIndexType label)
{
  assert(label < mArcs.size());

  mArcs[label].mVertices.push_back(v);

  return 1;
}

int MergeTree::splitBySize(LocalIndexType n)
{
  assert (n > 0);

  LocalIndexType i = 0;

  while (i < mArcs.size()) {
    if (mArcs[i].size() > n)
      splitArc(i,(LocalIndexType)mArcs[i].mVertices.size()/2);
    else
      i++;
  }

  return 1;
}

int MergeTree::splitByLength(FunctionType l)
{
  assert (l > 0);

  LocalIndexType i = 0;

  while (i < mArcs.size()) {
    if ((mArcs[i].size() > 1) && (mNodes[i].arcLength() > l)) {
      LocalIndexType k=1;

      for (k=1;k<mArcs.size();k++) {
        if (fabs(gData[mArcs[i].mVertices[0]] - gData[mArcs[i].mVertices[k]]) > mNodes[i].arcLength()/2)
          break;

        k++;
      }

      splitArc(i,k);
    }
    else
      i++;

  }

  return 1;
}

int MergeTree::splitArc(LocalIndexType a, LocalIndexType pos)
{
  // First we create a new node
  LocalIndexType label = addCriticalPoint(mArcs[a].mVertices[pos]);

  // Copy all the vertices
  mArcs[label].mVertices.insert(mArcs[label].mVertices.begin(),
                                mArcs[a].mVertices.begin()+pos+1,
                                mArcs[a].mVertices.end());

  // Remove the vertices from the original arc
  mArcs[a].mVertices.erase(mArcs[a].mVertices.begin()+pos,mArcs[a].mVertices.end());

  LocalIndexType down = mNodes[a].down();

  // Now add the new vertex into the tree
  if (down != LNULL)
    removeEdge(a,down);

  addEdge(a,label);

  if (down != LNULL)
    addEdge(label,down);

  return 1;
}

void MergeTree::constructFeature(LocalIndexType label, std::vector<GlobalIndexType>& feature) const
{
  feature.insert(feature.end(),mArcs[label].mVertices.begin(),mArcs[label].mVertices.end());

  if (mNodes[label].up() != LNULL) {
    LocalIndexType up = mNodes[label].up();
    do {
      constructFeature(up,feature);

      up = mNodes[up].next();
    } while (up != mNodes[label].up());

  }
}

int MergeTree::inflate()
{
  std::stack<LocalIndexType> front;
  LocalIndexType top,next;

  for (LocalIndexType i=0;i<size();i++) {
    // For all roots
    if (node(i).down() == LNULL) {
      front.push(i);

      while (!front.empty()) {
        top = front.top();
        front.pop();

        if (node(top).up() != LNULL) {
          next = node(top).up();
          do {
            front.push(next);

            if (node(next).metric() < node(top).metric())
              node(next).metric(node(top).metric());

            next = node(next).next();
          } while (next != node(top).up());
        }
      }
    }
  }

  return 1;
}

int MergeTree::deflate()
{


  return 1;
}

