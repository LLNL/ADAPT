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

#ifndef MERGETREE_H_
#define MERGETREE_H_

#include <vector>
#include <cmath>

#include "Definitions.h"

//! The node of a merge tree
class Node
{
public:

  //! A reference ot the vector of nodes
  static std::vector<Node>* sNodes;

  //! Default constructor
  Node(GlobalIndexType id, LocalIndexType i);

  //! Copy constructor
  Node(const Node& n);

  //! Return the index
  GlobalIndexType index() const {return mIndex;}

  //! Return the descendant
  LocalIndexType down() const {return mDown;}

  //! Return one parent
  LocalIndexType up() const {return mUp;}

  //! Return a sibling
  LocalIndexType next() const {return mNext;}

  //! Return the current representative
  LocalIndexType rep() const {return mRep;}

  //! Return the metric value
  FunctionType metric() const {return mMetric;}

  //! Set the descendant
  void down(LocalIndexType d) {mDown = d;}

  //! Set the parent
  void up(LocalIndexType u) {mUp = u;}

  //! Set the sibling
  void next(LocalIndexType n) {mNext = n;}

  //! Set the representative
  void rep(LocalIndexType r) {mRep = r;}

  //! Set the metric value
  void metric(FunctionType m) {mMetric = m;}

  //! Return the current length of the arc
  FunctionType arcLength() const;

private:

  //! The global index of vertex corresponding to this node
  GlobalIndexType mIndex;

  //! The local index of the descendant
  LocalIndexType mDown;

  //! The local index of one parent
  LocalIndexType mUp;

  //! The local index of the next sibling
  LocalIndexType mNext;

  //! The local index of the highest node in the subtree
  LocalIndexType mRep;

  //! The metric of the node
  FunctionType mMetric;
};

//! The arc of a merge tree storing all corresponding vertices
class Arc
{
public:

  //! Default constructor
  Arc(GlobalIndexType id) : mVertices(1,id) {}

  //! Copy constructor
  Arc(const Arc& a) : mVertices(a.mVertices) {}

  //! Return the current size of the arc
  GlobalIndexType size() const {return mVertices.size();}

  //! The set of vertices of this arc sorted in descending order
  std::vector<GlobalIndexType> mVertices;
};



//! Simple graph structure to store a merge tree
class MergeTree
{
public:

	//! Default constructor
	MergeTree();

	//! Destructor
	~MergeTree() {}

	//! Return the number of nodes/arcs
	LocalIndexType size() const {return (LocalIndexType)mNodes.size();}

  //! Return a reference to the i'th node
  Node& node(LocalIndexType i) {return mNodes[i];}

  //! Return a reference to the i'th node
  const Node& node(LocalIndexType i) const {return mNodes[i];}

	//! Return the minimum
	FunctionType minimum() const {return mMinimum;}

	//! Return the maximum
	FunctionType maximum() const {return mMaximum;}

	//! Set the minimum
	void minimum(FunctionType f) {mMinimum = f;}

	//! Set the maximum
	void maximum(FunctionType f) {mMaximum = f;}

	//! Add a critical point which adds both a new node and the corresponding arc
	LocalIndexType addCriticalPoint(GlobalIndexType id);

  //! Add edge
  int addEdge(LocalIndexType up, LocalIndexType down);

  //! Remove an edge
  int removeEdge(LocalIndexType up, LocalIndexType down);

	//! Add the given vertex to the arc with the given label
	int addVertex(GlobalIndexType v, LocalIndexType label);

	//! Split all arcs with more than n vertices
	int splitBySize(LocalIndexType n);

	//! Split all arcs with length longer than l
	int splitByLength(FunctionType l);

	//! Construct a feature by assembling all vertices that belong to it
	void constructFeature(LocalIndexType label, std::vector<GlobalIndexType>& feature) const;

	//! Inflate the metric values
	int inflate();

	//! Deflate the metric values
	int deflate();

private:

	//! The vector of nodes
	std::vector<Node> mNodes;

	//! The vector of arcs
	std::vector<Arc> mArcs;

	//! The highest function value in the tree
	FunctionType mMaximum;

	//! The lowest function value in the tree
	FunctionType mMinimum;


	//! Split the given arc with pos becoming the head of the new arc
	int splitArc(LocalIndexType a,LocalIndexType pos);
};



#endif /* MERGETREE_H_ */
