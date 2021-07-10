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

#ifndef METRIC_H
#define METRIC_H

#include <assert.h>
#include <cstddef.h>
#include <map>

#include "Definitions.h"
#include "MergeTree.h"

//! The basclass for all metrics
class Metric
{
public:

  //! Constructor indicating whether this metric will require explicit arcs
  Metric(bool explicit_arcs) : mExplicitArcs(explicit_arcs), mData(NULL), mTree(NULL) {}

  //! Destructor
  virtual ~Metric() {}

  virtual bool explicitArcs() const {return mExplicitArcs;}

  //! Return the default value
  virtual FunctionType fillValue() const {return mDefault;}

  //! Initialize the pointer to the data and the active tree
  virtual void initialize(const FunctionType* data, const MergeTree* tree) {
    mData = data;
    mTree = tree;
  }

  //! Evaluate the metric at vertex id with the given label
  virtual FunctionType eval(GlobalIndexType id, LocalIndexType label) const {assert(false);return 0;}

  //! Evaluate the metric for all nodes of the given tree
  virtual int eval(MergeTree& tree) const {assert(false);return 0;}

protected:

  //! Flag indicating whether we need explicit arcs
  const bool mExplicitArcs;

  //! Pointer to the data
  const FunctionType* mData;

  //! Pointer to the corresponding merge tree
  const MergeTree* mTree;

  //! Default value for vertices with no label
  FunctionType mDefault;
};




#endif /* METRIC_H_ */
