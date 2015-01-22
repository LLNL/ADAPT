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

#ifndef LOCALTHRESHOLD_H_
#define LOCALTHRESHOLD_H_

#include <stdint.h>
#include <vector>
#include <algorithm>
#include <set>
#include <cmath>

#include "Comparisons.h"
#include "Neighborhood.h"
#include "UnionFind.h"

/*!
 * Transform the given volume into relevance thresholds based on
 * computing either the merge or split tree of the data.
 * @param data The initial input field
 * @param transform The resulting relevance volume
 * @param dim The dimensions of the data
 * @param greater The metric to used to compute a merge tree
 * @param threshold An optional threshold to make the computation
 * more efficient. All data "below" the threshold is ignored
 * @return 1 if successful 0 otherwise
 */
int compute_local_thresholds(const FunctionType* data, FunctionType* transform, uint32_t dim[3],
                             Comparison<FunctionType>& greater,const FunctionType threshold)
{
  std::vector<IndexType> order; // A sorted list of vertex id's
  typename std::vector<IndexType>::iterator oIt;
  IndexType i;
  IndexType count = dim[0]*dim[1]*dim[2];
  IndexType* label = new IndexType[dim[0]*dim[1]*dim[2]]; // An array storing a label for each vertex
  FunctionType global_min = data[0];
  FunctionType tmp;

  fprintf(stderr,"Screening vertices\n");

  // First we collect all vertex indices above the threshold,
  // initialize the labels volume and keep track of the global
  // minimum
  for (i=0;i<count;i++) {
    label[i] = GNULL;
    if (greater(global_min,data[i]))
        global_min = data[i];
    if (greater(data[i],threshold)) {
      order.push_back(i);
    }
    //if (data[i] != 0)
    //  fprintf(stderr,"%f\n",data[i]);
  }

  fprintf(stderr,"Sorting\n");

  // Sort all the vertices above the threshold by descending order.
  IndexComp<IndexType,FunctionType> sort_comp(data,greater);
  std::sort(order.begin(),order.end(),sort_comp);

  // Initialize the neighborhood
  FullNeighborhood<IndexType> neighborhood(dim);
  typename FullNeighborhood<IndexType>::iterator it;
  typename std::set<IndexType>::iterator sIt;
  IndexType rep;
  UnionFind<IndexType> uf;

  // A map keeping track of the highest maximum for each subtree
  std::map<IndexType,FunctionType> local_maxima;

  int32_t progress = 0;
  int32_t next = 1;
  fprintf(stderr,"Processing  %03d%%\r",0);

  // For all vertices in descending order
  for (oIt=order.begin();oIt!=order.end();oIt++) {
    if (100*progress/order.size() >= next) {
      fprintf(stderr,"Processing  %03ld%%\r",100*progress/order.size());
      next++;
    }

    //fprintf(stderr,"%d of %d\n",progress,order.size());
    progress++;
    //if (*oIt == 612)
    //  fprintf(stderr,"%d %f\n",*oIt,data[*oIt]);

    // For all neighbors
    for (it=neighborhood.begin(*oIt);it!=neighborhood.end(*oIt);it++) {
      //fprintf(stderr,"\t%d\n",*it);
      if (label[*it] != GNULL) { // If the neighbor has already been labeled it is considered higher
        rep = uf.rep(label[*it]); // Find its current active label

        if (label[*oIt] == GNULL) // If this is the first label we see
          label[*oIt] = rep; // We pass on this label
        else if (rep != label[*oIt]) { // If we see a second label *oIt is a saddle

          if (label[*oIt] != *oIt) {// If we have not yet created a label for the saddle
            uf.addLabel(*oIt); // Create one
            local_maxima[*oIt] = data[*oIt]; // Initialize the local maximum

            uf.mergeLabel(label[*oIt],*oIt); // Merge the current label with the new one

            if (greater(local_maxima[label[*oIt]],local_maxima[*oIt]))
              local_maxima[*oIt] = local_maxima[label[*oIt]]; // Compute the actual maximum

            label[*oIt] = *oIt;// Set the new label
          }

          // Now merge the second label as well
          uf.mergeLabel(rep,*oIt);
          if (greater(local_maxima[rep],local_maxima[*oIt]))
            local_maxima[*oIt] = local_maxima[rep]; // Compute the actual maximum
        }
      }
    }

    if (label[*oIt] == GNULL) { // If we didn't have any valid neighbors
      uf.addLabel(*oIt); // We make a new label as maximum
      label[*oIt] = *oIt; // Set the label
      local_maxima[*oIt] = data[*oIt]; // Initialize the local maxima
    }

    // Finally, we compute the relevance values on the fly. Relevance
    // is the distance of the current vertex to its highest maximum
    // in the subtree divided by the height of its subtree
    tmp = local_maxima.find(label[*oIt])->second;
    transform[*oIt] = 1 - fabs(tmp - data[*oIt]) / fabs(tmp - global_min);
    //fprintf(stderr,"\t %f\n",transform[*oIt]);
  }


  fprintf(stderr,"Processing  100%%\n");

  delete[] label;

  return 1;

}





#endif /* LOCALTHRESHOLD_H_ */
