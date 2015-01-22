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

#include <vector>
#include <algorithm>
#include <set>
#include <cmath>

#include "MTAlgorithm.h"

extern FunctionType* gData;

extern GlobalIndexType gDim[3];


int merge_tree_sorted_sweep(Comparison& greater,
                            Neighborhood& neighborhood,
                            const FunctionType threshold,
                            MergeTree &tree, bool augmented,
                            LocalIndexType* label)
{
  std::vector<GlobalIndexType> order;
  std::vector<GlobalIndexType>::iterator oIt;

  GlobalIndexType i;
  GlobalIndexType count = gDim[0]*gDim[1]*gDim[2];
  FunctionType low = gData[0];


  fprintf(stderr,"Screening vertices\n");

  // First we collect all vertex indices above the threshold,
  // initialize the labels volume and keep track of the global
  // minimum
  for (i=0;i<count;i++) {
    label[i] = LNULL;
    if (greater(gData[i],threshold))
      order.push_back(i);

    if (greater(low,gData[i]))
      low = gData[i];
    //if (data[i] != 0)
    //  fprintf(stderr,"%f\n",data[i]);
  }

  fprintf(stderr,"Sorting\n");

  // Sort all the vertices above the threshold by descending order.
  IndexComp sort_comp(gData,greater);
  std::sort(order.begin(),order.end(),sort_comp);

  tree.maximum(gData[order[0]]);
  tree.minimum(low);

  // Get a neighborhood iterator
  Neighborhood::iterator it;

  // Create a local union find of labels
  UnionFind uf;
  LocalIndexType neigh_label;
  LocalIndexType new_label;

  // Setup some progress report
  uint32_t progress = 0;
  uint32_t next = 1;
  fprintf(stderr,"Processing  %03d%%\r",0);

  // For all vertices in descending order
  for (oIt=order.begin();oIt!=order.end();oIt++) {
    if (100*progress/order.size() >= next) {
      fprintf(stderr,"Processing  %03ld%%\r",100*progress/order.size());
      next++;
    }
    progress++;

    // For all neighbors
    for (it=neighborhood.begin(*oIt);it!=neighborhood.end(*oIt);it++) {
      if (label[*it] != LNULL) { // If the neighbor has already been labeled it is considered higher
        neigh_label = uf.rep(label[*it]); // Find its current active label

        if (label[*oIt] == LNULL) {// If this is the first label we see
          label[*oIt] = neigh_label; // We pass on this label
        }
        else if (neigh_label != label[*oIt]) { // If we see a second label *oIt is a saddle

          // If the node corresponding to our current label is not *oIt itself
          // then we have not yet created a critical point for *oIt
          if (tree.node(label[*oIt]).index() != *oIt) {

            // Add a new node into the tree and use its id as label
            new_label = tree.addCriticalPoint(*oIt);

            // Now set the pointer for the node corresponding to the current label
            tree.addEdge(label[*oIt],new_label);
            //tree.node(label[*oIt]).down(new_label);

            // And pass on the representative
            tree.node(new_label).rep(tree.node(label[*oIt]).rep());

            // Create a corresponding UF label
            uf.addLabel(new_label);

            // And merge the two labels making sure the later one survives
            uf.mergeLabel(label[*oIt],new_label);

            // And update our own label
            label[*oIt] = new_label;
          }

          // The above if statement took care of the first arc that reached *oIt.
          // Now we take care of the second arc with neigh_label

          // Set the appropriate down pointer
          tree.addEdge(neigh_label,label[*oIt]);
          //tree.node(neigh_label).down(label[*oIt]);

          // First, update the representative if necessary. Since we create the
          // labels (and nodes) in order of the sort, the rep id can be used to
          // determine which node is higher.
          if (tree.node(neigh_label).rep() < tree.node(label[*oIt]).rep()) {
            tree.node(label[*oIt]).rep(tree.node(neigh_label).rep());
          }

          // Now we merge the labels
          uf.mergeLabel(neigh_label,label[*oIt]);

        } // end-if we see a second/third/... label
      } // end-if we found a labeled neighbor
    } // end-for all neighbors

    if (label[*oIt] == LNULL) { // If we have not found a higher neighbor

      // Add a new node into the tree and use its id as label
      new_label = tree.addCriticalPoint(*oIt);

      // Initialize its representative
      tree.node(new_label).rep(new_label);

      // Add the label to the UF
      uf.addLabel(new_label);

      label[*oIt] = new_label;
    }

    // If we need the fully augmented tree
    if (augmented)
      tree.addVertex(*oIt,label[*oIt]);

  } // end-for all vertices in sorted order


  fprintf(stderr,"Processing  100%%\n");

  return 1;
}



