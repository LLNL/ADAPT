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

#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <stack>
#include <set>


#include "Definitions.h"
#include "Comparisons.h"
#include "FullNeighborhood.h"
#include "MergeTree.h"
#include "MTAlgorithm.h"
#include "Relevance.h"
#include "LocalThreshold.h"
#include "Threshold.h"
#include "ManPage.h"

#include "TopologyFileParser/DataHandle.h"
#include "TopologyFileParser/ValueElement.h"
#include "TopologyFileParser/FeatureElement.h"
#include "TopologyFileParser/ClanHandle.h"
#include "TopologyFileParser/StatHandle.h"
#include "TopologyFileParser/SimplificationHandle.h"

//!Number of available input options (size of gOptions)1
#define NUM_OPTIONS 9

//!Array with the list of all available input options
static const char* gOptions[NUM_OPTIONS] = {
    "--help",

    "--i",
    "--o",

    "--dim",

    "--tree-type",
    "--threshold",
    "--split-type",
    "--split",
    "--metric",
};

//! Name of the input file
const char* gInputFileName = NULL;

//! Name of the output file
const char* gOutputFileName = "output";

//! Global array of data (will be allocated to gDim[0]*gDim[1]*gDim[2]
FunctionType* gData = NULL;

//! Global array of dimensions
GlobalIndexType gDim[3] = {0,0,0};

//! Tree type 0 (merge tree), 1 (split tree)
int gTreeType = 0;

//! The lower (merge tree) or upper (split tree) threshold
FunctionType gThreshold = 0;

//! Bool to indicate whether the threshold has been set
bool gThresholdSet = false;

//! Number of different split types
#define NUM_SPLIT_TYPES 2
//! List of available split types
static const char* gSplitTypeOptions[NUM_SPLIT_TYPES] = {
    "length",
    "size",
};
//! Enum of split types
enum SplitType {
  SPLIT_LENGTH = 0,
  SPLIT_SIZE = 1,
};

//! The split type
SplitType gSplitType = SPLIT_LENGTH;

//! The splitting threshold
FunctionType gSplitLimit = -1;

//! Number of metrics
#define NUM_METRIC_TYPES 3
//! List of available metrics
static const char* gMetricTypeOptions[NUM_METRIC_TYPES] = {
    "threshold",
    "relevance",
    "local",
};
//! Enum of metrics
enum MetricType {
  METRIC_THRESHOLD = 0,
  METRIC_RELEVANCE = 1,
  METRIC_LOCAL = 2,
};
//! The metric used
MetricType gMetric = METRIC_THRESHOLD;

using namespace TopologyFileFormat;


void accumulateVolume(const MergeTree& tree, LocalIndexType root, Data<FunctionType>& data)
{
  std::stack<LocalIndexType> front;
  LocalIndexType top;
  std::set<LocalIndexType> processed;

  front.push(root);

  while (!front.empty()) {

    top = front.top();

    // If this is the first time we see this node
    if (processed.find(top) == processed.end()) {

      LocalIndexType up = tree.node(top).up();

      if (up == LNULL) {// If this is a leaf the volume it correct
        front.pop(); // and there is nothing to do
      }
      else { // Otherwise

        do { // add all the children
          front.push(up);
          up = tree.node(up).next();
        } while (up != tree.node(top).up());

        // Remember that we are done with this node
        processed.insert(top);
      }
    }
    else { // If we have already processed this node then the volume of my children is correct

      assert (tree.node(top).up() != LNULL); // We should not get here for a leaf

      LocalIndexType up = tree.node(top).up();
      do { // add the volume from all children

        data[top] += data[up];
        up = tree.node(up).next();
      } while (up != tree.node(top).up());

      front.pop();
    }
  }

  return;
}


/*! \brief Parse the command line input.
 *
 * This function parses the command line input containing the various
 * execution options and specifies the corresponding global variables used
 * during execution accordingly.
 * Note: The available execution options are defined in gOptions.
 * \param argc : The number of input arguments. (As given to main(...)).
 * \param argv : Array of lengths argc containing all input arguments.
 *               (As given to main(...)).
 * \return int : 0 in case of error and 1 in case of successs
 */
int parse_command_line(int argc, const char** argv)
{
  int i,j,option;

  for (i=1;i<argc;i++) {
    option = -1;
    for (j=0; j < NUM_OPTIONS;j++) {
      if(strcmp(gOptions[j],argv[i])==0)
        option= j;
    }

    switch (option) {

    case -1:  // Wrong input parameter
      fprintf(stderr,"\nError: Wrong input parameter \"%s\"\nTry %s --help\n\n",argv[i],argv[0]);
      return 0;
    case 0:   // --help
      print_help(stdout,argv[0]);
      return 0;
    case 1: // --i
      gInputFileName = argv[++i];
      break;
    case 2: // --o
      gOutputFileName = argv[++i];
      break;
    case 3: // --dim
      gDim[0] = atoi(argv[++i]);
      gDim[1] = atoi(argv[++i]);
      gDim[2] = atoi(argv[++i]);
      break;
    case 4: // --tree-type
      gTreeType = atoi(argv[++i]);
      break;
    case 5: // --threshold
      gThreshold = (FunctionType)atof(argv[++i]);
      break;
    case 6: // --split-type
      i++;
      for (j=0; j < NUM_SPLIT_TYPES;j++) {
        if(strcmp(gSplitTypeOptions[j],argv[i])==0) {
          gSplitType = (SplitType)j;
          break;
        }
      }
      if (j == NUM_SPLIT_TYPES) {
        fprintf(stderr,"Sorry, the split type \"%s\"is not recognized .....\n",argv[i]);
        return 0;
      }
      break;
    case 7: // --split
      gSplitLimit = (FunctionType)atof(argv[++i]);
      break;
    case 8: // --metric
      i++;
      for (j=0; j < NUM_METRIC_TYPES;j++) {
        if(strcmp(gMetricTypeOptions[j],argv[i])==0) {
          gMetric = (MetricType)j;
          break;
        }
      }
      if (j == NUM_METRIC_TYPES) {
        fprintf(stderr,"Sorry, the metric type \"%s\"is not recognized .....\n",argv[i]);
        return 0;
      }
      break;
    default:
      return 0;
    }
  }

  return 1;
}

int main(int argc, const char** argv)
{
  //Parse the command line input and define the execution settings
  if ((argc == 1) || (parse_command_line(argc,argv) == 0)) {
    print_help(stdout,argv[0]);
    return 0;
  }

  GlobalIndexType size = gDim[0]*gDim[1]*gDim[2];


  Metric* metric = NULL;
  switch (gMetric) {
    case METRIC_THRESHOLD: {
      metric = new Threshold();
      break;
    }
    case METRIC_RELEVANCE: {
      metric = new Relevance();
      break;
    }
    case METRIC_LOCAL: {
      metric = new LocalThreshold();
      break;
    }
    default:
      fprintf(stderr,"Error, unkonwn metric\n");
      return 0;
  }



  gData = new FunctionType[size];
  LocalIndexType* labels = new LocalIndexType[size];

  FILE* input = NULL;
  if (gInputFileName != NULL)
    input = fopen(gInputFileName,"rb");
  else {
    fprintf(stderr,"Error, no input filename given\n");
    return 0;
  }


  // Some systems do not gracefully handle very large files so we read
  // things by plane instead
  for (int i=0;i<gDim[2];i++)
    fread(gData+i*gDim[0]*gDim[1],sizeof(FunctionType),gDim[0]*gDim[1],input);
  fclose(input);


  MergeTree tree;
  FullNeighborhood neighborhood(gDim);
  bool augmented = metric->explicitArcs();

  if (gTreeType == 0) {
    MergeTreeComp comp;
    merge_tree_sorted_sweep(comp,neighborhood,gThreshold,tree,true,labels);
  }
  else {
    SplitTreeComp comp;
    merge_tree_sorted_sweep(comp,neighborhood,gThreshold,tree,true,labels);
  }

  LocalIndexType label;


  // For what we need we require explicit minima. So we take the last vertex in
  // root branch a make it into a critical point
  for (LocalIndexType i=0;i<tree.size();i++)  {

    // If the branch has a minimum to create and is a root branch
    if ((tree.arc(i).size() > 1) && (tree.node(i).down() == LNULL)) {

      // Create a new arc
      label = tree.addCriticalPoint(tree.arc(i).mVertices.back());

      // Attach it to the bottom
      tree.addEdge(i,label);

      /// Remove it from its current arc
      tree.arc(i).mVertices.pop_back();

      // Make sure that we pass on the representative
      tree.node(label).rep(tree.node(i).rep());
    }

  }



  // Now we potentially want to split the tree
  if (gSplitLimit > 0) { // FOr now assume we have no need for a negative split metric
    switch (gSplitType) {
      case SPLIT_LENGTH:
        tree.splitByLength(gSplitLimit);
        break;
      case SPLIT_SIZE:
        tree.splitBySize((LocalIndexType)gSplitLimit);
        break;
    }

  }


  metric->initialize(gData,&tree);

  // Evaluate the metric on all critical points
  for (LocalIndexType i=0;i<tree.size();i++)
    tree.node(i).metric(metric->eval(tree.node(i).index(),i));


  // Compute the volume
  Data<FunctionType> volume(tree.size());
  for (LocalIndexType i=0;i<tree.size();i++)
    volume[i] = tree.arc(i).size();

  // Accumulate the volume through the tree
  for (LocalIndexType i=0;i<tree.size();i++) {

    if (tree.node(i).down() == LNULL)  // For all roots
      accumulateVolume(tree,i,volume);
  }





  // Prepare the output
  GlobalIndexType progress = 0;
  GlobalIndexType next = 0;

  FILE* output;


  Data<FeatureElement> features;


  FunctionType life[2];
  FunctionType low,high;

  // Initialize the bounds
  low = 10e20;
  high = -10e20;

  // Create enough FeatureElements
  features.resize(tree.size(),FeatureElement(SINGLE_REPRESENTATIVE));


  for (LocalIndexType i=0;i<tree.size();i++) {

    life[1] = tree.node(i).metric();
    if (tree.node(i).down() == LNULL) { // local minimum
      life[0] = tree.node(i).metric();
    }
    else
      life[0] = tree.node(tree.node(i).down()).metric();


    if (life[0] > life[1])
      std::swap(life[0],life[1]);

    features[i].addLink(tree.node(i).down());


    if (gTreeType == 0)
      features[i].direction(0);
    else
      features[i].direction(1);


    features[i].lifeTime(life[0],life[1]);

    low = std::min(low,life[0]);
    high = std::max(high,life[1]);
  }

  fprintf(stderr,"GlobalIndexType %d ... LocalIndexType  %d \n",sizeof(GlobalIndexType),sizeof(LocalIndexType));

  // Create a handle to a family
  ClanHandle clan;
  FamilyHandle family;

  // Set the overall function range of the clan
  family.range(tree.minimum(),tree.maximum());

  // Now assemble the simplification sequence
  SimplificationHandle simp;
  FunctionType range[2];

  // Set the features as data from the simplification handle
  simp.setData(&features);

  // Set the range
  simp.setRange(low,high);

  // For the moment we assume that we have merge or split trees in which case we
  // always only have a single represetative
  simp.fileType(SINGLE_REPRESENTATIVE);

  switch (gMetric) {
    case METRIC_THRESHOLD: {
      simp.metric("Threshold");
      break;
    }
    case METRIC_RELEVANCE: {
      simp.metric("Relevance");
      break;
    }
    case METRIC_LOCAL: {
      simp.metric("LocalThreshold");
      break;
    }
    default:
      fprintf(stderr,"Error, unkonwn metric\n");
      return 0;
  }

  // Set the encoding
  simp.encoding(true);

  // Add the simplification handle to the family
  family.add(simp);


  // Create the volume handle
  StatHandle volume_handle;

  volume_handle.aggregated(true);
  volume_handle.stat("vertexCount");
  volume_handle.species("xray");
  volume_handle.encoding(true);
  volume_handle.setData(&volume);

  family.add(volume_handle);



  // Finally, we attach the family to the clan
  clan.add(family);


  // and write the file
  char full_name[200];
  sprintf(full_name,"%s.family",gOutputFileName);
  clan.write(full_name);


  ClanHandle seg_clan;
  FamilyHandle seg_family;

  // Set the overall function range of the clan
  seg_family.range(tree.minimum(),tree.maximum());

  SegmentationHandle seg_handle;
  seg_handle.domainType(REGULAR_GRID);

  char descriptor[100];
  sprintf(descriptor,"3 %d %d %d",gDim[0],gDim[1],gDim[2]);
  seg_handle.domainDescription(std::string(descriptor));
  seg_handle.encoding(false);

  std::vector<std::vector<GlobalIndexType> > segmentation(tree.size());

  for (LocalIndexType i=0;i<tree.size();i++) {
    segmentation[i] = std::move(tree.arc(i).mVertices);
    //segmentation[i] = tree.arc(i).mVertices;
  }

  seg_handle.setSegmentation(&segmentation);

  seg_family.add(seg_handle);
  seg_clan.add(seg_family);

  // and write the file
  sprintf(full_name,"%s.seg",gOutputFileName);
  seg_clan.write(full_name);



  delete[] gData;
  delete[] labels;
  delete metric;

  return 0;
}



