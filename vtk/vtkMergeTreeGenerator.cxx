/*
 * vtkMergeTree.cxx
 *
 *  Created on: Jan 8, 2016
 *      Author: bremer5
 */

#include <vector>
#include <algorithm>
#include <map>
#include <stack>

#include "vtkObjectFactory.h"
#include "vtkDataObject.h"
#include "vtkDataSet.h"
#include "vtkDataArray.h"
#include "vtkArrayData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkSmartPointer.h"
#include "vtkIntArray.h"
#include "vtkFloatArray.h"
#include "vtkPointSet.h"
#include "vtkPointData.h"
#include "vtkIdList.h"
#include "vtkMergeTreeGenerator.h"
#include "vtkMergeTree.h"
#include "vtkMergeTreeTransformation.h"




//! Standard union-find implementation
/*! This class implements a default union-find structure using an stl::map
 *  to maintain an index map from label values to local indices. This may
 *  not be the fastest implementation but it is convinient to maintain
 *  labels in a spase index space
 */
class UnionFind
{
public:

  //! Default constructor
  UnionFind() {}

  //! Default destructor
  ~UnionFind() {}

  //! Return the current representative of the given label
  vtkIdType rep(vtkIdType id);

  //! Add a label
  void addLabel(vtkIdType label);

  //! Combine the "from" label with the "to" label
  void mergeLabel(vtkIdType from, vtkIdType to);

private:

  //! The current representative of the i'th label
  std::vector<vtkIdType> mLabel;

  //! An index map to convert global label-indices into local mLabel indices
  std::map<vtkIdType,vtkIdType> mIndexMap;
};

vtkIdType UnionFind::rep(vtkIdType id)
{
  vtkIdType local;

  //! Sanity check to make sure we ask only for existsing labels
  assert(mIndexMap.find(id) != mIndexMap.end());

  //! Get the local index of the label in question
  local = mIndexMap.find(id)->second;

  //! Jump "upward" until you find the current representative
  std::stack<vtkIdType> s;
  while (mLabel[local] != id) {
    s.push(local);
    id = mLabel[local];

    assert(mIndexMap.find(id) != mIndexMap.end());
    local = mIndexMap.find(id)->second;
  }

  //! Shortcut the structure
  if (!s.empty()) {
    s.pop();
    while (!s.empty()) {
      mLabel[s.top()] = id;
      s.pop();
    }
  }
  return id;
}

void UnionFind::addLabel(vtkIdType label)
{
  mLabel.push_back(label);
  mIndexMap[label] = mLabel.size()-1;
}

void UnionFind::mergeLabel(vtkIdType from, vtkIdType to)
{
  assert(mIndexMap.find(from) != mIndexMap.end());
  assert(mIndexMap.find(to) != mIndexMap.end());

  // Make sure the "newer" label survives
  assert(from < to);

  mLabel[mIndexMap.find(from)->second] = to;
}

void vtkMergeTreeGenerator::DataSetIterator::initialize(vtkIdType v)
{
  static vtkSmartPointer<vtkIdList> cells = vtkSmartPointer<vtkIdList>::New();
  static vtkSmartPointer<vtkIdList> verts = vtkSmartPointer<vtkIdList>::New();


  Neighbors.clear();

  DataSet->GetPointCells(v,cells);

  for (vtkIdType i=0;i<cells->GetNumberOfIds();i++) {

    DataSet->GetCellPoints(cells->GetId(i),verts);
    for (vtkIdType j=0;j<verts->GetNumberOfIds();j++) {
      if (verts->GetId(j) != v)
        Neighbors.insert(verts->GetId(j));
    }
  }

  nIt = Neighbors.begin();
}





vtkStandardNewMacro(vtkMergeTreeGenerator);

vtkMergeTreeGenerator::vtkMergeTreeGenerator() : vtkPointSetAlgorithm()
{
  this->StoreSegmentation = true;
  this->ComputeMergeTree = true;
  this->Threshold = -10e34;
  this->SplitType = SPLIT_BY_NONE;

  this->SetNumberOfInputPorts(1);
  this->SetNumberOfOutputPorts(1);

}

vtkMergeTreeGenerator::~vtkMergeTreeGenerator()
{
}

int vtkMergeTreeGenerator::FillInputPortInformation(int port, vtkInformation* info)
{
  if (port == 0)
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");

  return 1;
}

int vtkMergeTreeGenerator::FillOutputPortInformation(int port, vtkInformation* info)
{
  switch (port) {
    case 0:
      //info->Set(vtkDataObject::DATA_TYPE_NAME(), "vtkPointSet");
      break;
    case 1:
      //info->Set(vtkDirectedGraph::DATA_TYPE_NAME(), "vtkMergeTree");
      break;
    default:
      break;
  }

  return 1;
}

vtkSmartPointer<vtkSegmentedMergeTree> vtkMergeTreeGenerator::GetTree()
{
  vtkSmartPointer<vtkSegmentedMergeTree> tmp = this->tree;

  this->tree = vtkSmartPointer<vtkSegmentedMergeTree>();

  return tmp;
}


int vtkMergeTreeGenerator::RequestData(
    vtkInformation* vtkNotUsed( request ),
    vtkInformationVector** inputVector,
    vtkInformationVector* outputVector)
{
  vtkIdType i;

  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkInformation *treeInfo = outputVector->GetInformationObject(1);

  // call ExecuteData
  vtkDataSet *data = vtkDataSet::SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));
  vtkPointSet *output = vtkPointSet::SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));



  //vtkMergeTree *tree = vtkMergeTree::SafeDownCast(treeInfo->Get(vtkDataObject::DATA_OBJECT()));

  tree = vtkSmartPointer<vtkSegmentedMergeTree>::New();

  // Create an array for the vertex labels
  vtkSmartPointer<vtkIntArray> labels = vtkSmartPointer<vtkIntArray>::New();
  labels->SetNumberOfComponents(1);
  labels->SetName("Segmentation");
  // For now we will create a full sized array
  labels->SetNumberOfTuples(data->GetNumberOfPoints());



  fprintf(stderr,"Number of points %d %d\n",data->GetNumberOfPoints(),labels->GetNumberOfTuples());

  vtkDataArray* function = data->GetPointData()->GetScalars();
  std::vector<vtkIdType> order;
  std::vector<vtkIdType>::const_iterator oIt;

  //fprintf(stderr,"%d %d %d\n",output->GetNumberOfPoints(),function->GetNumberOfTuples(),data->GetPointData()->GetScalars()->GetNumberOfTuples());

  if (ComputeMergeTree) {

    // Initialize the array (there really should be a faster way) and screen the vertices
    for (i=0;i<data->GetNumberOfPoints();i++) {
      labels->SetTuple1(i,LNULL);

      if (function->GetTuple1(i) >= Threshold)
        order.push_back(i);
    }

    Greater cmp(function);
    std::sort(order.begin(),order.end(),cmp);
  }
  else {
    // Initialize the array (there really should be a faster way) and screen the vertices
    for (i=0;i<labels->GetNumberOfTuples();i++) {
      labels->SetTuple1(i,LNULL);

      if (function->GetTuple1(i) <= Threshold)
        order.push_back(i);
    }

    Greater cmp(function);
    std::sort(order.rbegin(),order.rend(),cmp);
  }

  tree->SetMaximum(function->GetTuple1(order[0]));
  tree->SetMinimum(function->GetTuple1(order.back()));

  fprintf(stderr,"%D Max %f Min %f\n",ComputeMergeTree,tree->GetMaximum(),tree->GetMinimum());

  UnionFind uf;
  vtkIdType neigh_label;
  vtkIdType new_label,tmp;
  vtkIdList* cells;
  NeighborhoodIterator* it = new DataSetIterator(data);
  double fake_id;

  // For all vertices in descending order
  for (oIt=order.begin();oIt!=order.end();oIt++) {

    //fprintf(stderr,"Process vertex %d\n",*oIt);
    // For all neighbors (careful might touch the same neighbor
    // multiple times
    for (it->initialize(*oIt);!it->end();(*it)++) {
      if (labels->GetTuple1(it->id()) != LNULL) { // If the neighbor has already been labeled it is considered higher
        neigh_label = uf.rep(labels->GetTuple1(it->id())); // Find its current active label

        if (labels->GetTuple1(*oIt) == LNULL) {// If this is the first label we see
          labels->SetTuple1(*oIt,neigh_label); // We pass on this label
        }
        else if (neigh_label != labels->GetTuple1(*oIt)) { // If we see a second label *oIt is a saddle

          // If the node corresponding to our current label is not *oIt itself
          // then we have not yet created a critical point for *oIt
          if (tree->GetRep(labels->GetTuple1(*oIt)) != *oIt) {

            // Add a new node into the tree and use its id as label
            new_label = tree->AddNode(*oIt);

            // Now set the pointer for the node corresponding to the current label
            tree->AddEdge(labels->GetTuple1(*oIt),new_label);

            // And pass on the representative
            fake_id = tree->GetRep(labels->GetTuple1(*oIt));
            tree->SetRep(new_label,fake_id);

            // Create a corresponding UF label
            uf.addLabel(new_label);

            // And merge the two labels making sure the later one survives
            uf.mergeLabel(labels->GetTuple1(*oIt),new_label);

            // And update our own label
            labels->SetTuple1(*oIt,new_label);
          }

          // The above if statement took care of the first arc that reached *oIt.
          // Now we take care of the second arc with neigh_label

          // Set the appropriate down pointer
          tree->AddEdge(neigh_label,labels->GetTuple1(*oIt));

          // First, update the representative if necessary. Since we create the
          // labels (and nodes) in order of the sort, the rep id can be used to
          // determine which node is higher.
          if (tree->GetRep(neigh_label) < tree->GetRep(labels->GetTuple1(*oIt))) {
            tree->SetRep(labels->GetTuple1(*oIt),tree->GetRep(neigh_label));
          }

          // Now we merge the labels
          uf.mergeLabel(neigh_label,labels->GetTuple1(*oIt));

        } // end-if we see a second/third/... label
      } // end-if we found a labeled neighbor
    } // end-for all neighbors

    if (labels->GetTuple1(*oIt) == LNULL) { // If we have not found a higher neighbor

      // Add a new node into the tree and use its id as label
      new_label = tree->AddNode(*oIt);

      // Add the label to the UF
      uf.addLabel(new_label);

      // Set the new label
      labels->SetTuple1(*oIt,new_label);
    }

    // Finally, we have the actual label for *oIt and we add it to the respective branch
    tree->AddVertexToBranch(labels->GetTuple1(*oIt),*oIt);

    //fprintf(stderr,"Setting label[%d] to %d\n",*oIt,(int)labels->GetTuple1(*oIt));
  } // end-for all vertices in sorted order

  delete it;

  fprintf(stderr,"Done with MT computation returning 1\n");

  return 1;
}






