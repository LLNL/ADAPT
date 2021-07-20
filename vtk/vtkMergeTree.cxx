/*
 * vtkMergeTree.cxx
 *
 *  Created on: Jul 28, 2016
 *      Author: bremer5
 */

#include "vtkMergeTree.h"
#include "vtkObject.h"
#include "vtkObjectFactory.h"
#include "vtkDataObject.h"
#include "vtkIdTypeArray.h"

vtkStandardNewMacro(vtkMergeTree);

vtkMergeTree::vtkMergeTree() : vtkMutableDirectedGraph()
{
  Threshold = 0;

  // Create the id array as pedigree Ids
  GetVertexData()->AddArray(vtkSmartPointer<vtkTypeInt64Array>::New());

  // Create an array for representatives
  GetVertexData()->AddArray(vtkSmartPointer<vtkIdTypeArray>::New());
}


vtkIdType vtkMergeTree::AddNode(vtkTypeInt64 id)
{
  vtkIdType index = AddVertex();

  GetVertexData()->GetArray(vtkMergeTree::MESH_ID)->InsertTuple1(index,id);
  GetVertexData()->GetArray(vtkMergeTree::REP_ID)->InsertTuple1(index,id);

  return index;
}

vtkStandardNewMacro(vtkSegmentedMergeTree);

vtkSegmentedMergeTree::vtkSegmentedMergeTree() : vtkMergeTree()
{
}



vtkIdType vtkSegmentedMergeTree::AddNode(vtkTypeInt64 id)
{
  vtkIdType index = vtkMergeTree::AddNode(id);

  Branches.push_back(std::vector<vtkTypeInt64>());

  assert (Branches.size() == index+1);

  return index;
}


void vtkSegmentedMergeTree::AddVertexToBranch(vtkIdType branch, vtkTypeInt64 id)
{
  Branches[branch].push_back(id);
}
