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

vtkStandardNewMacro(vtkMergeTree);

vtkMergeTree::vtkMergeTree() : vtkMutableDirectedGraph()
{
  Threshold = 0;

  // Create the id array as pedigree Ids
  GetVertexData()->AddArray(vtkSmartPointer<vtkTypeInt64Array>::New());

  // Create an array for representatives
  GetVertexData()->AddArray(vtkSmartPointer<vtkIdType>::New());
}


vtkTypeInt64 vtkMergeTree::AddNode(vtkTypeInt64 id)
{
  vtkIdType index = AddVertex();

  GetVertexData()->GetArray(vtkMergeTree::MESH_ID)->SetTuple1(index,id);
  GetVertexData()->GetArray(vtkMergeTree::REP_ID)->SetTuple1(index,id);
}

vtkMergeTreeNode vtkMergeTree::GetNode(vtkTypeInt64 index)
{
  return vtkMergeTreeNode(this,index);
}

