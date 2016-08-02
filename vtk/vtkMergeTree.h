/*
 * vtkMergeTree.h
 *
 *  Created on: Jul 28, 2016
 *      Author: bremer5
 */

#ifndef VTK_VTKMERGETREE_H_
#define VTK_VTKMERGETREE_H_

#include <stdint.h>
#include "vtkCommonDataModelModule.h"
#include "vtkMutableDirectedGraph.h"
#include "vtkTypeInt64Array.h"
#include "vtkSmartPointer.h"
#include "vtkDataSetAttributes.h"

class vtkMergeTreeNode;


class VTKCOMMONDATAMODEL_EXPORT vtkMergeTree : public vtkMutableDirectedGraph
{
public:

  enum MergeTreeAttributes {
    MESH_ID = 0,
    REP_ID = 1,
  };

  vtkTypeMacro(vtkMergeTree,vtkMutableDirectedGraph);

  static vtkMergeTree *New();

  vtkGetMacro(Threshold,double);
  vtkSetMacro(Threshold,double);

  vtkGetMacro(Maximum,double);
  vtkSetMacro(Maximum,double);

  vtkGetMacro(Minimum,double);
  vtkSetMacro(Minimum,double);

  // Construct a new node
  vtkIdType AddNode(vtkTypeInt64 id);

  // Return a node for traversal
  vtkMergeTreeNode GetNode(vtkIdType index);

  // Return the representative for the given node
  vtkTypeInt64 GetId(vtkIdType index) {return GetVertexData()->GetArray(vtkMergeTree::MESH_ID)->GetTuple1(index);}

  // Return the representative for the given node
  vtkTypeInt64 GetRep(vtkIdType index) {return GetVertexData()->GetArray(vtkMergeTree::REP_ID)->GetTuple1(index);}

  // Set the representative for the given node
  void SetRep(vtkIdType index, vtkTypeInt64 rep) {return GetVertexData()->GetArray(vtkMergeTree::REP_ID)->SetTuple1(index,rep);}

protected:

  // Default constructor
  vtkMergeTree();

  // Destructor
  virtual ~vtkMergeTree() {}


private:

  // Description
  // Current threshold
  double Threshold;

  // Description
  // The global "maximum" function value
  double Maximum;

  // Description
  // The global "minimum" function value
  double Minimum;

  vtkMergeTree(const vtkMergeTree&) {}  // Not implemented.
  void operator=(const vtkMergeTree&) {}  // Not implemented.


};

// Description
// A MergeTreeVertex is a convenient way to traverse a merge tree
class vtkMergeTreeNode
{
public:

  // Description
  // Constructor that sets empty fields and creates an invalid vertex
  vtkMergeTreeNode();

  // Description
  // Constructor that creates a valid vertex
  vtkMergeTreeNode(vtkSmartPointer<vtkMergeTree> tree, vtkIdType index) : Tree(tree), Index(index) {}

  // Description
  // Copy-constructor
  vtkMergeTreeNode(const vtkMergeTreeNode& node) : Tree(node.Tree), Index(node.Index) {}

  // Description
  // Assignment operator
  vtkMergeTreeNode* operator=(const vtkMergeTreeNode& node);

  // Return the index of the vertex
  vtkIdType index() const {return Index;}

  // Return the vertex id (the index of the corresponding mesh vertex)
  vtkTypeInt64 id() const {return Tree->GetVertexData()->GetArray(vtkMergeTree::MESH_ID)->GetTuple1(Index);}

  // Return the representative of the current node
  vtkTypeInt64 rep() const {return Tree->GetVertexData()->GetArray(vtkMergeTree::REP_ID)->GetTuple1(Index);}

  // Set the representative
  void rep(vtkTypeInt64 id) {Tree->GetVertexData()->GetArray(vtkMergeTree::REP_ID)->SetTuple1(Index,id);}

private:

  // Description
  // Pointer to the local tree this vertex refers to
  vtkSmartPointer<vtkMergeTree> Tree;

  // Description
  // The index in the tree of this vertex
  vtkIdType Index;
};



#endif /* VTK_VTKMERGETREE_H_ */
