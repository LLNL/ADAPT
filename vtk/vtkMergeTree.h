/*
 * vtkMergeTree.h
 *
 *  Created on: Jul 28, 2016
 *      Author: bremer5
 */

#ifndef __vtkMergeTree_h
#define __vtkMergeTree_h

#include <stdint.h>
#include <vector>

#include <vtkCommonDataModelModule.h>
#include <vtkMutableDirectedGraph.h>
#include <vtkTypeInt64Array.h>
#include <vtkSmartPointer.h>
#include <vtkDataSetAttributes.h>


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
  virtual vtkIdType AddNode(vtkTypeInt64 id);

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

  // Description
  // Current threshold
  double Threshold;

  // Description
  // The global "maximum" function value
  double Maximum;

  // Description
  // The global "minimum" function value
  double Minimum;

private:

  vtkMergeTree(const vtkMergeTree&) {}  // Not implemented.
  void operator=(const vtkMergeTree&) {}  // Not implemented.


};

class vtkSegmentedMergeTree : public vtkMergeTree
{
public:

  vtkTypeMacro(vtkSegmentedMergeTree,vtkMergeTree);

  static vtkSegmentedMergeTree *New();

  // Construct a new node
  virtual vtkIdType AddNode(vtkTypeInt64 id);

  // Add a vertex to a branch
  void AddVertexToBranch(vtkIdType branch, vtkTypeInt64 id);

  vtkIdType GetNumberOfBranches() const {return Branches.size();}

  const std::vector<vtkTypeInt64>& GetBranch(vtkIdType branch) const {return Branches[branch];}

protected:

    // Description
    // A collection of mesh indices for each branch
    std::vector<std::vector<vtkTypeInt64> > Branches;

    // Default constructor
    vtkSegmentedMergeTree();

    // Destructor
    virtual ~vtkSegmentedMergeTree() {}

private:

    vtkSegmentedMergeTree(const vtkSegmentedMergeTree&) {}  // Not implemented.
  void operator=(const vtkSegmentedMergeTree&) {}  // Not implemented.


};



#endif /* VTK_VTKMERGETREE_H_ */
