/*
 * vtkMergeTreeGenerator.h
 *
 *  Created on: Jan 8, 2016
 *      Author: bremer5
 */

#ifndef VTK_VTKMERGETREEGENERATOR_H_
#define VTK_VTKMERGETREEGENERATOR_H_

#include <set>

#include "vtkCommonExecutionModelModule.h" // For export macro
#include "vtkDataSetAlgorithm.h"
#include "vtkPointSetAlgorithm.h"
#include "vtkMergeTree.h"

enum MergeTreeSplitType {
  SPLIT_BY_SIZE = 0,
  SPLIT_BY_LENGTH = 1,
  SPLIT_BY_NONE = 2
};


class VTKCOMMONEXECUTIONMODEL_EXPORT vtkMergeTreeGenerator : public vtkPointSetAlgorithm
{
public:

  // Description
  // Invalid segmentation id
  static const vtkIdType LNULL = -1;

  // Description:
  // Constructs with initial values of zero.
  static vtkMergeTreeGenerator *New();

  vtkTypeMacro(vtkMergeTreeGenerator,vtkPointSetAlgorithm);

  void SetSegmentation(bool seg) {StoreSegmentation = seg;}

  void SetMergeTree(bool mt) {ComputeMergeTree = mt;}

  void SetThreshold(double t) {Threshold = t;}

  void SetSplitType(MergeTreeSplitType t) {SplitType = t;}

  void PrintSelf(ostream& os, vtkIndent indent) {}

  int FillInputPortInformation(int port, vtkInformation* info);

  int FillOutputPortInformation(int port, vtkInformation* info);

  vtkSmartPointer<vtkMergeTree> GetTree();

protected:

  class Greater {

  public:
    vtkDataArray* Function;

    Greater(vtkDataArray* f) : Function(f) {}

    bool operator()(vtkIdType u, vtkIdType v) {return this->Function->GetTuple1(u) > this->Function->GetTuple1(v);}
  };

  // Description
  // Abstract base class for different dataset types
  class NeighborhoodIterator
  {
  public:
    NeighborhoodIterator(vtkDataSet* data) {}
    virtual ~NeighborhoodIterator() {}

    virtual void initialize(vtkIdType v) = 0;
    virtual void operator++(int i) = 0;
    virtual vtkIdType id() = 0;
    virtual bool end() = 0;
  };

  class DataSetIterator : public NeighborhoodIterator
  {
  public:
    DataSetIterator(vtkDataSet* data) : NeighborhoodIterator(data), DataSet(data) {}
    virtual void initialize(vtkIdType v);
    virtual void operator++(int i) {nIt++;}
    virtual vtkIdType id() {return *nIt;}
    virtual bool end() {return (nIt == Neighbors.end());}
  private:
    vtkDataSet* DataSet;
    std::set<vtkIdType> Neighbors;
    std::set<vtkIdType>::iterator nIt;
 };

  // Description
  // Flag to determine whether the segmentation id's should be stored
  bool StoreSegmentation;

  // Description
  // Flag to indicate whether to compute merge (true) or split (false) trees
  bool ComputeMergeTree;

  // Description
  // Threshold to apply before computing the tree
  double Threshold;

  // Description
  // Flag to encode whether and how to split the tree
  MergeTreeSplitType SplitType;

  vtkSmartPointer<vtkMergeTree> tree;

  vtkMergeTreeGenerator();
  virtual ~vtkMergeTreeGenerator();

  virtual int RequestData(vtkInformation* request,
                          vtkInformationVector** inputVector,
                          vtkInformationVector* outputVector);

private:

  vtkMergeTreeGenerator(const vtkMergeTreeGenerator&) {}  // Not implemented.
  void operator=(const vtkMergeTreeGenerator&) {}  // Not implemented.


  NeighborhoodIterator* neighborhood(vtkDataSet* data, vtkIdType v);
};


#endif /* VTK_VTKMERGETREEGENERATOR_H_ */
