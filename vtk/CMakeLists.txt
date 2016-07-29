
SET(MT_SOURCES
  vtkMergeTree.cxx 
  vtkMergeTreeGenerator.cxx
)

ADD_LIBRARY(vtkMergeTree ${MT_SOURCES})

ADD_EXECUTABLE(merge_tree_test merge_tree_test.cxx)
TARGET_LINK_LIBRARIES(merge_tree_test vtkMergeTree ${VTK_LIBRARIES})



