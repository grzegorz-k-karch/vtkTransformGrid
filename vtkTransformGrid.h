#ifndef __vtkTransformGrid_h
#define __vtkTransformGrid_h

#include "vtkUnstructuredGridAlgorithm.h"

class vtkTransformGrid : public vtkUnstructuredGridAlgorithm
{
 public:
  static vtkTransformGrid *New();
  vtkTypeMacro(vtkTransformGrid, vtkUnstructuredGridAlgorithm);
  void PrintSelf(ostream &os, vtkIndent indent);

 protected:

  int FillInputPortInformation(int port, vtkInformation* info);
  int RequestData(vtkInformation*, 
		  vtkInformationVector**, 
		  vtkInformationVector*);
};
#endif
