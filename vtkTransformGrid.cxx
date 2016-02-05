#include "vtkObjectFactory.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkSmartPointer.h"
#include "vtkRectilinearGrid.h"
#include "vtkUnstructuredGrid.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkDataArray.h"
#include "vtkCellArray.h"
#include "vtkDataSet.h"
#include "vtkImageData.h"

#include "vtkTransformGrid.h"

vtkStandardNewMacro(vtkTransformGrid);

//----------------------------------------------------------------------------
int vtkTransformGrid::FillInputPortInformation(int port, vtkInformation* info)
{
  if (!this->Superclass::FillInputPortInformation(port, info)) {
    return 0;
  }
  if (port == 0) {
    info->Set(vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet");
    return 1;
  }
  return 0;
}

void processRectilinearGrid(vtkRectilinearGrid *input, vtkUnstructuredGrid *output)
{
  vtkDataArray *coords[3] = {input->GetXCoordinates(),
			     input->GetYCoordinates(),
			     input->GetZCoordinates()};
  int nodeRes[3] = {coords[0]->GetNumberOfTuples(),
		    coords[1]->GetNumberOfTuples(),
		    coords[2]->GetNumberOfTuples()};
  
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  const int numPoints = nodeRes[0]*nodeRes[1]*nodeRes[2];
  points->SetNumberOfPoints(numPoints);

  int idx = 0;
  for (int k = 0; k < nodeRes[2]; ++k) {
    double z = coords[2]->GetComponent(k,0);
    for (int j = 0; j < nodeRes[1]; ++j) {
      double y = coords[1]->GetComponent(j,0);
      for (int i = 0; i < nodeRes[0]; ++i) {
	double x = coords[0]->GetComponent(i,0);
	points->SetPoint(idx, x, y, z);
	++idx;
      }
    }
  }
  output->SetPoints(points);

  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  const vtkIdType cellType = VTK_HEXAHEDRON;
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  const int numCells = cellRes[0]*cellRes[1]*cellRes[2];
  cells->SetNumberOfCells(numCells);
  vtkIdTypeArray *cellData = cells->GetData();
  cellData->SetNumberOfComponents(9);
  cellData->SetNumberOfTuples(numCells);

  idx = 0;
  const int n = 8;
  for (int k = 0; k < cellRes[2]; ++k) {
    for (int j = 0; j < cellRes[1]; ++j) {
      for (int i = 0; i < cellRes[0]; ++i) {
	
	int ids[8] = {i   + j    *nodeRes[0] + k    *nodeRes[0]*nodeRes[1],
		      i+1 + j    *nodeRes[0] + k    *nodeRes[0]*nodeRes[1],
		      i+1 + (j+1)*nodeRes[0] + k    *nodeRes[0]*nodeRes[1],
		      i   + (j+1)*nodeRes[0] + k    *nodeRes[0]*nodeRes[1],
		      i   + j    *nodeRes[0] + (k+1)*nodeRes[0]*nodeRes[1],
		      i+1 + j    *nodeRes[0] + (k+1)*nodeRes[0]*nodeRes[1],
		      i+1 + (j+1)*nodeRes[0] + (k+1)*nodeRes[0]*nodeRes[1],
		      i   + (j+1)*nodeRes[0] + (k+1)*nodeRes[0]*nodeRes[1]};
	
	cellData->SetTuple9(idx,n,ids[0],ids[1],ids[2],ids[3],ids[4],ids[5],ids[6],ids[7]);
	++idx;
      }
    }
  }

  output->SetCells(cellType, cells);
  output->GetPointData()->PassData(input->GetPointData());
}

void processImageData(vtkImageData *input, vtkUnstructuredGrid *output)
{
  double origin[3];
  double spacing[3];
  int extent[6];
  int nodeRes[3];

  input->GetOrigin(origin);
  input->GetSpacing(spacing);
  input->GetExtent(extent);
  input->GetDimensions(nodeRes);
  
  vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
  const int numPoints = nodeRes[0]*nodeRes[1]*nodeRes[2];
  points->SetNumberOfPoints(numPoints);
  
  int idx = 0;
  for (int k = extent[4]; k <= extent[5]; ++k) {
    double z = origin[2] + k*spacing[2];
    for (int j = extent[2]; j <= extent[3]; ++j) {
      double y = origin[1] + j*spacing[1];
      for (int i = extent[0]; i <= extent[1]; ++i) {	
	double x = origin[0] + i*spacing[0];
	points->SetPoint(idx, x, y, z);
	++idx;
      }
    }
  }
  output->SetPoints(points);

  vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
  const vtkIdType cellType = VTK_HEXAHEDRON;
  int cellRes[3] = {nodeRes[0]-1, nodeRes[1]-1, nodeRes[2]-1};
  const int numCells = cellRes[0]*cellRes[1]*cellRes[2];
  cells->SetNumberOfCells(numCells);
  vtkIdTypeArray *cellData = cells->GetData();
  cellData->SetNumberOfComponents(9);
  cellData->SetNumberOfTuples(numCells);

  idx = 0;
  const int n = 8;
  for (int k = 0; k < cellRes[2]; ++k) {
    for (int j = 0; j < cellRes[1]; ++j) {
      for (int i = 0; i < cellRes[0]; ++i) {
	
	int ids[8] = {i   + j    *nodeRes[0] + k    *nodeRes[0]*nodeRes[1],
		      i+1 + j    *nodeRes[0] + k    *nodeRes[0]*nodeRes[1],
		      i+1 + (j+1)*nodeRes[0] + k    *nodeRes[0]*nodeRes[1],
		      i   + (j+1)*nodeRes[0] + k    *nodeRes[0]*nodeRes[1],
		      i   + j    *nodeRes[0] + (k+1)*nodeRes[0]*nodeRes[1],
		      i+1 + j    *nodeRes[0] + (k+1)*nodeRes[0]*nodeRes[1],
		      i+1 + (j+1)*nodeRes[0] + (k+1)*nodeRes[0]*nodeRes[1],
		      i   + (j+1)*nodeRes[0] + (k+1)*nodeRes[0]*nodeRes[1]};
	
	cellData->SetTuple9(idx,n,ids[0],ids[1],ids[2],ids[3],ids[4],ids[5],ids[6],ids[7]);
	++idx;
      }
    }
  }

  output->SetCells(cellType, cells);
  output->GetPointData()->PassData(input->GetPointData());
}

//----------------------------------------------------------------------------
int vtkTransformGrid::RequestData(vtkInformation *request,
				  vtkInformationVector **inputVector,
				  vtkInformationVector *outputVector)
{
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet *input = vtkDataSet::
    SafeDownCast(inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);
  vtkUnstructuredGrid *output = vtkUnstructuredGrid::
    SafeDownCast(outInfo->Get(vtkDataObject::DATA_OBJECT()));

  if (input->IsA("vtkRectilinearGrid")) {
    processRectilinearGrid(dynamic_cast<vtkRectilinearGrid*>(input), output);
  }
  if (input->IsA("vtkImageData")) {
    processImageData(dynamic_cast<vtkImageData*>(input), output);
  }
  
  return 1;
}

////////// External Operators /////////////
void vtkTransformGrid::PrintSelf(ostream &os, vtkIndent indent)
{
}
