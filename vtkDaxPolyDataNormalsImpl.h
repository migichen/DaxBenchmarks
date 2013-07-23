/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkMarchingCubes.cxx

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
#ifndef vtkDaxPolyDataNormalsImpl_H
#define vtkDaxPolyDataNormalsImpl_H

#include <iostream>
#include "vtkCellArray.h"
#include "vtkCharArray.h"
#include "vtkDoubleArray.h"
#include "vtkFloatArray.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkIntArray.h"
#include "vtkLongArray.h"
#include "vtkMarchingCubesTriangleCases.h"
#include "vtkMath.h"
#include "vtkMergePoints.h"
#include "vtkObjectFactory.h"
#include "vtkPoints.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkShortArray.h"
#include "vtkStructuredPoints.h"
#include "vtkUnsignedCharArray.h"
#include "vtkUnsignedIntArray.h"
#include "vtkUnsignedLongArray.h"
#include "vtkUnsignedShortArray.h"
#include "vtkStreamingDemandDrivenPipeline.h"
#include "vtkIncrementalPointLocator.h"
#include "vtkUnstructuredGrid.h"

#include "vtkAcceleratorsDaxModule.h" //required for correct implementation
#include "vtkDaxConfig.h"

#include <dax/cont/Scheduler.h>
#include <dax/cont/ArrayHandle.h>
#include <dax/cont/UnstructuredGrid.h>
//#include <Containers.h>
#include <vtkToDax/DataSetConverters.h>
#include <vtkToDax/DataSetTypeToType.h>
#include <vtkToDax/CellTypeToType.h>
#include <daxToVtk/DataSetConverters.h>
#include <daxToVtk/CellTypeToType.h>

#include "vtkDaxPolyDataNormals.h"

#include "SharedStatus.h"
#include "tlog/tlog.h"
#include "tlog/cp_time.h"

namespace vtkDax { namespace worklet {

struct Normals: public dax::exec::WorkletMapCell
{
typedef void ControlSignature(Topology, Field(In,Point),
                  Field(Out));
typedef _3 ExecutionSignature(_2);

DAX_EXEC_EXPORT
dax::Vector3 operator()(const dax::exec::CellField<
          dax::Vector3,
          dax::CellTagTriangle >& coordinates) const
  {
#ifdef _PROFILE
	static long long pretime=0;
	if (Timer::getTimeMS()-pretime > 10) {
		tlog->logMessage("vtkdax worklet");
		pretime = Timer::getTimeMS();
	}
#endif

  const dax::Vector3 a = coordinates[1] - coordinates[0];
  const dax::Vector3 b = coordinates[2] - coordinates[0];
  return dax::math::Normal( dax::math::Cross(a,b) );
  }

};


} } // namespace


//convert an unstructured grid type
namespace vtkToDax{

template<> struct CellTypeAndDataType<VTK_POLY_DATA,VTK_HEXAHEDRON>{enum{Valid=1};};

template<typename CellTypeToTypeDef> struct DataSetTypeToType<CellTypeToTypeDef,vtkPolyData>
{
  typedef CellTypeToTypeDef CellTypeToType;
  typedef typename CellTypeToTypeDef::DaxCellType DaxCellType;
  enum {VTKDataSetType=VTK_POLY_DATA};
  enum {Valid=(CellTypeAndDataType<VTK_POLY_DATA,CellTypeToTypeDef::VTKCellType>::Valid)};
  typedef dax::cont::UnstructuredGrid<DaxCellType,
          vtkToDax::vtkTopologyContainerTag<CellTypeToType>,
          vtkToDax::vtkPointsContainerTag>
          DaxDataSetType;
};

template<typename VTKDataSetType>
inline typename VTKDataSetType::DaxDataSetType dataSetConverter(
		vtkPolyData* input,
    VTKDataSetType)
  {
  //we convert to a hexahedron unstruction grid
  //this uses a vtkCellArrayContainer to get the needed topology information
  typedef typename VTKDataSetType::DaxDataSetType DataSet;
  typedef typename VTKDataSetType::CellTypeToType CellTypeToType;

  const int NUM_POINTS = VTKDataSetType::CellTypeToType::NUM_POINTS;

  dax::cont::ArrayHandle<dax::Vector3,vtkToDax::vtkPointsContainerTag>
      pointsHandle(vtkToDax::vtkPointsPortal<dax::Vector3>(input->GetPoints(),
                                                           input->GetNumberOfPoints()));

  //
  dax::cont::ArrayHandle<dax::Id,vtkToDax::vtkTopologyContainerTag<CellTypeToType> >
      topoHandle(vtkToDax::vtkTopologyPortal<dax::Id, NUM_POINTS >(input->GetPolys(),
                                              input->GetNumberOfCells()*NUM_POINTS));

  return DataSet(topoHandle,pointsHandle);
  }
} // namespace



vtkStandardNewMacro(vtkDaxPolyDataNormals);

vtkDaxPolyDataNormals::vtkDaxPolyDataNormals() {

}

vtkDaxPolyDataNormals::~vtkDaxPolyDataNormals() {

}




//
// Contouring filter specialized for volumes and "short int" data values.
//
int vtkDaxPolyDataNormals::RequestData(
  vtkInformation *vtkNotUsed(request),
  vtkInformationVector **inputVector,
  vtkInformationVector *outputVector)
{
	static int evtMemHostToDev = tlog->createEventID("host->dev", 255,255,0);
	static int evtMemDevToHost = tlog->createEventID("dev->host", 0, 0, 255);

	printf("vtkDaxPolyDataNormals::RequestData (Time from creation: %lf)\n", this->timer.GetElapsedTime() );

	  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
	  vtkInformation *outInfo = outputVector->GetInformationObject(0);
	  vtkPolyData *input = vtkPolyData::SafeDownCast(
	    inInfo->Get(vtkDataObject::DATA_OBJECT()));

	  vtkPolyData* output = vtkPolyData::SafeDownCast(
	                                  outInfo->Get(vtkDataObject::DATA_OBJECT()));
	  int result = 0;
	  if (input)
	  {
		  //
		  // convert vtk dataset to dax
		  //
		  vtkPolyData *data = input;
		  if (data==NULL) {
			  printf("Data not unstructured grid\n");
		  }

	      typedef vtkToDax::CellTypeToType<vtkTriangle> VTKCellTypeStruct;
	      typedef vtkToDax::DataSetTypeToType<VTKCellTypeStruct, vtkPolyData> DataSetTypeToTypeStruct;

	      typedef DataSetTypeToTypeStruct::DaxDataSetType InputDataSetType;

		  InputDataSetType inputDaxData =
				  vtkToDax::dataSetConverter(data, DataSetTypeToTypeStruct() );

		  dax::cont::Scheduler<> scheduler;

		  dax::cont::Timer<> timer;

		  vtkDax::worklet::Normals normWorklet;
		  dax::cont::ArrayHandle<dax::Vector3, vtkToDax::vtkArrayContainerTag<vtkFloatArray> > normals;

		  //
		  // manually move array from control to execution env
		  //
		  tlog->startEvent( evtMemHostToDev );
		  normals.PrepareForOutput( inputDaxData.GetNumberOfCells() );

		  dax::cont::Timer<> timer1;
		  inputDaxData.GetPointCoordinates().PrepareForInput();
		  SharedStatus::getInstance()->norm_copyPointsToDev_time.push_back( timer1.GetElapsedTime() );

		  timer1.Reset();
		  inputDaxData.GetCellConnections().PrepareForInput();
		  SharedStatus::getInstance()->norm_copyCellsToDev_time.push_back( timer1.GetElapsedTime() );

		  tlog->endEvent( evtMemHostToDev );
		  //end moving data to execution env

		  scheduler.Invoke(normWorklet, inputDaxData, inputDaxData.GetPointCoordinates(), normals);

		  //std::cout << "inputDaxData c: " << inputDaxData.GetNumberOfCells() << std::endl;
		  //std::cout << "inputDaxData p: " << inputDaxData.GetNumberOfPoints() << std::endl;
		  //std::cout << "normals: " << normals.GetNumberOfValues() << std::endl;
		  //std::cout << "normals vtk tuples: " << normals.GetPortalControl().GetVtkData()->GetNumberOfTuples() << std::endl;
		  //std::cout << "normals vtk comp: " << normals.GetPortalControl().GetVtkData()->GetNumberOfComponents() << std::endl;


		  output->ShallowCopy(input);


		  //calling GetPortalControl brings memory back from execution to control env
		  //so we should time this too.
		  tlog->startEvent( evtMemDevToHost );
		  timer1.Reset();
		  vtkDataArray *dataArray = normals.GetPortalControl().GetVtkData();
		  SharedStatus::getInstance()->norm_copyToMem_time.push_back( timer1.GetElapsedTime() );
		  tlog->endEvent( evtMemDevToHost );

		  dataArray->SetName("Normals");
		  output->GetPointData()->AddArray(dataArray);

		  double time = timer.GetElapsedTime();
		  printf("Done with time: %lf\n", time);

		  result = 1;
	 } else
	   	 printf("Point data is empty\n");

#if 0  // for benchmarking we don't ask VTK for help
	  if(!result)
	    {
	    result = this->Superclass::RequestData(NULL,inputVector,outputVector);
	    }
#endif
	 return result;
}


// convenience method
int vtkDaxPolyDataNormals::RequestInformation(vtkInformation* request,
                               vtkInformationVector** inputVector,
                               vtkInformationVector* outputVector)
{
	printf("vtkDaxPolyDataNormals::RequestInformation (Time from creation: %lf)\n", this->timer.GetElapsedTime() );
	Superclass::RequestInformation(request, inputVector, outputVector);
	return 1;
}

// Description:
// This is called by the superclass.
// This is the method you should override.
int vtkDaxPolyDataNormals::RequestUpdateExtent(vtkInformation* request,
                                vtkInformationVector** inputVector,
                                vtkInformationVector* outputVector)
{
	printf("vtkDaxPolyDataNormals::RequestUpdateExtent (Time from creation: %lf)\n", this->timer.GetElapsedTime() );
	Superclass::RequestUpdateExtent(request, inputVector, outputVector);
	return 1;
}

#endif // vtkDaxPolyDataNormalsImpl_H
