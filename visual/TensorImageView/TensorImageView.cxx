#include <vtkSmartPointer.h>
#include <vtkLine.h>
#include <vtkMath.h>
#include <vtkCellArray.h>
#include <vtkTubeFilter.h>
#include <vtkLineSource.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyDataWriter.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkProperty.h>
#include <vtkMath.h>
#include <vtkPoints.h>
#include <vtkProbeFilter.h>
#include <itkImageToVTKImageFilter.h>
#include <vtkParametricTorus.h>
#include <vtkParametricBoy.h>
#include <vtkParametricConicSpiral.h>
#include <vtkParametricCrossCap.h>
#include <vtkParametricDini.h>
#include <vtkParametricEllipsoid.h>
#include <vtkParametricEnneper.h>
#include <vtkParametricFigure8Klein.h>
#include <vtkParametricKlein.h>
#include <vtkVersion.h>
#include <vtkImageData.h>
#include <vtkPoints.h>
#include <vtkImageDataToPointSet.h>
#include <vtkSmartPointer.h>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkImageReader2.h>
#include <vtkImageReader2Factory.h>
#include <vtkStructuredGrid.h>
#include <vtkParametricMobius.h>
#include <vtkParametricRandomHills.h>
#include <vtkParametricRoman.h>
#include <vtkParametricSpline.h>
#include <vtkParametricSuperEllipsoid.h>
#include <vtkParametricSuperToroid.h>
#include <vtkParametricTorus.h>
#include <vtkVersion.h>
#include <vtkSmartPointer.h>
#include <vtkImageData.h>
#include <vtkImageMapper3D.h>
#include <vtkExtractVOI.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkImageCast.h>
#include <vtkImageMandelbrotSource.h>
#include <vtkImageData.h>
#include <vtkMetaImageReader.h>
#include <vtkSmartPointer.h>
#include <vtkInteractorStyleImage.h>
#include <vtkRenderer.h>
#include <vtkImageActor.h>
#include <vtkImageMapper3D.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkImageMapToColors.h>
#include "vtkBoxWidget.h"
#include "vtkCamera.h"
#include "vtkCommand.h"
#include "vtkColorTransferFunction.h"
#include "vtkDICOMImageReader.h"
#include "vtkImageData.h"
#include "vtkImageResample.h"
#include "vtkMetaImageReader.h"
#include "vtkPiecewiseFunction.h"
#include "vtkPlanes.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkVolume.h"
#include "vtkVolumeProperty.h"
#include "vtkXMLImageDataReader.h"
#include "vtkSmartVolumeMapper.h"
#include <vtkVersion.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkRenderWindow.h>
#include <vtkRenderer.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkPolyData.h>
#include <vtkSphereSource.h>
#include <vtkOrientationMarkerWidget.h>
#include <vtkAxesActor.h>
#include <vtkPropAssembly.h>
#include <vtkSmartPointer.h>
#include <vtkActor.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>
#include <vtkMath.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPointSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkPolyData.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkSmartPointer.h>
#include <vtkCubeSource.h>
#include <vtkTensorGlyph.h>
#include <vtkVertexGlyphFilter.h>
#include <vtkSphereSource.h>
#include <vtkCubeSource.h>
#include <sstream>
#include <vtkStructuredPointsReader.h>
#include <fstream>
#include <iostream>
#include <vtkGlyph3D.h>
#include "itkImageRegion.h"
#include "itkIndex.h"
#include "itkSize.h"
#include <math.h>
#include "itkMetaImageIO.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkStatisticsImageFilter.h"
#include "itkMaskImageFilter.h"
#include "itkResampleImageFilter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkAffineTransform.h"
#include "itkCastImageFilter.h"
#include "itkNumericSeriesFileNames.h"
#include "itkCenteredTransformInitializer.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkTranslationTransform.h"
#include "itkAndImageFilter.h"
#include "itkImageMomentsCalculator.h"
#include <iostream>
#include <fstream>
#include <itkExtractImageFilter.h>
#include "itkBinaryThresholdImageFilter.h"
#include <vtkArrowSource.h>
#include <vtkScalarBarActor.h>
#include <vtkParametricFunctionSource.h>
#include <vtkPolyDataConnectivityFilter.h>
#include <vtkVersion.h>
#include <vtkActor.h>
#include <vtkFloatArray.h>
#include <vtkLookupTable.h>
#include <vtkPointData.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkScalarBarActor.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkMarchingCubes.h>
#include <vtkImageActor.h>
#include <vtkJPEGReader.h>
#include "itkExtractImageFilter.h"
#include <vtkPNGWriter.h>
#include <vtkRenderLargeImage.h>
#include <itkChangeInformationImageFilter.h>
#include "itkRescaleIntensityImageFilter.h"
#include <vtkWindowToImageFilter.h>
#include <vtkImageReslice.h>
#include <vtkTransform.h>
#include <vtkConeSource.h>
#include <itkVTKImageToImageFilter.h>
#include <vtkImageLuminance.h>
#include <vtkGPUVolumeRayCastMapper.h>
#include <vtkArrowSource.h>
#include <iostream>
#include <iomanip>
#include <itkBinaryMask3DMeshSource.h>
#include <itkMesh.h>
#include <itkMeshSpatialObject.h>
#include <vtkPolyData.h>
#include <vtkCellArray.h>
#include <vtkSphereSource.h>
#include <vtkCellLinks.h>
#include <itkImage.h>
#include <vtkImageStencil.h>
#include <vtkSTLWriter.h>
#include <vtkImageExport.h>
#include <vtkPlaneSource.h>
#include <vtkPolyDataMapper.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkRenderer.h>
//#include <vtkStreamLine.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkStructuredGridOutlineFilter.h>
#include <vtkStructuredGrid.h>
#include <vtkProperty.h>
 #include <iostream>
 #include <iomanip>
 #include <vtkSmartPointer.h>
 #include <itkSpatialObjectToImageFilter.h>
 #include <itkImageFileWriter.h>
 #include  <vtkPolyDataToImageStencil.h>
#include <vtkMetaImageWriter.h>
#include "vtkHedgeHog.h"

#include <vtkContourFilter.h>
#include <vtkImageGradientMagnitude.h>
#include <vtkAppendFilter.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkSphereSource.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
#include <vtkXMLUnstructuredGridWriter.h>
#include <vtkRTAnalyticSource.h>
#include <vtkXMLImageDataWriter.h>
 #if defined(_MSC_VER)
 #pragma warning ( disable : 4786 )
 #endif

 #ifdef __BORLANDC__
 #define ITK_LEAN_AND_MEAN
 #endif






#define VTI_FILETYPE 1
#define MHA_FILETYPE 2
using namespace std;
	typedef itk::Image<short, 3>  ImageType;
	typedef itk::Image<unsigned char, 3>  ImageType_uchar;
	 #ifndef vtkFloatingPointType
 #define vtkFloatingPointType float
 #endif

typedef  short        InputPixelType;
  typedef  short        OutputPixelType;
  typedef itk::Image< InputPixelType,  3 >    InputImageType;
  typedef itk::Image< OutputPixelType, 2 >    OutputImageType;
  typedef itk::ImageFileReader< InputImageType  >  ReaderType;
  typedef itk::ImageFileWriter< OutputImageType >  WriterType;


// Callback for moving the planes from the box widget to the mapper
class vtkBoxWidgetCallback : public vtkCommand
{
public:
  static vtkBoxWidgetCallback *New()
    { return new vtkBoxWidgetCallback; }
  virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
      if (this->Mapper)
        {
        vtkPlanes *planes = vtkPlanes::New();
        widget->GetPlanes(planes);
        this->Mapper->SetClippingPlanes(planes);
        planes->Delete();
        }
    }
  void SetMapper(vtkSmartVolumeMapper* m)
    { this->Mapper = m; }

protected:
  vtkBoxWidgetCallback()
    { this->Mapper = 0; }

  vtkSmartVolumeMapper *Mapper;
};
	typedef itk::Vector<double,3>  VectorType;

bool SetupEnvironmentForDepthPeeling(
  vtkSmartPointer<vtkRenderWindow> renderWindow,
  vtkSmartPointer<vtkRenderer> renderer, int maxNoOfPeels,
  double occlusionRatio)
{
  if (!renderWindow || !renderer)
    return false;
 
  // 1. Use a render window with alpha bits (as initial value is 0 (false)):
  renderWindow->SetAlphaBitPlanes(true);
 
  // 2. Force to not pick a framebuffer with a multisample buffer
  // (as initial value is 8):
  renderWindow->SetMultiSamples(0);
 
  // 3. Choose to use depth peeling (if supported) (initial value is 0 (false)):
  renderer->SetUseDepthPeeling(true);
 
  // 4. Set depth peeling parameters
  // - Set the maximum number of rendering passes (initial value is 4):
  renderer->SetMaximumNumberOfPeels(maxNoOfPeels);
  // - Set the occlusion ratio (initial value is 0.0, exact image):
  renderer->SetOcclusionRatio(occlusionRatio);
 
  return true;
}

bool IsDepthPeelingSupported(vtkSmartPointer<vtkRenderWindow> renderWindow,
                             vtkSmartPointer<vtkRenderer> renderer,
                             bool doItOffScreen)
{
  if (!renderWindow || !renderer)
    {
    return false;
    }
 
  bool success = true;
 
  // Save original renderer / render window state
  bool origOffScreenRendering = renderWindow->GetOffScreenRendering() == 1;
  bool origAlphaBitPlanes = renderWindow->GetAlphaBitPlanes() == 1;
  int origMultiSamples = renderWindow->GetMultiSamples();
  bool origUseDepthPeeling = renderer->GetUseDepthPeeling() == 1;
  int origMaxPeels = renderer->GetMaximumNumberOfPeels();
  double origOcclusionRatio = renderer->GetOcclusionRatio();
 
  // Activate off screen rendering on demand
  renderWindow->SetOffScreenRendering(doItOffScreen);
 
  // Setup environment for depth peeling (with some default parametrization)
  success = success && SetupEnvironmentForDepthPeeling(renderWindow, renderer,
                                                       100, 0.1);
 
  // Do a test render
  renderWindow->Render();
 
  // Check whether depth peeling was used
  success = success && renderer->GetLastRenderingUsedDepthPeeling();
 
  // recover original state
  renderWindow->SetOffScreenRendering(origOffScreenRendering);
  renderWindow->SetAlphaBitPlanes(origAlphaBitPlanes);
  renderWindow->SetMultiSamples(origMultiSamples);
  renderer->SetUseDepthPeeling(origUseDepthPeeling);
  renderer->SetMaximumNumberOfPeels(origMaxPeels);
  renderer->SetOcclusionRatio(origOcclusionRatio);
 
  return success;
}
 

int main(int argc, char *argv[])
{
  vtkRenderer *renderer = vtkRenderer::New();
  renderer->SetBackground(1,1,1);
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(renderer);
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  iren->GetInteractorStyle()->SetDefaultRenderer(renderer);
  vtkAlgorithm *reader1=0;
  vtkAlgorithm *reader2=0;
  vtkAlgorithm *reader3=0;
  vtkImageData *input1=0;
  vtkImageData *input2=0;
  vtkImageData *input3=0;



//typedef itk::ImageFileReader<ImageType> ReaderType;
  ReaderType::Pointer filterchange3 = ReaderType::New();
  ReaderType::Pointer filterchange2 = ReaderType::New();
  ReaderType::Pointer filterchange1 = ReaderType::New();

  
typedef itk::ImageToVTKImageFilter<ImageType> ConnectorType;
  ConnectorType::Pointer metaReader1 = ConnectorType::New();
  ConnectorType::Pointer metaReader2 = ConnectorType::New();
  ConnectorType::Pointer metaReader3 = ConnectorType::New();
 

typedef itk::ChangeInformationImageFilter< ImageType > FilterType1;
  FilterType1::Pointer readeritk3 = FilterType1::New();
  FilterType1::Pointer readeritk2 = FilterType1::New();
  FilterType1::Pointer readeritk1 = FilterType1::New();

    std::string result_file1;
	std::stringstream sstm_file1;
	//sstm_file1 <<argv[1]<< "new_modified_choped.nii";
	sstm_file1 <<argv[1];
	result_file1 = sstm_file1.str();

  filterchange1->SetFileName(result_file1);
  filterchange1->Update();

 	typedef itk::BinaryThresholdImageFilter <ImageType, ImageType> BinaryThresholdImageFilterType;
  BinaryThresholdImageFilterType::Pointer thresholdFilter= BinaryThresholdImageFilterType::New();
  thresholdFilter->SetInput(filterchange1->GetOutput());
  thresholdFilter->SetLowerThreshold(0);
  thresholdFilter->SetUpperThreshold(0);
  thresholdFilter->SetInsideValue(0);
  thresholdFilter->SetOutsideValue(255);
  thresholdFilter->Update();

  ImageType::Pointer rval2 = thresholdFilter->GetOutput();
//  
//      std::string result_bw;
//	std::stringstream sstm_bw;
//	sstm_bw <<argv[2]<< "result_bw.nii";
//	result_bw = sstm_bw.str();
//
//typedef itk::ImageFileWriter<ImageType>  WriterType;
//    WriterType::Pointer writer_red = WriterType::New();
//  writer_red->SetFileName( result_bw );
//  writer_red->SetInput( rval2);
//  writer_red->Update();



  std::cout<<"\n FileName 1:" <<result_file1;


	readeritk1->SetInput(filterchange1->GetOutput());
	readeritk1->SetOutputOrigin(filterchange1->GetOutput()->GetOrigin());
	readeritk1->ChangeOriginOn();
	readeritk1->Update();

	std::cout<<"\n Origin 1:" <<readeritk1->GetOutputOrigin();

	metaReader1->SetInput(readeritk1->GetOutput());
	metaReader1->Update();

	input1=metaReader1->GetOutput();
	std::cout<<"\n Components:"<<input1->GetNumberOfScalarComponents();

vtkImageResample *resample = vtkImageResample::New();
vtkImageResample *resample2 = vtkImageResample::New();
	float reductionFactor=1;

	resample->SetInputData(metaReader1->GetOutput());
	resample->SetAxisMagnificationFactor(0, reductionFactor);
    resample->SetAxisMagnificationFactor(1, reductionFactor);
    resample->SetAxisMagnificationFactor(2, reductionFactor);
	resample->Update();

vtkVolume *volume1 = vtkVolume::New();
vtkGPUVolumeRayCastMapper *mapper1 = vtkGPUVolumeRayCastMapper::New();

	mapper1->SetInputConnection( resample->GetOutputPort() );
	mapper1->Update();


vtkColorTransferFunction *colorFun = vtkColorTransferFunction::New();
vtkPiecewiseFunction *opacityFun = vtkPiecewiseFunction::New();
vtkVolumeProperty *property = vtkVolumeProperty::New();
  
  property->SetIndependentComponents(true);
  property->SetColor( colorFun );
  property->SetScalarOpacity( opacityFun );
  property->SetInterpolationTypeToLinear();
  
  volume1->SetProperty( property );
  volume1->SetMapper( mapper1 );
  
  colorFun->AddHSVSegment(1,1.,1.,1.,1.,0.,1.,1.);
  //colorFun->AddHSVSegment(0.6666,0.6666,1.,1.,1.,0.,1.,1.);
  colorFun->AddRGBPoint( 4071, .83, .66, 1, 0.5, 1 );
  //colorFun->AddRGBPoint( 300, .83, .66, 1, 0.5, 0.5 );

  opacityFun->AddPoint(84, 0); 
  opacityFun->AddPoint(151, .3); 
  opacityFun->AddPoint(255, 0.3);
  
  mapper1->SetBlendModeToMaximumIntensity();

    double origin[3];  
    origin[0] = filterchange1->GetOutput()->GetOrigin()[0];
    origin[1] = filterchange1->GetOutput()->GetOrigin()[1];
    origin[2] = filterchange1->GetOutput()->GetOrigin()[2];

	volume1->SetMapper( mapper1);
	volume1->SetOrigin(origin);
	volume1->Update();

	std::cout<<"\n Origin 1:" <<volume1->GetOrigin()[0]<<"\t"<<volume1->GetOrigin()[1]<<"\t"<<volume1->GetOrigin()[2];
	std::cout<<"\n XRange 1:" <<volume1->GetXRange()[0]<<"\t" <<volume1->GetXRange()[1];
	std::cout<<"\n YRange 1:" <<volume1->GetYRange()[0]<<"\t" <<volume1->GetYRange()[1];
	std::cout<<"\n ZRange 1:" <<volume1->GetZRange()[0]<<"\t" <<volume1->GetZRange()[1];
	std::cout<<"\n ";

    renderer->AddVolume(volume1);

vtkSmartPointer<vtkDoubleArray> tensors_MIL = vtkSmartPointer<vtkDoubleArray>::New();
    tensors_MIL->SetNumberOfTuples(1);
    tensors_MIL->SetNumberOfComponents(9);

vtkSmartPointer<vtkDoubleArray> vector_MIL = vtkSmartPointer<vtkDoubleArray>::New();
    vector_MIL->SetNumberOfTuples(1);
    vector_MIL->SetNumberOfComponents(3);

	vtkSmartPointer<vtkDoubleArray> vector_structuredgrid = vtkSmartPointer<vtkDoubleArray>::New();
    vector_structuredgrid->SetNumberOfTuples(1);
    vector_structuredgrid->SetNumberOfComponents(3);

vtkSmartPointer<vtkDoubleArray> tensors_MIL1 = vtkSmartPointer<vtkDoubleArray>::New();
    tensors_MIL1->SetNumberOfTuples(1);
    tensors_MIL1->SetNumberOfComponents(9);

std::string result_MIL;
std::stringstream sstm_MIL;
	sstm_MIL <<argv[2];
	result_MIL = sstm_MIL.str();


std::string filename_MIL = result_MIL;
std::ifstream fin_MIL(filename_MIL.c_str());


	 cout<<"MILDATA :"<<filename_MIL;


		 std::string line_MIL;
		 std::string line_MIL1;
vtkSmartPointer<vtkPoints> points1_MIL = vtkSmartPointer<vtkPoints>::New();
vtkSmartPointer<vtkPoints> points1_MIL1 = vtkSmartPointer<vtkPoints>::New();
	 int counter_MIL=0;
	 int counter_structuredgrid=0;
	 int counter_MIL1=0;
     
vtkSmartPointer<vtkUnsignedCharArray> colors = 
    vtkSmartPointer<vtkUnsignedCharArray>::New();
	colors->SetName("colors");
	colors->SetNumberOfComponents(3);
	  resample->Update();
	
	//
	//vtkSmartPointer<vtkImageData> whiteImage = vtkSmartPointer<vtkImageData>::New();    
 // whiteImage->DeepCopy(resample->GetOutput());
 // whiteImage->AllocateScalars(VTK_FLOAT, 3);
 //   vtkIdType count = whiteImage->GetNumberOfPoints();
 //   for (vtkIdType i = 0; i < count; ++i)
 //   {
 //   whiteImage->GetPointData()->GetScalars()->SetTuple3(i, 0 ,0, 0);
 //   }
//
//vtkSmartPointer<vtkStructuredGrid> structuredGrid = vtkSmartPointer<vtkStructuredGrid>::New();
//structuredGrid->SetDimensions(25,29,17);
//vtkSmartPointer<vtkPoints> points_structuredGrid = vtkSmartPointer<vtkPoints>::New();
//
//for(int i=0;i<17;i++)
//{
//for(int j=0;j<22;j++)
//{
//for(int k=0;k<11;k++)
//{
//
//	points_structuredGrid->InsertNextPoint(i,j,k);
//	vector_structuredgrid->InsertTuple3(counter_structuredgrid,0,0,0);
//	counter_structuredgrid++;
//}
//}
//}




 //vtkSmartPointer<vtkDoubleArray> weights =vtkSmartPointer<vtkDoubleArray>::New();
 //weights->SetNumberOfValues(516);
 // vtkSmartPointer<vtkDoubleArray> weights1 =vtkSmartPointer<vtkDoubleArray>::New();
 //weights1->SetNumberOfValues(516);
 //vtkIdType cellid;
  while(std::getline(fin_MIL, line_MIL))
    {
    std::stringstream linestream_MIL;
    linestream_MIL << line_MIL;
    //linestream_MIL >> x_MIL>>y_MIL>>z_MIL>>x1_MIL;
    double x1_MIL,x2_MIL,x3_MIL, y1_MIL,y2_MIL,y3_MIL, z1_MIL,z2_MIL,z3_MIL, x_coord,y_coord,z_coord;
    double eval1, eval2, eval3;
    double own_DA, BVTV, VarDA,VarMetJaccom,VarMetRotcom,VarMetRotglo,VarAngJaccom,VarAngRotcom,VarAngRotglo,VarDAJaccom,VarDARotcom,VarDARotglo;
    // read a line in the input file
	  linestream_MIL >> x_coord>>y_coord>>z_coord>> x1_MIL>>x2_MIL>>x3_MIL ;
	  // linestream_MIL >> x_coord>>y_coord>>z_coord>> eval1>> eval2>> eval3 >> x1_MIL>>x2_MIL>>x3_MIL >>y1_MIL >> y2_MIL >>y3_MIL >>z1_MIL>>z2_MIL>>z3_MIL;
    //own_DA = eval3/eval1;
	//x_MIL=x_MIL-3;
	//y_MIL=y_MIL-3;
	//z_MIL=z_MIL-3;
	//cout<<"\n\n\n[x,y,z]:"<<"\t"<<x_MIL<<"\t"<<y_MIL<<"\t"<<z_MIL<<"\t"<<x1_MIL<<"\t"<<x2_MIL<<"\t"<<x3_MIL<<"\t"<<y1_MIL<<"\t"<<y2_MIL<<"\t"<<y3_MIL<<"\t"<<z1_MIL<<"\t"<<z2_MIL<<"\t"<<z3_MIL;
    //points2_MIL->InsertNextPoint(x_MIL,y_MIL,z_MIL);
	//ImageType::IndexType pixelIndex;
	//pixelIndex[0]=(x_MIL/0.082);
	//pixelIndex[1]=(y_MIL/0.082);
	//pixelIndex[2]=(z_MIL/0.082);

	//int values_tensor[3];
	//values_tensor[0]=std::abs(x_MIL/0.082);
	//values_tensor[1]=std::abs(y_MIL/0.082);
	//values_tensor[2]=std::abs(z_MIL/0.082);

	//
	//if(rval2->GetPixel(pixelIndex)>0)
	//{
	//whiteImage->GetPointData()->GetScalars()->SetTuple3(whiteImage->ComputeCellId(values_tensor),x1_MIL,x2_MIL, x3_MIL);
	//points_structuredGrid->InsertNextPoint(x_MIL,y_MIL,z_MIL);
	points1_MIL->InsertNextPoint(x_coord,y_coord,z_coord);
	vector_MIL->InsertTuple3(counter_MIL, x1_MIL,x2_MIL, x3_MIL);
    //tensors_MIL->InsertTuple9(counter_MIL, x1_MIL,x2_MIL, x3_MIL, y1_MIL,y2_MIL, y3_MIL,z1_MIL,z2_MIL, z3_MIL);
	//weights->SetValue(counter_MIL,own_DA);
	//weights1->SetValue(counter_MIL,BVTV);
	counter_MIL++;
	//}

    }
  //weights1->SetName("Metric");
  //weights->SetName("DA");
  fin_MIL.close();


  //structuredGrid->SetPoints(points_structuredGrid);
  //structuredGrid->GetPointData()->SetVectors(vector_MIL);
	int a[3]={0,0,0};
	//whiteImage->ComputeCellId(a);


std::cout<<"\n Counter_MIL:" <<counter_MIL;
std::cout<<"\n  argv[1]:" <<argv[1];
std::cout<<"\n  argv[2]:" <<argv[2];
std::cout<<"\n  argv[3]:" <<argv[3];



 // while(std::getline(fin_MIL1, line_MIL1))
 //   {
 //   std::stringstream linestream_MIL1;
 //   linestream_MIL1 << line_MIL1;
 //   //linestream_MIL >> x_MIL>>y_MIL>>z_MIL>>x1_MIL;
	//double x1_MIL,x2_MIL,x3_MIL,y1_MIL,y2_MIL,y3_MIL,z1_MIL,z2_MIL,z3_MIL,x_MIL,y_MIL,z_MIL;
	//linestream_MIL1 >> x_MIL>>y_MIL>>z_MIL>>x1_MIL>>x2_MIL>>x3_MIL>>y1_MIL>>y2_MIL>>y3_MIL>>z1_MIL>>z2_MIL>>z3_MIL;
	////cout<<"\n\n\n[x,y,z]:"<<"\t"<<x_MIL<<"\t"<<y_MIL<<"\t"<<z_MIL<<"\t"<<x1_MIL<<"\t"<<x2_MIL<<"\t"<<x3_MIL<<"\t"<<y1_MIL<<"\t"<<y2_MIL<<"\t"<<y3_MIL<<"\t"<<z1_MIL<<"\t"<<z2_MIL<<"\t"<<z3_MIL;
 //   //points2_MIL->InsertNextPoint(x_MIL,y_MIL,z_MIL);
	//ImageType::IndexType pixelIndex;
	//pixelIndex[0]=(x_MIL/0.082);
	//pixelIndex[1]=(y_MIL/0.082);
	//pixelIndex[2]=(z_MIL/0.082);
	//if(rval2->GetPixel(pixelIndex)>0)
	//{
	//points1_MIL1->InsertNextPoint(x_MIL,y_MIL,z_MIL);
 //   tensors_MIL1->InsertTuple9(counter_MIL1, x1_MIL,x2_MIL, x3_MIL, y1_MIL,y2_MIL, y3_MIL,z1_MIL,z2_MIL, z3_MIL);
	//counter_MIL1++;
	//}
 //   }
 // fin_MIL1.close();
  //std::cout<<"\n Counter_MIL:" <<counter_MIL1;


  
  //std::cout<<"\n No of compo:" <<whiteImage->GetNumberOfPoints();


 

  //std::cin>>m;

//
vtkSmartPointer<vtkPolyData> polyDataTen_MIL = vtkSmartPointer<vtkPolyData>::New();
polyDataTen_MIL->SetPoints(points1_MIL);
//polyDataTen_MIL->GetPointData()->SetTensors(tensors_MIL);
polyDataTen_MIL->GetPointData()->SetVectors(vector_MIL);
//polyDataTen_MIL->GetPointData()->SetScalars(weights);
//polyDataTen_MIL->GetPointData()->AddArray(weights1);
	
//	
////vtkSmartPointer<vtkPolyData> polyDataTen_MIL1 = vtkSmartPointer<vtkPolyData>::New();
////    polyDataTen_MIL1->SetPoints(points1_MIL1);
////	polyDataTen_MIL1->GetPointData()->SetTensors(tensors_MIL1);
//
vtkSmartPointer<vtkArrowSource> arrowSource = vtkSmartPointer<vtkArrowSource>::New();
arrowSource->SetTipRadius(0.0);
arrowSource->SetShaftRadius(0.05);
arrowSource->Update();
//
//vtkSmartPointer<vtkSphereSource> sphereSource = vtkSmartPointer<vtkSphereSource>::New();
//        sphereSource->SetThetaResolution(36);
//        sphereSource->SetPhiResolution(36);
//        sphereSource->Update();
//
//vtkSmartPointer<vtkCubeSource> CubeSource = vtkSmartPointer<vtkCubeSource>::New();
//		CubeSource->SetXLength(6);
//		CubeSource->SetYLength(6);
//		CubeSource->SetZLength(6);
//		CubeSource->Update();
//
//
////
vtkSmartPointer<vtkGlyph3D> glyph3D = vtkSmartPointer<vtkGlyph3D>::New();
glyph3D->SetInputData(polyDataTen_MIL);
glyph3D->SetSourceConnection(arrowSource->GetOutputPort());
glyph3D->SetScaleFactor(3);
//glyph3D->SetScaleModeToScaleByScalar();
glyph3D->SetVectorModeToUseVector();
glyph3D->Update();
//
//
////vtkSmartPointer<vtkTensorGlyph> tensorGlyph_MIL = vtkSmartPointer<vtkTensorGlyph>::New();
////    tensorGlyph_MIL->SetInputData(polyDataTen_MIL);
////    tensorGlyph_MIL->SetSourceConnection(arrowSource->GetOutputPort());
////	//tensorGlyph_MIL->SetColorModeToScalars();
////    tensorGlyph_MIL->ThreeGlyphsOff();
////	//tensorGlyph_MIL->SetColorMode(1);
////	//tensorGlyph_MIL->ColorGlyphsOn();
////    tensorGlyph_MIL->ExtractEigenvaluesOn();
////	tensorGlyph_MIL->GetExtractEigenvalues();
////	//tensorGlyph_MIL->ScalingOn();
////	tensorGlyph_MIL->SetScaling(1);
////	//tensorGlyph_MIL->SetScaleFactor(3);
////	tensorGlyph_MIL->Update();
//
////
////
//vtkSmartPointer<vtkProbeFilter> probe = vtkSmartPointer<vtkProbeFilter>::New();
//probe->SetInputData(metaReader1->GetOutput());
////probe->SetInputData(resample->GetOutput());
////probe->SetSourceData(glyph3D->GetOutput());
//
//probe->SetSourceData(polyDataTen_MIL);
//   probe->SpatialMatchOn();
//
//   probe->Update();
//
////      resample
//   
//
//  vtkSmartPointer<vtkImageDataToPointSet> imageDataToPointSet =
//    vtkSmartPointer<vtkImageDataToPointSet>::New();
//  imageDataToPointSet->SetInputData(probe->GetOutput());
//  imageDataToPointSet->Update();
 


  //vtkSmartPointer<vtkXMLStructuredGridWriter> writer =
  //  vtkSmartPointer<vtkXMLStructuredGridWriter>::New();
  //writer->SetInputData(structuredGrid);
  //writer->SetFileName(argv[3]);
  //writer->Write();



   //vtkSmartPointer<vtkMetaImageWriter> writer =
   //   vtkSmartPointer<vtkMetaImageWriter>::New();
   //writer->SetInputData(probe->GetOutput());
   //writer->SetFileName(argv[5]);
   //writer->SetRAWFileName(argv[6]);
   //writer->Write();
     std::cout<<"\n\n\n\n Done:";
/*
  vtkSmartPointer<vtkSphereSource> sphereSource1 =
    vtkSmartPointer<vtkSphereSource>::New();
  sphereSource1->SetCenter(0.0, 0.0, 0.0);
  sphereSource1->SetRadius(5.0);
  sphereSource1->Update();

  vtkSmartPointer<vtkAppendFilter> appendFilter =
    vtkSmartPointer<vtkAppendFilter>::New();
  appendFilter->AddInputData(glyph3D->GetOutput());
  appendFilter->Update();


	vtkSmartPointer<vtkUnstructuredGrid> unstructuredGrid =
     vtkSmartPointer<vtkUnstructuredGrid>::New();
  unstructuredGrid->ShallowCopy(appendFilter->GetOutput());*/
 
  //// Write the unstructured grid
  //vtkSmartPointer<vtkXMLUnstructuredGridWriter> writer =
  //  vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
  //writer->SetFileName("UnstructuredGrid.vtu");
  //  writer->SetInputData(unstructuredGrid);
  //writer->Write();



//vtkXMLImageDataWriter *imgWriter = vtkXMLImageDataWriter::New();
//std::ostringstream oss;
// oss << prefix << "." << imgWriter->GetDefaultFileExtension();
// imgWriter->SetFileName( oss.str().c_str() );
// imgWriter->SetInputData( g );
// imgWriter->Write();
 //imgWriter->Delete();




// 
//vtkSmartPointer<vtkTensorGlyph> tensorGlyph_MIL1 = vtkSmartPointer<vtkTensorGlyph>::New();
//    tensorGlyph_MIL1->SetInputData(polyDataTen_MIL1);
//    tensorGlyph_MIL1->SetSourceConnection(arrowSource->GetOutputPort());
//	//tensorGlyph_MIL1->SetColorModeToScalars();
//	//tensorGlyph_MIL1->SetColorMode(1);
//	tensorGlyph_MIL1->SetColorModeToEigenvalues();
//    tensorGlyph_MIL1->ThreeGlyphsOff();
//	tensorGlyph_MIL1->SetScaleFactor(3);
//	tensorGlyph_MIL1->ColorGlyphsOn();
//    tensorGlyph_MIL1->ExtractEigenvaluesOn();
//	tensorGlyph_MIL1->GetExtractEigenvalues();
//	tensorGlyph_MIL1->Update();

//  // Create a lookup table to share between the mapper and the scalarbar
//
//vtkSmartPointer<vtkLookupTable> hueLut = vtkSmartPointer<vtkLookupTable>::New();
//	hueLut->SetRange(0.0,255.0);
//    hueLut->SetTableRange (0, 1);
//    hueLut->SetHueRange (0, 1);
//    hueLut->SetSaturationRange (1, 1);
//    hueLut->SetValueRange (1, 1);
//    hueLut->Build();
// 
//vtkSmartPointer<vtkPolyDataMapper> mapper1_MIL = vtkSmartPointer<vtkPolyDataMapper>::New();
//mapper1_MIL->SetInputData(glyph3D->GetOutput());
//    mapper1_MIL->Update();
//
////
vtkSmartPointer<vtkPolyDataWriter> polyDataVec_MIL_writer = vtkSmartPointer<vtkPolyDataWriter>::New();
polyDataVec_MIL_writer->SetFileName(argv[3]);
//polyDataVec_MIL_writer->SetInputData(polyDataTen_MIL);
polyDataVec_MIL_writer->SetInputConnection(glyph3D->GetOutputPort());
polyDataVec_MIL_writer->Write();

//////
//////vtkSmartPointer<vtkPolyDataMapper> mapper1_MIL1 = vtkSmartPointer<vtkPolyDataMapper>::New();
//////    mapper1_MIL1->SetInputConnection(tensorGlyph_MIL1->GetOutputPort());
//////    mapper1_MIL1->Update();
////
////
////	//vtkSmartPointer<vtkActor> actor1_MIL1 = vtkSmartPointer<vtkActor>::New();
//// //   actor1_MIL1->SetMapper(mapper1_MIL);
////	////actor1_MIL1->GetProperty()->SetColor(0,0,1);
////	//actor1_MIL1 ->GetProperty()->SetOpacity(1);
////	//	  
////	//renderer->AddActor(actor1_MIL1);
////
////	vtkSmartPointer<vtkActor> actor1_MIL2 = vtkSmartPointer<vtkActor>::New();
////    actor1_MIL2->SetMapper(mapper1_MIL);
////	//actor1_MIL2->GetProperty()->SetColor(1,0,0);
////	actor1_MIL2 ->GetProperty()->SetOpacity(1);
////
////	//vtkSmartPointer<vtkImageData> input1new=vtkSmartPointer<vtkImageData>::New();
////resample->Update();
//////
//////vtkImageData *input_new;
//////input_new->AllocateScalars(VTK_FLOAT, 3);
//////input_new->DeepCopy(resample->GetOutput());
////
////	
////	 //vtkSmartPointer<vtkImageCast> castFilter = 
////  //    vtkSmartPointer<vtkImageCast>::New();
////  // castFilter->SetOutputScalarTypeToUnsignedChar();
////  // castFilter->SetInputData(input_new);
////  // castFilter->Update();
//// 
////
////
////
////
////
////	renderer->AddActor(actor1_MIL2);
//
//    renderer->ResetCamera();
//	renderer->GetActiveCamera()->Zoom(2.2);
//	renWin->AddRenderer(renderer);
//    renWin->Render();
//
//
//    iren->Start();
//    iren->ReInitialize();
//    renderer->Delete();
//    renWin->Delete();
//    iren->Delete();
//  

return 0;

}
