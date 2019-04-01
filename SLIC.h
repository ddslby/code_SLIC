// SLIC.h: interface for the SLIC class.

#if !defined(_SLIC_H_INCLUDED_)
#define _SLIC_H_INCLUDED_

#include <vector>
#include <string>
#include <algorithm>
using namespace std;

class SLIC
{
public:
	SLIC();
	virtual ~SLIC();
	//============================================================================
	// Superpixel segmentation for a given step size (superpixel size ~= step*step)
	//给定步长的超像素分割（超像素大小~=step*step） 
	//============================================================================
	void DoSuperpixelSegmentation_ForGivenSuperpixelSize(
		const unsigned int*                            ubuff,//Each 32 bit unsigned int contains ARGB pixel values.(每个32位无符号int都包含argb像素值)
		const int					width,
		const int					height,
		int*&						klabels,
		int&						numlabels,
		const int&					superpixelsize,
		const double&                                   compactness);
	//============================================================================
	// Superpixel segmentation for a given number of superpixels
	//给定数量的超像素的超像素分割 
	//============================================================================
	void DoSuperpixelSegmentation_ForGivenNumberOfSuperpixels(
		const unsigned int*                             ubuff,
		const int					width,
		const int					height,
		int*&						klabels,
		int&						numlabels,
		const int&					K,//required number of superpixels
		const double&                                   compactness);//10-20 is a good value for CIELAB space
	//============================================================================
	// Supervoxel segmentation for a given step size (supervoxel size ~= step*step*step)
	//给定步长的超体素分割（超体素大小~=step*step*step） 
	//============================================================================
	void DoSupervoxelSegmentation(
		unsigned int**&		ubuffvec,
		const int&					width,
		const int&					height,
		const int&					depth,
		int**&						klabels,
		int&						numlabels,
		const int&					supervoxelsize,
		const double&                                   compactness);
	//============================================================================
	// Save superpixel labels in a text file in raster scan order
	//以光栅扫描顺序将超像素标签保存在文本文件中
	//============================================================================
	void SaveSuperpixelLabels(
		const int*&					labels,
		const int&					width,
		const int&					height,
		const string&				filename,
		const string&				path);
	//============================================================================
	// Save supervoxel labels in a text file in raster scan, depth order
	//在光栅扫描的文本文件中保存超体素标签，深度顺序
	//============================================================================
	void SaveSupervoxelLabels(
		const int**&				labels,
		const int&					width,
		const int&					height,
		const int&					depth,
		const string&				filename,
		const string&				path);
	//============================================================================
	// Function to draw boundaries around superpixels of a given 'color'.
	// Can also be used to draw boundaries around supervoxels, i.e layer by layer.
	//在给定“颜色”的超像素周围绘制边界的函数也可用于在超像素周围绘制边界，即逐层绘制。 
	//============================================================================
	void DrawContoursAroundSegments(
		unsigned int*&				segmentedImage,
		int*&						labels,
		const int&					width,
		const int&					height,
		const unsigned int&			color);

private:
	//============================================================================
	// The main SLIC algorithm for generating superpixels
	//超像素生成的主要SLIC算法
	//============================================================================
	void PerformSuperpixelSLIC(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		int*&						klabels,
		const int&					STEP,
		const vector<double>&		edgemag,
		const double&				m = 10.0);
	//============================================================================
	// The main SLIC algorithm for generating supervoxels
	//产生超体素的主要SLIC算法
	//============================================================================
	void PerformSupervoxelSLIC(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		vector<double>&				kseedsz,
		int**&						klabels,
		const int&					STEP,
		const double&				compactness);
	//============================================================================
	// Pick seeds for superpixels when step size of superpixels is given.
	//当给定超级像素的步长时，选择超级像素的种子。 
	//============================================================================
	void GetLABXYSeeds_ForGivenStepSize(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const int&					STEP,
		const bool&					perturbseeds,
		const vector<double>&		edgemag);
	//============================================================================
	// Pick seeds for supervoxels
	//为超体素选择种子
	//============================================================================
	void GetKValues_LABXYZ(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		vector<double>&				kseedsz,
		const int&					STEP);
	//============================================================================
	// Move the superpixel seeds to low gradient positions to avoid putting seeds
	// at region boundaries.
	//将超像素种子移动到低渐变位置以避免将种子放置在区域边界。
	//============================================================================
	void PerturbSeeds(
		vector<double>&				kseedsl,
		vector<double>&				kseedsa,
		vector<double>&				kseedsb,
		vector<double>&				kseedsx,
		vector<double>&				kseedsy,
		const vector<double>&		edges);
	//============================================================================
	// Detect color edges, to help PerturbSeeds() 彩色边缘检测
	//============================================================================
	void DetectLabEdges(
		const double*				lvec,
		const double*				avec,
		const double*				bvec,
		const int&					width,
		const int&					height,
		vector<double>&				edges);
	//============================================================================
	// sRGB to XYZ conversion; helper for RGB2LAB()
	//============================================================================
	void RGB2XYZ(
		const int&					sR,
		const int&					sG,
		const int&					sB,
		double&						X,
		double&						Y,
		double&						Z);
	//============================================================================
	// sRGB to CIELAB conversion (uses RGB2XYZ function)
	//============================================================================
	void RGB2LAB(
		const int&					sR,
		const int&					sG,
		const int&					sB,
		double&						lval,
		double&						aval,
		double&						bval);
	//============================================================================
	// sRGB to CIELAB conversion for 2-D images
	//============================================================================
	void DoRGBtoLABConversion(
		const unsigned int*&		ubuff,
		double*&					lvec,
		double*&					avec,
		double*&					bvec);
	//============================================================================
	// sRGB to CIELAB conversion for 3-D volumes
	//============================================================================
	void DoRGBtoLABConversion(
		unsigned int**&				ubuff,
		double**&					lvec,
		double**&					avec,
		double**&					bvec);
	//============================================================================
	// Post-processing of SLIC segmentation, to avoid stray labels.
	//============================================================================
	void EnforceLabelConnectivity(
		const int*					labels,
		const int					width,
		const int					height,
		int*&						nlabels,//input labels that need to be corrected to remove stray labels
		int&						numlabels,//the number of labels changes in the end if segments are removed
		const int&					K); //the number of superpixels desired by the user
	//============================================================================
	// Post-processing of SLIC supervoxel segmentation, to avoid stray labels.
	//============================================================================
	void EnforceSupervoxelLabelConnectivity(
		int**&						labels,//input - previous labels, output - new labels
		const int&					width,
		const int&					height,
		const int&					depth,
		int&						numlabels,
		const int&					STEP);

private:
	int										m_width;
	int										m_height;
	int										m_depth;

	double*									m_lvec;
	double*									m_avec;
	double*									m_bvec;

	double**								m_lvecvec;
	double**								m_avecvec;
	double**								m_bvecvec;
};

#endif // !defined(_SLIC_H_INCLUDED_)
