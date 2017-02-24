#ifndef SLIC_H
#define SLIC_H

/* slic.h.
 *
 * Written by: Pascal Mettes.
 *
 * This file contains the class elements of the class Slic. This class is an
 * implementation of the SLIC Superpixel algorithm by Achanta et al. [PAMI'12,
 * vol. 34, num. 11, pp. 2274-2282].
 *
 * This implementation is created for the specific purpose of creating
 * over-segmentations in an OpenCV-based environment.
 */

#include <opencv2\opencv.hpp>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <float.h>
#include <utility>

using namespace std;
using namespace cv;

/* 2d matrices are handled by 2d vectors. */
#define vec2dd vector<vector<double> >
#define vec2di vector<vector<int> >
#define vec2db vector<vector<bool> >
/* The number of iterations run by the clustering algorithm. */
#define NR_ITERATIONS 10

/*
 * class Slic.
 *
 * In this class, an over-segmentation is created of an image, provided by the
 * step-size (distance between initial cluster locations) and the colour
 * distance parameter.
 */
class Slic {
    private:
        /* The cluster assignments and distance values for each pixel. */
        vec2di clusters;
        vec2dd distances;
        
        /* The LAB and xy values of the centers. */
        vec2dd centers;
        /* The number of occurences of each center. */
        vector<int> center_counts;
        
        /* The step size per cluster, and the colour (nc) and distance (ns)
         * parameters. */
		int step = 0, nc = 0, ns = 0;
        
        /* Compute the distance between a center and an individual pixel. */
        double compute_dist(int ci, Point pixel, Scalar colour);
        /* Find the pixel with the lowest gradient in a 3x3 surrounding. */
        Point find_local_minimum( Point center);
        
        /* Remove and initialize the 2d vectors. */
        void clear_data();
        void init_data();

		/* Replacement for cvGet2D function*/
		Scalar get_Scal(int row, int col);
		void set_Scal(Mat &img, int row, int col, Vec3d colour);

		Mat img;
		
		int neutral = 1e4;
    public:
		
        /* Class constructors and deconstructors. */
		Slic(Mat image, int nc, int step, int ns);
        ~Slic();
        
        /* Generate an over-segmentation for an image. */
        void generate_superpixels();
        /* Enforce connectivity for an image. */
        void create_connectivity();
        
        /* Draw functions. Resp. displayal of the centers and the contours. */
        void display_center_grid(Mat &image, Vec3d colour);
        void display_contours(Mat &image, Vec3d colour);
		vector<cv::Vec3d> colour_with_cluster_means();

		/*Generate pixel assignment Mat file*/
		cv::Mat pixel_assignments();

		/*Generate a list of neighbouring superpixels*/
		//std::vector< std::pair<int, int> > SPneighbours(cv::Mat pxl_wise_assmts);
		std::vector<Vec3d> SPneighbours(cv::Mat pxl_wise_assmts);
		

		/* Get number of pixels */
		int get_nSP();

		/* Set parameters for superpixelation*/
		void set_param(int nc, int step, int ns);
};

#endif
