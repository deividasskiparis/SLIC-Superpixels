//#define _CRT_SECURE_NO_WARNINGS

#include "slic.h"

/*
 * Constructor. Nothing is done here.
 */
Slic::Slic(Mat image, int nc, int step, int ns) {
	this->step = step;
	this->nc = nc;
	this->ns = ns;
	img = image.clone();
	init_data();
}

/*
 * Destructor. Clear any present data.
 */
Slic::~Slic() {
    clear_data();
}

/*
 * Clear the data as saved by the algorithm.
 *
 * Input : -
 * Output: -
 */
void Slic::clear_data() {
    //clusters.clear();
    //distances.clear();
	//centers.clear();
	//center_counts.clear();
}

/*
 * Initialize the cluster centers and initial values of the pixel-wise cluster
 * assignment and distance values.
 *
 * Input : The image (IplImage*).
 * Output: -
 */
void Slic::init_data() {
	/* Initialize the cluster and distance matrices. */
	vector<int> cr(img.rows, neutral);
	vector<double> dr(img.rows, FLT_MAX);
	clusters = vec2di(img.cols, cr);
	distances = vec2dd(img.cols, dr);

    
    /* Initialize the centers and counters. */
	for (int i = step; i < img.cols - step / 2; i += step) {
		for (int j = step; j < img.rows - step / 2; j += step) {
            vector<double> center;
            /* Find the local minimum (gradient-wise). */
			Point nc = find_local_minimum(Point(i, j));
			Scalar colour = get_Scal(nc.y, nc.x);
            
            /* Generate the center vector. */
            center.push_back(colour.val[0]);
            center.push_back(colour.val[1]);
            center.push_back(colour.val[2]);
            center.push_back(nc.x);
            center.push_back(nc.y);
            
            /* Append to vector of centers. */
            centers.push_back(center);
            center_counts.push_back(0);
        }
    }

}

/*
 * Compute the distance between a cluster center and an individual pixel.
 *
 * Input : The cluster index (int), the pixel (CvPoint), and the Lab values of
 *         the pixel (CvScalar).
 * Output: The distance (double).
 */
double Slic::compute_dist(int ci, Point pixel, Scalar colour) {


    double dc = sqrt(pow(centers[ci][0] - colour.val[0], 2) + pow(centers[ci][1]
            - colour.val[1], 2) + pow(centers[ci][2] - colour.val[2], 2));
    double ds = sqrt(pow(centers[ci][3] - pixel.x, 2) + pow(centers[ci][4] - pixel.y, 2));
    
    return sqrt(pow(dc / nc, 2) + pow(ds / ns, 2));
    
    //double w = 1.0 / (pow(ns / nc, 2));
    //return sqrt(dc) + sqrt(ds * w);

}

/*
 * Find a local gradient minimum of a pixel in a 3x3 neighbourhood. This
 * method is called upon initialization of the cluster centers.
 *
 * Input : The image (IplImage*) and the pixel center (CvPoint).
 * Output: The local gradient minimum (CvPoint).
 */
Point Slic::find_local_minimum(Point center) {

    double min_grad = FLT_MAX;
    Point loc_min = Point(center.x, center.y);
    
    for (int i = center.x-1; i < center.x+2; i++) {
        for (int j = center.y-1; j < center.y+2; j++) {
            Scalar c1 = get_Scal(j+1, i);
			Scalar c2 = get_Scal(j, i + 1);
			Scalar c3 = get_Scal(j, i);
            /* Convert colour values to grayscale values. */
            double i1 = c1.val[0];
            double i2 = c2.val[0];
            double i3 = c3.val[0];
            /*double i1 = c1.val[0] * 0.11 + c1.val[1] * 0.59 + c1.val[2] * 0.3;
            double i2 = c2.val[0] * 0.11 + c2.val[1] * 0.59 + c2.val[2] * 0.3;
            double i3 = c3.val[0] * 0.11 + c3.val[1] * 0.59 + c3.val[2] * 0.3;*/
            
            /* Compute horizontal and vertical gradients and keep track of the
               minimum. */
            if (sqrt(pow(i1 - i3, 2)) + sqrt(pow(i2 - i3,2)) < min_grad) {
                min_grad = fabs(i1 - i3) + fabs(i2 - i3);
                loc_min.x = i;
                loc_min.y = j;
            }
        }
    }
    return loc_min;
}

/*
 * Compute the over-segmentation based on the step-size and relative weighting
 * of the pixel and colour values.
 *
 * Input : The Lab image (IplImage*), the stepsize (int), and the weight (int).
 * Output: -
 */
void Slic::generate_superpixels() {


    
    /* Clear previous data (if any), and re-initialize it. */
    //clear_data();
    //init_data(image);
    
    /* Run EM for 10 iterations (as prescribed by the algorithm). */
    for (int i = 0; i < NR_ITERATIONS; i++) {
        /* Reset distance values. */
		vector<double> dist(img.rows, FLT_MAX);
		distances = vec2dd(img.cols, dist);
	
        for (int j = 0; j < (int) centers.size(); j++) {
            /* Only compare to pixels in a 2 x step by 2 x step region. */
            for (int k = centers[j][3] - step; k < centers[j][3] + step; k++) {
                for (int l = centers[j][4] - step; l < centers[j][4] + step; l++) {
                
					if (k >= 0 && k < img.cols && l >= 0 && l < img.rows) {
						Scalar colour = get_Scal(l, k);
                        double d = compute_dist(j, Point(k,l), colour);
                        
                        /* Update cluster allocation if the cluster minimizes the
                           distance. */
						//cout << "d: " << d << "\tk: " << k << "\tl: " << l << endl;
	
						if (d < distances[k][l]) {
							distances[k][l] = d;
							clusters[k][l] = j;
						}
                    }
                }
            }
        }
        
        /* Clear the center values. */
        for (int j = 0; j < (int) centers.size(); j++) {
            centers[j][0] = centers[j][1] = centers[j][2] = centers[j][3] = centers[j][4] = 0;
            center_counts[j] = 0;
        }
        
        /* Compute the new cluster centers. */
		for (int j = 0; j < img.cols; j++) {
			for (int k = 0; k < img.rows; k++) {
                int c_id = clusters[j][k];
                
				if (c_id != neutral) {
					Scalar colour = get_Scal(k, j);
                    
                    centers[c_id][0] += colour.val[0];
                    centers[c_id][1] += colour.val[1];
                    centers[c_id][2] += colour.val[2];
                    centers[c_id][3] += j;
                    centers[c_id][4] += k;
                    
                    center_counts[c_id] += 1;
                }
            }
        }

        /* Normalize the clusters. */
        for (int j = 0; j < (int) centers.size(); j++) {
            centers[j][0] /= center_counts[j];
            centers[j][1] /= center_counts[j];
            centers[j][2] /= center_counts[j];
            centers[j][3] /= center_counts[j];
            centers[j][4] /= center_counts[j];
        }
    }

}

/*
 * Enforce connectivity of the superpixels. This part is not actively discussed
 * in the paper, but forms an active part of the implementation of the authors
 * of the paper.
 *
 * Input : The image (IplImage*).
 * Output: -
 */
void Slic::create_connectivity() {

    int label = 0, adjlabel = 0;
	const int lims = (img.cols * img.rows) / ((int)centers.size());
    
    const int dx4[4] = {-1,  0,  1,  0};
	const int dy4[4] = { 0, -1,  0,  1};
    
    /* Initialize the new cluster matrix. */
 //   vec2di new_clusters;
	//for (int i = 0; i < img.cols; i++) {
 //       vector<int> nc;
	//	for (int j = 0; j < img.rows; j++) {
	//		nc.push_back(neutral);
 //       }
 //       new_clusters.push_back(nc);
 //   }
	vector<int> nc(img.rows,neutral);
	vec2di new_clusters(img.cols, nc);

	for (int i = 0; i < img.cols; i++) {
		for (int j = 0; j < img.rows; j++) {
			if (new_clusters[i][j] == neutral) {
                vector<Point> elements;
                elements.push_back(Point(i, j));
            
                /* Find an adjacent label, for possible use later. */
                for (int k = 0; k < 4; k++) {
                    int x = elements[0].x + dx4[k], y = elements[0].y + dy4[k];
                    
					if (x >= 0 && x < img.cols && y >= 0 && y < img.rows) {
                        if (new_clusters[x][y] =! neutral) {
                            adjlabel = new_clusters[x][y];
                        }
                    }
                }
                
                int count = 1;
                for (int c = 0; c < count; c++) {
                    for (int k = 0; k < 4; k++) {
                        int x = elements[c].x + dx4[k], y = elements[c].y + dy4[k];
                        
						if (x >= 0 && x < img.cols && y >= 0 && y < img.rows) {
							if (new_clusters[x][y] == neutral && clusters[i][j] == clusters[x][y]) {
                                elements.push_back(Point(x, y));
                                new_clusters[x][y] = label;
                                count += 1;
                            }
                        }
                    }
                }
                
                /* Use the earlier found adjacent label if a segment size is
                   smaller than a limit. */
                if (count <= lims >> 2) {
                    for (int c = 0; c < count; c++) {
                        new_clusters[elements[c].x][elements[c].y] = adjlabel;
                    }
                    label -= 1;
                }
                label += 1;
            }
        }
    }

}

/*
 * Display the cluster centers.
 *
 * Input : The image to display upon (IplImage*) and the colour (CvScalar).
 * Output: -
 */
void Slic::display_center_grid(Mat &image, Vec3d colour) {
	Scalar col = Scalar(colour[0], colour[1], colour[2]);
    for (int i = 0; i < (int) centers.size(); i++) {
		circle(image, Point(centers[i][3], centers[i][4]), 2, col, 2);
    }

}

/*
 * Display a single pixel wide contour around the clusters.
 *
 * Input : The target image (IplImage*) and contour colour (CvScalar).
 * Output: -
 */
void Slic::display_contours(Mat &image, Vec3d colour) {

    const int dx8[8] = {-1, -1,  0,  1, 1, 1, 0, -1};
	const int dy8[8] = { 0, -1, -1, -1, 0, 1, 1,  1};
	
	/* Initialize the contour vector and the matrix detailing whether a pixel
	 * is already taken to be a contour. */
	vector<Point> contours;
	vector<bool> nb(image.rows, false);
	vec2db istaken(image.cols, nb);

    
    /* Go through all the pixels. */
    for (int i = 0; i < image.cols; i++) {
        for (int j = 0; j < image.rows; j++) {
            int nr_p = 0;
            
            /* Compare the pixel to its 8 neighbours. */
            for (int k = 0; k < 8; k++) {
                int x = i + dx8[k], y = j + dy8[k];
                
                if (x >= 0 && x < image.cols && y >= 0 && y < image.rows) {
                    if (istaken[x][y] == false && clusters[i][j] != clusters[x][y]) {
                        nr_p += 1;
                    }
                }
            }
            
            /* Add the pixel to the contour list if desired. */
            if (nr_p >= 2) {
                contours.push_back(Point(i,j));
                istaken[i][j] = true;
            }
        }
    }
    
    /* Draw the contour pixels. */
    for (int i = 0; i < (int)contours.size(); i++) {
        set_Scal(image, contours[i].y, contours[i].x, colour);
    }

}

/*
 * Generate matrix of pixel-wise cluster assignments
 *
 * Input - Target image
 * Output - Matrix of pixel-wise allocations to cluster centres
 *
 */
cv::Mat Slic::pixel_assignments(){

	cv::Mat_<int> px_assign = Mat(img.rows, img.cols, CV_8UC1, Scalar(0));
	for (int i = 0; i < px_assign.rows; ++i){
		for (int j = 0; j < px_assign.cols; ++j){
			//cout << "i: " << i << ", j: " << j << ", value: " << clusters.at(i).at(j) << endl;
			int clust_id = (int)clusters[j][i];
			px_assign(i, j) = clust_id;
		}
	}

	return(px_assign);
}

/*
 * Give the pixels of each cluster the same colour values. The specified colour
 * is the mean RGB colour per cluster.
 *
 * Input : The target image (Mat).
 * Output: -
 */
vector<cv::Vec3d> Slic::colour_with_cluster_means() {
	vector<cv::Vec3d> colours(centers.size(), Vec3d(0, 0, 0));
	vector<int> color_counts(centers.size(), 0);

    /* Gather the colour values per cluster. */
	for (int i = 0; i < img.cols; i++) {
		for (int j = 0; j < img.rows; j++) {
            int index = clusters[i][j];
			if (index != neutral){
				Scalar colour = get_Scal(j, i);
				colours[index].val[0] += colour.val[0];
				colours[index].val[1] += colour.val[1];
				colours[index].val[2] += colour.val[2];
				color_counts[index] += 1;
			}
        }
    }

    /* Divide by the number of pixels per cluster to get the mean colour. */
    for (int i = 0; i < (int)colours.size(); i++) {
		int cnt = color_counts[i];
		if (cnt == 0) continue;
		colours[i].val[0] /= cnt;
		colours[i].val[1] /= cnt;
		colours[i].val[2] /= cnt;
    }

    /* Fill in. */
	for (int i = 0; i < img.cols; i++) {
		for (int j = 0; j < img.rows; j++) {
			int cluster_idx = clusters[i][j];
			if (cluster_idx != neutral){
				Vec3d ncolour = colours[cluster_idx];
				set_Scal(img, j, i, ncolour);
			}

        }
    }

	return colours;
}
/* 
 * Generate a list of neighbouring superpixels
 * Makes a non-repeating list, meaning if 1 and 2 are neighbours the list will
 * contain either (1,2) or (2,1) but not both
 */
std::vector<Vec3d> Slic::SPneighbours(cv::Mat pxl_wise_assmts){

	std::vector<Vec3d> result;
	cv::Mat_<int> sps(pxl_wise_assmts);
	int w = pxl_wise_assmts.size().width, h = pxl_wise_assmts.size().height;
	
	int N_SP = get_nSP();
	cv::Mat_<int> assimts_temp = cv::Mat::zeros(cv::Size(N_SP,N_SP), CV_8UC1);


	for (int i = 0; i < h - 1; i++){ // for every row
		for (int j = 0; j < w - 1; j++){ // for every column
			//std::pair<int, int> p1(0, 0), p2(0, 0);
			Vec3d p1(0, 0, 1), p2(0, 0, 1);

			if ((i < h - 1) && (j < w - 1)){


				int p_ref = sps(i, j), p_dx = sps(i, j + 1), p_dy = sps(i + 1, j), p_dxy = sps(i + 1, j + 1);
				if (p_ref >= 0 && p_dx >= 0 && p_ref < N_SP && p_dx < N_SP)
					p_ref < p_dx ? p1[0] = p_ref, p1[1] = p_dx : p1[0] = p_dx, p1[1] = p_ref;
				if (p_ref >= 0 && p_dy >= 0 && p_ref < N_SP && p_dy < N_SP)
					p_ref < p_dx ? p2[0] = p_ref, p2[1] = p_dy : p2[0] = p_dy, p2[1] = p_ref;

					//p_ref < p_dy? p2 = std::make_pair(p_ref, p_dy) : p2 = std::make_pair(p_dy, p_ref);
			}

			if ((assimts_temp(p1[0], p1[1]) != 1) && (p1[0] != p1[1])){ // value has not been added to the list yet
				result.push_back(p1);
				assimts_temp(p1[0], p1[1]) = 1;
			}

			if ((assimts_temp(p2[0], p2[1]) != 1) && (p2[0] != p2[1])){ // value has not been added to the list yet
				result.push_back(p2);
				assimts_temp(p2[0], p2[1]) = 1;
			}
		}
		
	}
	return result;
}
/*
Get number of superpixels

*/
int Slic::get_nSP(){
	return (int)centers.size();
}

/*
Helper functions for opencv 1.x to opencv 2.x conversion
*/
Scalar Slic::get_Scal(int row, int col){
	Scalar v = img.at<Vec3b>(row, col);
	return v;
}
void Slic::set_Scal(Mat &in_img, int row, int col, Vec3d colour){

	Vec3b v(colour);
	in_img.at<Vec3b>(row, col) = v;

}
void Slic::set_param(int nc, int step, int ns){
	this->step = step;
	this->nc = nc;
	this->ns = ns;
}