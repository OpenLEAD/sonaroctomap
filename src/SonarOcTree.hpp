#ifndef SONAR_OCTREE_HPP
#define SONAR_OCTREE_HPP

#include "octomap/OcTree.h"
#include "matio.h"
#include <string>
#include <base/samples/SonarBeam.hpp>
#include <base/samples/RigidBodyState.hpp>

#include <stdexcept>
#include <sys/stat.h>
#include <boost/concept_check.hpp>
#include <octomap/math/Utils.h>

namespace octomap {

class SonarOcTree;
  
typedef void (SonarOcTree::*pUpdateMethod)(const octomap::OcTreeKey&,const double& );
  

/** 
 *Extension of the Octree class with methods to use sonars
 *
 */

class SonarOcTree: public OcTree {
  
	double rbinlimits[500];
	matvar_t *matvar;
	mat_t *openmatfp;

private:
	//variables to use in the evaluateSonarBeam method
        double sum_gain_occ;
	double sum_gain2;
	
	double bin_max_occ;

	
	const float * logitprob;

public:
	void updater(const octomap::OcTreeKey& stepkey, const double&  gain);

	/// Default constructor, sets resolution of leafs
	SonarOcTree(double resolution,std::string filename = "ResizeRR.mat", std::string varname = "ResizeRR") :
			OcTree(resolution) {
	for(int i=0; i<=373; i++)
	  rbinlimits[i]=(0.2*i);
	
	struct stat buffer; 
	if(stat (filename.c_str(), &buffer))
	  throw std::runtime_error( "No "+filename+" file detected in the current folder." );
	
	
	openmatfp = Mat_Open(filename.c_str(),MAT_ACC_RDONLY);
	matvar = Mat_VarRead(openmatfp,varname.c_str());
	this->setClampingThresMax(0.9999);
	this->setClampingThresMin(0.0001);
	}
	;

	/**
	 * Reads an OcTree from a binary file
	 * @param _filename
	 *
	 */
	SonarOcTree(std::string _filename) :
			OcTree(_filename) {
	}
	;

	~SonarOcTree() {
	Mat_VarFree(matvar);
	Mat_Close(openmatfp);
	}
	;

	/// virtual constructor: creates a new object of same type
	/// (Covariant return type requires an up-to-date compiler)
	SonarOcTree* create() const {
		return new SonarOcTree(resolution,"ResizeRR.mat","ResizeRR");
	}

	std::string getTreeType() const {
		return "OcTree";
	}

	bool updateOccupancyNodeBinRay(const point3d& origin, const point3d& end,
			float log_odd_update, bool lazy_eval);

	/*	bool insertBinRay(const point3d& origin, const point3d& mid,
	 const point3d& end, double maxrange, float log_odd_update,
	 bool lazy_eval);*/

	bool createBin( int bin, double bearing, float poweroffset, base::samples::RigidBodyState sonar_state );
        bool createBin(int bin, double bearing = 0, float offset=0, double alpha = 0, double beta = 0);


	bool binShape( int bin, double bearing,   pUpdateMethod fnode, Eigen::Affine3d sonar_orientation );
	

	bool insertBinsRay(std::vector<uint8_t> beam_vector,
			octomath::Vector3 origin, octomath::Vector3 ray_direction,
			double length);

	bool insertOccupancyBin(std::vector<uint8_t> beam_vector,
			octomath::Vector3 origin, octomath::Vector3 ray_direction,
			double length);

	bool insertRealBeam(const base::samples::SonarBeam& beam,
			base::samples::RigidBodyState& sonar_state);
	
	bool insertBeam(const base::samples::SonarBeam& beam,
			base::samples::RigidBodyState sonar_state);
	
	bool mergeTrees(octomap::SonarOcTree &other,octomap::point3d offset);
	
	double compareTrees(const SonarOcTree& tree, base::Vector3d& tree_position, base::Quaterniond& tree_orientation);
	
	double evaluateSonarBeam( const Eigen::Affine3d& sonar_pose, const base::samples::SonarBeam& sonar_beam);
	
	double evaluateBinBinary(int current_bin, const base::samples::SonarBeam& sonar_beam,const Eigen::Affine3d& sonar_pose);
	
	double evaluateBinLQ(int current_bin, const base::samples::SonarBeam& sonar_beam, Eigen::Affine3d& sonar_pose);
	
	void projectIntensityBinary(const octomap::OcTreeKey& key,const double& alpha);
	
	void projectIntesityLQ(const octomap::OcTreeKey& key,const double& alpha);
	
	

};



} //end namespace

#endif
