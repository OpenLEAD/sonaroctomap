#ifndef SONAR_OCTREE_HPP
#define SONAR_OCTREE_HPP

#include "octomap/OcTree.h"
#include <string>
#include "base/samples/SonarBeam.hpp"
#include <base/samples/RigidBodyState.hpp>

namespace octomap {

/** 
 *Extension of the Octree class with methods to use sonars
 *
 */

class SonarOcTree: public OcTree {
  
	double rbinlimits[500];
public:

	/// Default constructor, sets resolution of leafs
	SonarOcTree(double resolution) :
			OcTree(resolution) {
	for(int i=0; i<=373; i++)
	  rbinlimits[i]=(0.2*i+0.3);
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
	}
	;

	/// virtual constructor: creates a new object of same type
	/// (Covariant return type requires an up-to-date compiler)
	SonarOcTree* create() const {
		return new SonarOcTree(resolution);
	}

	std::string getTreeType() const {
		return "OcTree";
	}

	bool updateOccupancyNodeBinRay(const point3d& origin, const point3d& end,
			float log_odd_update, bool lazy_eval);

	/*	bool insertBinRay(const point3d& origin, const point3d& mid,
	 const point3d& end, double maxrange, float log_odd_update,
	 bool lazy_eval);*/

	bool CreateBin(std::string filename, std::string varname, int bin, 
		       double bearing, float offset, base::Matrix3d& Pantilt );
	
	bool CreateBin(std::string filename, std::string varname, int bin, 
		       double bearing = 0, float offset=0, double alpha = 0, double beta = 0);

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
	
	bool mergeTrees(octomap::SonarOcTree &tree2,octomap::point3d offset);
	
	double compareTrees(const SonarOcTree& tree, base::Vector3d& treePosition, base::Quaterniond& treeOrientation);

};



} //end namespace

#endif
