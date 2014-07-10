#include "octomap/OcTree.h"
#include <string>
#include "base/samples/sonar_beam.h"
#include <base/samples/rigid_body_state.h>

namespace octomap {

/** 
 *Extension of the Octree class with methods to use sonars
 *
 */

class SonarOcTree: public OcTree {

public:

	/// Default constructor, sets resolution of leafs
	SonarOcTree(double resolution) :
			OcTree(resolution) {
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

	bool insertBinsRay(std::vector<uint8_t> beam_vector,
			octomath::Vector3 origin, octomath::Vector3 ray_direction,
			double length);

	bool insertBeam(const base::samples::SonarBeam& beam,
			base::samples::RigidBodyState sonar_state);

};

} //end namespace
