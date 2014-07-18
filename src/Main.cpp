#include <stdio.h>
#include <octomap/math/Utils.h>
#include <base/Time.hpp>
#include <sonaroctomap/SonarOcTree.hpp>

#include <octomap_wrapper/OctomapWrapper.hpp>
#include <octomap_wrapper/conversion.hpp>

void print_query_info(octomath::Vector3 query, octomap::OcTreeNode* node) {
	if (node != NULL) {
		std::cout << "occupancy probability at " << query << ":\t "
				<< node->getOccupancy() << std::endl;
	} else
		std::cout << "occupancy probability at " << query << ":\t is unknown"
				<< std::endl;
}

int main(int argc, char** argv) {

	base::samples::SonarBeam beam;

	beam.time = base::Time::now();

	beam.bearing.rad = 0.0;

	beam.sampling_interval = 1.0;

	beam.speed_of_sound = 1.0;

	beam.beamwidth_vertical = 0.610865238;

	beam.beamwidth_horizontal = 0.052359878;

	std::vector < uint8_t > bins(10);

	for (int i = 0; i < 10; ++i) {
		bins[i] = i;

	}

	beam.beam = bins;

	octomap::SonarOcTree sonarTree(0.05);

	octomap_wrapper::OctomapWrapper* wrapper =
			new octomap_wrapper::OctomapWrapper();

	octomath::Vector3 origin(0.0, 0.0, 0.0);

	base::samples::RigidBodyState sonar_state;
	bool ok = sonarTree.insertBeam(beam, sonar_state);

	if (ok) {

		octomap_wrapper::binaryMapToMsg < octomap::OcTree
				> ((octomap::OcTree)sonarTree, *wrapper);

		std::cout << "tree - wrapper" << std::endl;

		octomap::OcTree* new_tree = octomap_wrapper::binaryMsgToMap(*wrapper);

		std::cout << "wrapper - tree" << std::endl;

		new_tree->writeBinary("new_tree.bt");

		sonarTree.writeBinary("teste.bt");

	} else

		return 0;
}

