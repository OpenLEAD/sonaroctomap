#include <stdio.h>
#include <octomap/math/Utils.h>
#include <octomap/octomap.h>
#include <base/Time.hpp>
#include <sonaroctomap/SonarOcTree.hpp>

#include <octomap_wrapper/OctomapWrapper.hpp>
#include <octomap_wrapper/Conversion.hpp>

void print_query_info(octomath::Vector3 query, octomap::OcTreeNode* node) {
	if (node != NULL) {
		std::cout << "occupancy probability at " << query << ":\t "
				<< node->getOccupancy() << std::endl;
	} else
		std::cout << "occupancy probability at " << query << ":\t is unknown"
				<< std::endl;
}

void testWrappers (octomap::SonarOcTree &sonarTree){
		
		octomap_wrapper::OctomapWrapper* wrapper =
			new octomap_wrapper::OctomapWrapper();
			
		sonarTree.write("sonar.ot");	

		octomap::OcTree* tree = dynamic_cast<octomap::OcTree*>(&sonarTree);
		
		octomap_wrapper::fullMapToMsg(*tree, *wrapper);

		octomap::AbstractOcTree* abstree = octomap_wrapper::fullMsgToMap(*wrapper);
		
		if (abstree){
                    octomap::OcTree* new_tree = dynamic_cast<octomap::OcTree*>(abstree);
		    new_tree->write("new_tree.ot");
		    
		     for (octomap::OcTree::tree_iterator it = new_tree->begin_tree(0), end =
			tree->end_tree(); it != end; ++it) {

		        if (it.isLeaf()) {
		           std::cout << "occupancy = " << it->getOccupancy()<<  std::endl;
		        }
		    }
		        
                }
}

void mergeTrees(octomap::OcTree &rootTree,octomap::OcTree &localTree,octomap::point3d offset){
  
  rootTree.expand();
  localTree.expand();
  //verify the consequences of diferent depths
  unsigned int depthRoot = rootTree.getTreeDepth();
  unsigned int depthLocal = localTree.getTreeDepth();
  
  for (octomap::OcTree::leaf_iterator it = localTree.begin(depthLocal), end = localTree.end(); it != end; ++it)
  {
        
	octomap::point3d localPoint = it.getCoordinate();
	
	octomap::point3d rootPoint = offset + localPoint; 
	
	octomap::OcTreeNode* rootNode = rootTree.search(rootPoint);
	
        // check occupancy prob:
        double localProb = it->getOccupancy();
	double rootProb; 
	
	if(!rootNode){
	  rootProb = -1;
	  std::cout << "point only in localTree with prob = "<< localProb << std::endl;
	  rootTree.updateNode(rootPoint.x(),rootPoint.y(),rootPoint.z(),octomap::logodds(localProb),true);
	}
	else
	{
	  rootProb = rootNode->getOccupancy();
	  std::cout << "occupancy difference = "<< rootProb << " - " << localProb << " = " << rootProb - localProb << std::endl;
	  double diffProb = rootProb - localProb;
	  rootTree.updateNode(rootPoint.x(),rootPoint.y(),rootPoint.z(),octomap::logodds(diffProb),true);
	}
	delete rootNode;
  }  
}

bool beamTreeCreator(octomap::SonarOcTree &tree)
{
    base::samples::SonarBeam beam;
    beam.time = base::Time::now();
    beam.bearing.rad = 0.0;
    beam.sampling_interval = 1.0;
    beam.speed_of_sound = 1.0;
    beam.beamwidth_vertical = 0.610865238;
    beam.beamwidth_horizontal = 0.052359878;

    std::vector < uint8_t > bins(10);
    for (int i = 0; i < 10; ++i) 
    {
        bins[i] = 9.3*i;

    }
    beam.beam = bins;

    base::samples::RigidBodyState sonar_state;
    return tree.insertBeam(beam, sonar_state);
	
}

void cubeTreeCreator(octomap::SonarOcTree &tree,octomap::point3d origin, double prob, int size)
{
    for (int i = -size; i < size; i++) {
        for (int j = -size; j < size; j++) {
	    for (int k = -size; k < size; k++) {
		float x = i*0.049999;
		float y = j*0.049999;
		float z = k*0.049999;
	        octomap::point3d currentPoint(x,y,z);
		octomap::point3d endpoint = currentPoint + origin;
		float logOdd = octomap::logodds(prob);
		tree.updateNode(endpoint,logOdd,true); 
	    }
	}
    }
}

int main(int argc, char** argv) {
        
	octomap::SonarOcTree* sonarCube1 = new octomap::SonarOcTree(0.05);
	octomap::point3d origin1(-0.01,0.01,0.01);
	cubeTreeCreator(*sonarCube1,origin1,0.9,10);
	octomap::OcTree* cube1 = dynamic_cast<octomap::OcTree*>(sonarCube1);
		
	octomap::SonarOcTree* sonarCube2 = new octomap::SonarOcTree(0.05);
	octomap::point3d origin2(1.0,0.0,0.0);
	cubeTreeCreator(*sonarCube2,origin2,0.1,10);
	octomap::OcTree* cube2 = dynamic_cast<octomap::OcTree*>(sonarCube2);
	
	octomap::SonarOcTree* sonarCube3 = new octomap::SonarOcTree(0.05);
	octomap::point3d origin3(0.0,1.0,0.0);
	cubeTreeCreator(*sonarCube3,origin3,0.9,10);
	octomap::OcTree* cube3 = dynamic_cast<octomap::OcTree*>(sonarCube3);
	
	octomap::SonarOcTree* sonarCube4 = new octomap::SonarOcTree(0.05);
	octomap::point3d origin4(0.0,1.5,0.0);
	cubeTreeCreator(*sonarCube4,origin4,0.1,10);
	octomap::OcTree* cube4 = dynamic_cast<octomap::OcTree*>(sonarCube4);
	
	octomap::point3d zero(0.0,0.0,0.0);
	mergeTrees(*cube1,*cube2,zero);
	mergeTrees(*cube1,*cube3,zero);
	mergeTrees(*cube1,*cube4,zero);

	std::cout << "performing some queries:" << std::endl;
  
        octomap::point3d query (0., 0., 0.);
        octomap::OcTreeNode* result = cube1->search (query);
        print_query_info(query, result);
	
	query = octomap::point3d(1.,0.,0.);
        result = cube1->search (query);
        print_query_info(query, result);

        query = octomap::point3d(0.,.6,0.);
        result = cube1->search (query);
        print_query_info(query, result);

        query = octomap::point3d(0.,1.1,0.);
        result = cube1->search (query);
        print_query_info(query, result);
	
	query = octomap::point3d(0.,1.8,0.);
        result = cube1->search (query);
        print_query_info(query, result);
	
        query = octomap::point3d(0.,10,0.);
        result = cube1->search (query);
        print_query_info(query, result);


        std::cout << std::endl;
	
	
	
	
	cube1->write("tree.ot");
	
 	return 0;
}
		
		
		
