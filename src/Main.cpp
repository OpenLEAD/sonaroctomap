#include <stdio.h>
#include <boost/concept_check.hpp>
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
		float x = i*0.05;
		float y = j*0.05;
		float z = k*0.05;
	        octomap::point3d currentPoint(x,y,z);
		octomap::point3d endpoint = currentPoint + origin;
		float logOdd = octomap::logodds(prob);
		tree.updateNode(endpoint,logOdd,true); 
	    }
	}
    }
}

int main(int argc, char** argv) {
        
	octomap::SonarOcTree* sonarCube1 = new octomap::SonarOcTree(0.2);
	
	for(int i=0;i<150;i++){
	  std::cout << i << std::endl;
	  sonarCube1->CreateBin("ResizeRR.mat","ResizeRR",i);}
	  
// 	  sonarCube1->CreateBinPointCloud(0.2,"ResizeRR.mat","ResizeRR",50);
// 	  sonarCube1->CreateBinPointCloud(0.2,"ResizeRR.mat","ResizeRR",52);
// 	  
	  
// 	sonarCube1->CreateBinPointCloud(0.1,"ResizeRR.mat","ResizeRR",3,45*M_PI/180,45*M_PI/180);
// 	sonarCube1->CreateBinPointCloud(0.1,"ResizeRR.mat","ResizeRR",4,45*M_PI/180,45*M_PI/180,45*M_PI/180);
// 	sonarCube1->CreateBinPointCloud(0.1,"ResizeRRnorm.mat","ResizeRR",4,90*M_PI/180);
// 	sonarCube1->CreateBinPointCloud(0.1,"ResizeRRnorm.mat","ResizeRR",5);
// 	sonarCube1->CreateBinPointCloud(0.1,"ResizeRR.mat","ResizeRR",6);
	
// 	octomap::point3d origin1(0,0,0);
// 	cubeTreeCreator(*sonarCube1,origin1,0.9,10);
// 
// 		
// 	octomap::SonarOcTree* sonarCube2 = new octomap::SonarOcTree(0.1);
// 	octomap::point3d origin2(1.0,0.0,0.0);
// 	cubeTreeCreator(*sonarCube2,origin2,0.1,10);
// 	
// 	octomap::SonarOcTree* sonarCube3 = new octomap::SonarOcTree(0.1);
// 	octomap::point3d origin3(0.0,1.1,0.0);
// 	cubeTreeCreator(*sonarCube3,origin3,0.9,10);
// 	
// 	octomap::SonarOcTree* sonarCube4 = new octomap::SonarOcTree(0.1);
// 	octomap::point3d origin4(0.0,1.6,0.0);
// 	cubeTreeCreator(*sonarCube4,origin4,0.1,10);
// 	
// 	octomap::point3d zero(0.0,0.0,0.0);
// 	octomap::SonarOcTree* sonarCube0 = new octomap::SonarOcTree(0.1);
//         
// 	
// 	sonarCube1->mergeTrees(*sonarCube2,zero);
// 	sonarCube1->mergeTrees(*sonarCube3,zero);
// 	sonarCube1->mergeTrees(*sonarCube4,zero);	
// 
// 	std::cout << "performing some queries:" << std::endl;
//   
//         octomap::point3d query (0., 0., 0.);
//         octomap::OcTreeNode* result = sonarCube1->search (query);
//         print_query_info(query, result);
// 	
// 	query = octomap::point3d(1.,0.,0.);
//         result = sonarCube1->search (query);
//         print_query_info(query, result);
// 
//         query = octomap::point3d(0.,.6,0.);
//         result = sonarCube1->search (query);
//         print_query_info(query, result);
// 
//         query = octomap::point3d(0.,2.1,0.);
//         result = sonarCube1->search (query);
//         print_query_info(query, result);
// 	
// 	query = octomap::point3d(0.,2.8,0.);
//         result = sonarCube1->search (query);
//         print_query_info(query, result);
// 	
//         query = octomap::point3d(0.,10,0.);
//         result = sonarCube1->search (query);
//         print_query_info(query, result);


        std::cout << std::endl;
	
	sonarCube1->write("treeExtra.ot");
	
 	return 0;
}
		
		
		

