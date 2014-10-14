#include <stdio.h>
#include <boost/concept_check.hpp>
#include <octomap/math/Utils.h>
#include <octomap/octomap.h>
#include <base/Time.hpp>
#include <sonaroctomap/SonarOcTree.hpp>
#include <string>

#include <octomap_wrapper/OctomapWrapper.hpp>
#include <octomap_wrapper/Conversion.hpp>

#include <boost/random/mersenne_twister.hpp>
#include <boost/random/normal_distribution.hpp>
#include <boost/random/variate_generator.hpp>
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */


/* WORKAROUND */
#include <pthread.h>
void junk() {	
  int i;
  i=pthread_getconcurrency();
};
/* WORKAROUND */




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

void evaluateBeamTester()
{    octomap::SonarOcTree* tree = new octomap::SonarOcTree(0.2);
     octomap::OcTreeKey startkey = tree->coordToKey(-20,10,-20);
     octomap::OcTreeKey endkey = tree->coordToKey(20,11,20);	
     for(int stepX=startkey[0]; stepX<=endkey[0]; stepX++)
	{
               for(int stepY=startkey[1]; stepY<=endkey[1]; stepY++)
	       {
	           for(int stepZ=startkey[2]; stepZ<=endkey[2]; stepZ++)
		   {
		        octomap::OcTreeKey stepkey = octomap::OcTreeKey(stepX,stepY,stepZ);
			tree->updateNode(stepkey,octomap::logodds(0.9),false);
		   }
	       }
	}

     octomap::OcTreeKey startkey2 = tree->coordToKey(20,11,-20);
     octomap::OcTreeKey endkey2 = tree->coordToKey(21,30,20);	
     for(int stepX=startkey2[0]; stepX<=endkey2[0]; stepX++)
	{
               for(int stepY=startkey2[1]; stepY<=endkey2[1]; stepY++)
	       {
	           for(int stepZ=startkey2[2]; stepZ<=endkey2[2]; stepZ++)
		   {
		        octomap::OcTreeKey stepkey2 = octomap::OcTreeKey(stepX,stepY,stepZ);
			tree->updateNode(stepkey2,octomap::logodds(0.9),false);
		   }
	       }
	}	
	
    std::cout << "map created" << std::endl;
    
    base::samples::RigidBodyState rbs;
    rbs.initUnknown();
    base::samples::SonarBeam* sonar_beam =  new base::samples::SonarBeam();
    sonar_beam->bearing.rad = 0.785398163;
    sonar_beam->beam = std::vector<u_int8_t>(150);
    sonar_beam->sampling_interval = 1;
    sonar_beam->speed_of_sound = 0.2;

   std::vector<u_int8_t> projected_beam(150);
    
    tree->evaluateSonarBeam(rbs,*sonar_beam, projected_beam);
    
    sonar_beam->beam.swap(projected_beam);
    
    tree->insertRealBeam(*sonar_beam,rbs);
    
    std::cout<< "beam inserted" << std::endl;
    
    tree->write("cornerProj45.ot");
    
    std::cout << "tree saved" << std::endl;
    
    delete tree;
  
}


int main(int argc, char** argv) {
  
        evaluateBeamTester();
// 	typedef boost::mt19937                     ENG;    // Mersenne Twister
// 	typedef boost::normal_distribution<double> DIST;   // Normal Distribution
// 	typedef boost::variate_generator<ENG,DIST> GEN;    // Variate generator
//   
// 	ENG  eng;
// 	DIST dist(0,0.2);
// 	GEN  gen(eng,dist);
// 	
// 	std::srand (std::time(NULL));
// 
// 	
//         
// 	octomap::SonarOcTree* sonarCube1 = new octomap::SonarOcTree(0.2,"ResizeRR.mat","ResizeRR");
// 	base::samples::RigidBodyState sonar_state;
// 	base::Pose sonar_pose;
// 	sonar_pose.orientation.Identity();
// 	sonar_state.setPose(sonar_pose);
// 	
// 	int p;
// 	
// 	for(int voltas=0;voltas<4;voltas++){
// 	
// 	  struct timeval start, end;
// 	  gettimeofday(&start, NULL);
// 	  
// 	for(int i=voltas*10;i<10*(voltas+1);i++){
// 	  
//   
// 	  sonar_state.orientation = Eigen::AngleAxisd(rand()%360, Eigen::MatrixBase<base::Vector3d>::UnitZ()) * Eigen::AngleAxisd(base::Angle::deg2Rad(rand()%360), Eigen::MatrixBase<base::Vector3d>::UnitY());
// 	  
// 	  sonar_state.position = sonar_state.orientation.toRotationMatrix() *base::Position(10+(p=rand()%5)+gen(),0,0);
// 	  
// 	  //std::cout << i;	
// 	  
// 	  #pragma omp parallel num_threads(6)
// 	  {
// 	  #pragma omp for schedule(dynamic)
// 	  for(int bin=0;bin<(40+5*p);bin++)
// 	    sonarCube1->createBin(bin,base::Angle::deg2Rad(180),15,sonar_state);
// 	  }
// 	  sonarCube1->createBin(40+5*p,base::Angle::deg2Rad(180),32,sonar_state);
// 	  
// 	  
// 	  gettimeofday(&end, NULL);
// 	  unsigned long long useconds =end.tv_sec - start.tv_sec;
// 	  useconds *= 1000000;
// 	  useconds +=  end.tv_usec - start.tv_usec;
// 	  
// 	  std::cout<< (end.tv_sec - start.tv_sec) << "; "<< std::flush;
// 
// 	  //std::cout<< " - " << ((double) (i+1)*1000000)/useconds << " beams per second in " << end.tv_sec - start.tv_sec << "s" << std::endl;
// 	}
// 	std::stringstream ss;
// 	ss << voltas;
// 	  
// 	sonarCube1->write("treeExtraRandRPProgressiveNewMptm"+ss.str()+".ot"); //treeExtraRandRPProgressiveNewM
//         std::cout << "treeExtraRandRPProgressiveNewMptm"+ss.str()+".ot"<< " written" << std::endl;
// 	  
// 	}
	
		
		
		
		
		
		
/*		
	for(int voltas=0;voltas<10;voltas++){
	
	  struct timeval start, end;
	  gettimeofday(&start, NULL);
	  
	for(int i=voltas*360;i<360*(voltas+1);i++){
	  
  
	  
	  
	  sonar_state.position = 10*base::Position(cos(base::Angle::deg2Rad(i)),sin(base::Angle::deg2Rad(i)),0);
	  
	  //std::cout << i;	
	  
	  #pragma omp parallel num_threads(6)
	  {
	  #pragma omp for schedule(dynamic)
	  for(int bin=0;bin<40;bin++)
	    sonarCube1->CreateBin(bin,base::Angle::deg2Rad(180+i),10,sonar_state);
	  }
	  sonarCube1->CreateBin(40,base::Angle::deg2Rad(180+i),32,sonar_state);
	  
	  
	  gettimeofday(&end, NULL);
	  unsigned long long useconds =end.tv_sec - start.tv_sec;
	  useconds *= 1000000;
	  useconds +=  end.tv_usec - start.tv_usec;
	  
	  std::cout<< (end.tv_sec - start.tv_sec) << "; "<< std::flush;

	  //std::cout<< " - " << ((double) (i+1)*1000000)/useconds << " beams per second in " << end.tv_sec - start.tv_sec << "s" << std::endl;
	  
	}
	std::stringstream ss;
	ss << (voltas+5);
	  
	sonarCube1->write("treeExtraRand"+ss.str()+".ot");
        std::cout << "treeExtraRand"+ss.str()+".ot"<< " written" << std::endl;
	  
	}*/
	  











// 	for(int i=15;i<30;i++){
// 	  sonar_state.position = base::Position(0,0,(int)1);
// 	  std::cout << i << std::endl;
// 	  sonarCube1->CreateBin("ResizeRR.mat","ResizeRR",i,base::Angle:int:deg2Rad(i),32,sonar_state);}
	  
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

/*
        std::cout << std::endl;
	
	sonarCube1->write("treeExtra.ot");
	
 	return 0;*/
}
		
		
		

