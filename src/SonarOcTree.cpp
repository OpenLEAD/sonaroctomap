#include <sonaroctomap/SonarOcTree.hpp>
#include "matio.h"
#include <octomap/math/Utils.h>

using namespace octomap;

bool SonarOcTree::updateOccupancyNodeBinRay(const point3d& origin,
		const point3d& end, float log_odd_update, bool lazy_eval) {

	if (!this->computeRayKeys(origin, end, this->keyrays.at(0))) {
		return false;
	}

	for (KeyRay::iterator it = this->keyrays[0].begin();
			it != this->keyrays[0].end(); it++) {
		updateNode(*it, log_odd_update, lazy_eval); // insert the log odd of the entire BinRay
	}

	return true;
}

const double ThetaMin = -M_PI*686/1800;
const double ThetaMax = M_PI*874/1800;
const double PhiMin = -M_PI*315/1800;
const double PhiMax = M_PI*290/1800;

bool SonarOcTree::CreateBinPointCloud(double octo_resolution, std::string filename, std::string varname, double length){
  mat_t *openmatfp;
  matvar_t *matvar;
  openmatfp = Mat_Open(filename.c_str(),MAT_ACC_RDONLY);
  matvar = Mat_VarRead(openmatfp,varname.c_str());
  octomap::Pointcloud bin;
  double rhomin = 8;
  double rhomax = 12;
  for(double rho = DBL_EPSILON+rhomin; rho<rhomax;rho+=octo_resolution){
    
    
  }
    
  std::cout << ((double*)matvar->data)[0];
    
  //std::cout<< "MEUS TESTE -> " << ((double*)matvar->data)[matvar->dims[0]*2+1] << std::endl;
  Mat_VarPrint(matvar,1);
  Mat_VarFree(matvar);
  Mat_Close(openmatfp);
  return true;
}

bool SonarOcTree::insertBinsRay(std::vector<uint8_t> beam_vector,
		octomath::Vector3 origin, octomath::Vector3 ray_direction,
		double length) {

	octomath::Vector3 current_origin = origin;
	octomath::Vector3 current_endpoint;

	for (int i = 0; i < beam_vector.size(); i++) {
		int multiplier = i + 1;
		current_endpoint = ray_direction * multiplier * length;
		int bin_value = beam_vector[i - 1];
		
		//TODO conversion between the value of the bins and a log_odd, for now /100

		double prob = bin_value / 100.0;

		if(prob>0.1){

		if (!this->updateOccupancyNodeBinRay(current_origin, current_endpoint,
				 octomap::logodds(prob), false)) {
			//TODO exception
			

		}
                }
		current_origin = current_endpoint;
	}

	return true;
}

bool SonarOcTree::insertBeam(const base::samples::SonarBeam& beam,
		base::samples::RigidBodyState sonar_state) {

	octomath::Vector3 origin(0.0, 0.0, 0.0);
	//calculate the direction of the beam
	octomath::Vector3 direction(1.0, 0.0, 0.0);
	
	double length = beam.speed_of_sound * (beam.sampling_interval);

	double limit_horizontal = beam.beamwidth_horizontal * 180 / 3.14159265;
	double limit_vertical = beam.beamwidth_vertical * 180 / 3.14159265;

	octomath::Vector3 ray_direction = direction;

	octomath::Vector3 sonar_pos(sonar_state.position.x(),
			sonar_state.position.y(), sonar_state.position.z());

	octomath::Quaternion sonar_ori(sonar_state.orientation.w(),
			sonar_state.orientation.x(), sonar_state.orientation.y(),
			sonar_state.orientation.z());

	octomath::Pose6D sonar_tf(sonar_pos, sonar_ori);

	//TODO verificar se é necessario iterar em 0.5 de passo
	for (double i = -limit_horizontal / 2; i < limit_horizontal / 2; i++) {

		for (double j = -limit_vertical / 2; j < limit_vertical / 2; j++) {
			ray_direction = direction;
			ray_direction.rotate_IP(0, DEG2RAD(j), DEG2RAD(i)+beam.bearing.rad);
			insertBinsRay(beam.beam, origin, ray_direction, length);

		}

	}

	return true;
}

bool SonarOcTree::mergeTrees(octomap::SonarOcTree &tree2,octomap::point3d offset){
  tree2.expand();
  
  for (octomap::OcTree::leaf_iterator it = tree2.begin_leafs(), end = tree2.end_leafs(); it != end; ++it)
  {
	octomap::point3d tree2Point = it.getCoordinate();
	
	octomap::point3d rootPoint = offset + tree2Point; 

        // check occupancy prob:
        float tree2PointLogOdd = it->getLogOdds();
	this->updateNode(rootPoint,tree2PointLogOdd,false);
  }  
  
  return true;
  
}

double compareTrees(octomap::OcTree &rootTree,octomap::OcTree &localTree,octomap::point3d local2root)
{
    rootTree.expand();
    localTree.expand();

    double kld_sum = 0;
    for (octomap::OcTree::leaf_iterator it = localTree.begin_leafs(), end = localTree.end_leafs(); it != end; ++it)
    {
	octomap::point3d localPoint = it.getCoordinate();
	octomap::point3d rootPoint = local2root + localPoint; 
	octomap::OcTreeNode* rootNode = rootTree.search(rootPoint);
	if(!rootNode)
	    continue;
	
	double rootPointProb = rootPointProb = rootNode->getOccupancy();
	double localPointProb = it->getOccupancy();
	
	double kld = 0;
	if (localPointProb < 0.0001)
	    kld =log((1-rootPointProb)/(1-localPointProb))*(1-rootPointProb);
	else if (localPointProb > 0.9999)
	    kld =log(rootPointProb/localPointProb)*rootPointProb;
	else
	    kld +=log(rootPointProb/localPointProb)*rootPointProb + log((1-rootPointProb)/(1-localPointProb))*(1-rootPointProb);

	kld_sum+=kld;
    }
    
    if (isnan(kld_sum))
	throw std::logic_error("found NaN when comparing octomaps");
    return kld_sum;

}
