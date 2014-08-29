#include <sonaroctomap/SonarOcTree.hpp>
#include "matio.h"
#include <Eigen/Dense>
#include <octomap/math/Utils.h>

using namespace octomap;



inline void sph2cart(base::Vector3d& cartvec, double radius, double theta, double phi){
	cartvec[0]=radius*cos(phi)*cos(theta);
	cartvec[1]=radius*cos(phi)*sin(theta);
	cartvec[2]=radius*sin(phi);
// 	std::cout << "carttovec: " << std::endl << cartvec << "..." <<std::endl;
}

template<typename T>	
void coutarray(T s, int i, int j)
{
    for (int ii = 0; ii < i; ii++){
    for (int jj = 0; jj < j; jj++) 
    std::cout << s[ii][jj] << " "; 
    std::cout << std::endl;
    }
}


inline void updaterange(double range[][2],const base::Vector3d& source){
for(int axis = 0; axis < 3; axis++){
  range[axis][0]=std::min(range[axis][0],source(axis));
  range[axis][1]=std::max(range[axis][1],source(axis));
}
}
  
  
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

const double PhiMin= -M_PI*686/1800;
const double PhiMax= M_PI*874/1800;
const double ThetaMin = -M_PI*315/1800;
const double ThetaMax = M_PI*290/1800;

bool SonarOcTree::CreateBinPointCloud(double octo_resolution, std::string filename, std::string varname, int bin){
  mat_t *openmatfp;
  matvar_t *matvar;
  openmatfp = Mat_Open(filename.c_str(),MAT_ACC_RDONLY);
  matvar = Mat_VarRead(openmatfp,varname.c_str());
  std::cout << ((double*)matvar->data)[0] << std::endl;
  
//   Mat_VarPrint(matvar,0);
  
  
  
  
  base::Matrix3d Pantilt;
  base::Vector3d Backcorners[4];
  base::Vector3d Frontcorners[4];
  double tangentplanegain;
  double alpha = 0;
  double beta = 0;
  int Nrow = matvar->dims[0];
  int Ncol = matvar->dims[1];
  
  Pantilt = Eigen::AngleAxisd(alpha, Eigen::MatrixBase<base::Vector3d>::UnitZ())
	    * Eigen::AngleAxisd(beta, Eigen::MatrixBase<base::Vector3d>::UnitY());

  std::cout << Pantilt << std::endl;
  //TODO conferir essas constantes
  double thetalimits[]={ThetaMin,ThetaMax};
  double philimits[]={PhiMin,PhiMax};
  double rbinlimits[]={1,4,6,8,10,12,14,16};
  
  double thetaM = (thetalimits[1]+thetalimits[0])/2;
  double phiM = (philimits[1]+philimits[0])/2;

  std::cout<< "first loop"<<std::endl 
			  <<"rbinlimits[bin:bin+1] "<<  rbinlimits[bin] << ":" << rbinlimits[bin+1] << std::endl
			  <<"thetalimits[0:1] "<<  thetalimits[0] << ":" << thetalimits[1] << std::endl
			  <<"philimits[0:1] "<<  philimits[0] << ":" << philimits[1] << std::endl
			  <<"thetaM "<<  thetaM << std::endl
			  <<"phiM "<<  phiM << std::endl << std::endl;

  for(int theta_index=0; theta_index<2; theta_index++)
  for(int phi_index=0; phi_index<2; phi_index++)
  {
  
    //From sonar perpective
    sph2cart(Backcorners[(theta_index << 1) | phi_index],
    rbinlimits[bin],
    thetalimits[theta_index],
    philimits[phi_index]);	
/*
    std::cout<< "dT=" << thetaM-thetalimits[theta_index] << ", dP=" << phiM-philimits[phi_index] << std::endl;
    std::cout<< "abs(cos(dT))=" << abs(cos(thetaM-thetalimits[theta_index])) << ", cos(dP)=" << cos(phiM-philimits[phi_index]) << std::endl;*/
    sph2cart(Frontcorners[(theta_index << 1) | phi_index],
    (rbinlimits[bin+1]/cos(thetaM-thetalimits[theta_index]))/cos(phiM-philimits[phi_index]),
    thetalimits[theta_index],
    philimits[phi_index]);

    //from inertial perspective
    Backcorners[(theta_index << 1) | phi_index] = Pantilt * Backcorners[(theta_index << 1) | phi_index];
    //std::cout << theta_index << "," <<phi_index << " Backcorners["<< ((theta_index << 1) + phi_index) <<"] " << std::endl << Backcorners[(theta_index << 1) + phi_index] << std::endl << std::endl;
    
    
    Frontcorners[(theta_index << 1) | phi_index] = Pantilt * Frontcorners[(theta_index << 1) | phi_index];
    
  }
    std::cout << "Backcorners: [" << std::endl;
    for (int k = 0; k < 4; k++) 
    std::cout << Backcorners[k] << std::endl << std::endl;
    std::cout << "]" << std::endl;
    
        std::cout << "Frontcorners: [" << std::endl;
    for (int k = 0; k < 4; k++) 
    std::cout << Frontcorners[k] << std::endl <<  std::endl;
    std::cout << "]" << std::endl;
  

  
  
  double rangexyz[3][2];


  std::cout<< "first box guess" << std::endl;
  std::cout<< "min(Frontcorners[0](1),Backcorners[0](1))" << std::endl;
  std::cout<< "min("<<Frontcorners[0](1)<<","<<Backcorners[0](1)<<") = "<< std::min(Frontcorners[0](1),Backcorners[0](1)) << std::endl;
  
  //Initial bounding box guess
  for(int axis=0; axis<3; axis++){
    rangexyz[axis][0]=std::min(Frontcorners[0](axis),Backcorners[0](axis));
    rangexyz[axis][1]=std::max(Frontcorners[0](axis),Backcorners[0](axis));
  }
  std::cout << "rangexyz = " << std::endl;
  coutarray(rangexyz,3,2);
  std::cout << std::endl;

  std::cout<< "computing box" << std::endl;
  //Bounding box computation
  for(int corner = 1; corner < 4; corner++){
    updaterange(rangexyz,Frontcorners[corner]);
    updaterange(rangexyz,Backcorners[corner]);
  }
  
  std::cout << "rangexyz = " << std::endl;
  coutarray(rangexyz,3,2);
  std::cout << std::endl;



  std::cout<< "side discretization" << std::endl;
  //Origin (minor corner of the bounding box) and side discretization (number of box between corners)

  octomap::OcTreeKey startkey = this->coordToKey(rangexyz[0][0],rangexyz[1][0],rangexyz[2][0]);
  octomap::OcTreeKey endkey = this->coordToKey(rangexyz[0][1],rangexyz[1][1],rangexyz[2][1]);
  
  /* filling octomap & angle to matrix index convertions */

  const double kstep = 1800/M_PI;
  bool dummy=false;
  double old=0;
  
  std::cout<< "(0,0) - " << ((double*)matvar->data)[((int)trunc((0.001-PhiMin)*kstep)) +  Nrow*(int)trunc((0.001-ThetaMin)*kstep)] << std::endl;
  std::cout<< "(0,0)[686,315] - " << ((double*)matvar->data)[686 + Nrow*315] << std::endl;
  
  std::cout<< "Octofill (" << Nrow << "x" << Ncol << ") from "<< this->keyToCoord(startkey) << " and to "<< this->keyToCoord(endkey) << std::endl;
  
  for(int stepX=startkey[0]; stepX<=endkey[0]; stepX++)
  for(int stepY=startkey[1]; stepY<=endkey[1]; stepY++)
  for(int stepZ=startkey[2]; stepZ<=endkey[2]; stepZ++){
    octomap::OcTreeKey stepkey = octomap::OcTreeKey(stepX,stepY,stepZ);
    octomap::point3d coord = this->keyToCoord(stepkey);
    
    
    double radius = coord.norm();

    
    if(radius < rbinlimits[bin] || radius > rbinlimits[bin+1])
      continue;

    double phi = atan2( coord.y(), coord.x());

    double theta = asin(coord.z()/radius);

    int col = trunc((theta-ThetaMin)*kstep);
    if(col < 0 || col >= Ncol)	
      continue;

    int row = trunc((phi-PhiMin)*kstep);
    if(row < 0 || row >= Nrow)
      continue; 

    
    if(((double*)matvar->data)[row + Nrow*col]==0.0)
      continue;
    
    this->updateNode(stepkey, (float) 	((double*)matvar->data)[row + Nrow*col]);
  }
    
  //std::cout<< "MEUS TESTE -> " << ((double*)matvar->data)[matvar->dims[0]*2+1] << std::endl;
  //Mat_VarPrint(matvar,1);
  std::cout<< "DONE" << std::endl;
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
