#include "SonarOcTree.hpp"
#include <Eigen/Dense>
#include <sys/time.h>
#include <octomap/math/Vector3.h>
#include <octomap/math/Quaternion.h>
#include <octomap/math/Pose6D.h>

using namespace octomap;


const double value_trh = 30;
const double empty_prob = 0.1;
const double full_prob = 0.8;

const float logit_full_prob = logodds(full_prob);
const float logit_empty_prob = logodds(empty_prob);

inline void sph2cart(base::Vector3d& cartvec, double radius, double theta, double phi){
	cartvec[0]=radius*cos(phi)*cos(theta);
	cartvec[1]=radius*cos(phi)*sin(theta);
	cartvec[2]=radius*sin(phi);
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
  
  
octomath::Vector3 eigen2octoVector(base::Vector3d vec)
{
  return octomath::Vector3(vec[0],vec[1],vec[2]);
}

octomath::Quaternion eigen3octoQuaternion(base::Quaterniond q)
{
  return octomath::Quaternion(q.w(),q.x(),q.y(),q.z());
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


bool SonarOcTree::CreateBin(int bin, double bearing, float offset, double alpha, double beta  ){
  base::samples::RigidBodyState sonar_state;
  base::Pose sonar_pose;
  //base::Matrix3d Pantilt;
  
  sonar_pose.orientation = Eigen::AngleAxisd(alpha, Eigen::MatrixBase<base::Vector3d>::UnitZ())
	    * Eigen::AngleAxisd(beta, Eigen::MatrixBase<base::Vector3d>::UnitY())
	    * Eigen::AngleAxisd(bearing, Eigen::MatrixBase<base::Vector3d>::UnitZ());
	    
  
  
  sonar_state.setPose(sonar_pose);
	    
  return this->CreateBin( bin, 0, &octomap::SonarOcTree::updater, sonar_state );
}

bool SonarOcTree::CreateBin( int bin, double bearing,  pUpdateMethod fnode, base::samples::RigidBodyState sonar_state ){
    
  sonar_state.orientation = sonar_state.orientation * Eigen::AngleAxisd(bearing, Eigen::MatrixBase<base::Vector3d>::UnitZ());
  
  base::Vector3d Backcorners[4];
  base::Vector3d Frontcorners[4];
  
  int Nrow = matvar->dims[0];
  int Ncol = matvar->dims[1];
  
  //TODO conferir essas constantes
  double thetalimits[]={ThetaMin,ThetaMax};
  double philimits[]={PhiMin,PhiMax};
  
  double thetaM = (thetalimits[1]+thetalimits[0])/2;
  double phiM = (philimits[1]+philimits[0])/2;
  


			  
  for(int theta_index=0; theta_index<2; theta_index++)
  for(int phi_index=0; phi_index<2; phi_index++)
  {
  
    //From sonar perpective
    sph2cart(Backcorners[(theta_index << 1) | phi_index],
    rbinlimits[bin],
    thetalimits[theta_index],
    philimits[phi_index]);	
    
    sph2cart(Frontcorners[(theta_index << 1) | phi_index],
    (rbinlimits[bin+1]/cos(thetaM-thetalimits[theta_index]))/cos(phiM-philimits[phi_index]),
    thetalimits[theta_index],
    philimits[phi_index]);

    //from inertial perspective
    Backcorners[(theta_index << 1) | phi_index] = sonar_state.getTransform() * Backcorners[(theta_index << 1) | phi_index];
    
    
    Frontcorners[(theta_index << 1) | phi_index] = sonar_state.getTransform() * Frontcorners[(theta_index << 1) | phi_index];
    
  }
  
  double rangexyz[3][2];
  //Initial bounding box guess
  for(int axis=0; axis<3; axis++){
    rangexyz[axis][0]=std::min(Frontcorners[0](axis),Backcorners[0](axis));
    rangexyz[axis][1]=std::max(Frontcorners[0](axis),Backcorners[0](axis));
  }
  
  


  //Bounding box computation
  for(int corner = 1; corner < 4; corner++){
    updaterange(rangexyz,Frontcorners[corner]);
    updaterange(rangexyz,Backcorners[corner]);
  }

  
  
  //Origin (minor corner of the bounding box) and side discretization (number of box between corners)

  octomap::OcTreeKey startkey = this->coordToKey(rangexyz[0][0],rangexyz[1][0],rangexyz[2][0]);
  octomap::OcTreeKey endkey = this->coordToKey(rangexyz[0][1],rangexyz[1][1],rangexyz[2][1]);
  
  /* filling octomap & angle to matrix index convertions */

  const double kstep = 1800/M_PI;
  Eigen::Affine3d PantiltInverse = sonar_state.getTransform().inverse();  
  
  
/*  
 * double normalizer=0;
  struct timeval start, end;
  gettimeofday(&start, NULL);*/

  OcTreeNode* actualnode;
  
  for(int stepX=startkey[0]; stepX<=endkey[0]; stepX++)
  for(int stepY=startkey[1]; stepY<=endkey[1]; stepY++)
  for(int stepZ=startkey[2]; stepZ<=endkey[2]; stepZ++){
    octomap::OcTreeKey stepkey = octomap::OcTreeKey(stepX,stepY,stepZ);
    octomap::point3d coord = this->keyToCoord(stepkey);
    base::Vector3d coordfromsonar = PantiltInverse * base::Vector3d(coord.x(),coord.y(),coord.z());
    
    double radius = coordfromsonar.norm();

    
    if(radius < rbinlimits[bin] || radius > rbinlimits[bin+1])
      continue;

    double theta = atan2( coordfromsonar[1], coordfromsonar[0]);
    double phi = asin(coordfromsonar[2]/radius);

    int col = trunc((theta-ThetaMin)*kstep);
    if(col < 0 || col >= Ncol)	
      continue;

    int row = trunc((phi-PhiMin)*kstep);
    if(row < 0 || row >= Nrow)
      continue; 
    
    
    if(((double*)matvar->data)[row + Nrow*col]==0.0)
      continue;
    
//     normalizer+=((double*)matvar->data)[row + Nrow*col];
    
    if((actualnode = this->search(stepkey)) == NULL){
      this->updateNode(stepkey, (float) 0 );
      actualnode = this->search(stepkey);
      
    }
  //  else
  
    (this->*fnode)(actualnode,&(((double*)matvar->data)[row + Nrow*col]));
    
     // actualnode->addValue((float) (*logitprob * ((double*)matvar->data)[row + Nrow*col]/30.0) );      //logodds(((double*)matvar->data)[row + Nrow*col]));
    //this->updateNode(stepkey, (float) (*logitprob + (((double*)matvar->data)[row + Nrow*col]/30.0 - 1)*2.2));
  }
    
  /*  
  std::cout<< "Normalizing" << std::endl;
  
  float updater = logodds(1.0/(1.0+normalizer));
  
   for (OcTree::leaf_iterator it = this->begin_leafs(), end = this->end_leafs(); it != end; ++it)
     this->updateNode(it.getKey(),updater);
  
  */
  
/*  
  gettimeofday(&end, NULL);
  double voxels = (endkey[0] - startkey[0])*(endkey[1] - startkey[1])*(endkey[2] - startkey[2]);
  unsigned long long useconds = end.tv_sec - start.tv_sec;
  useconds*=1000000;
  useconds+=end.tv_usec - start.tv_usec;

  std::cout<< voxels/useconds << " Mega voxels per second in " << end.tv_sec - start.tv_sec << "s" << std::endl;
  std::cout<< "DONE" << std::endl << std::endl;
  */
   
  return true;
}

inline void SonarOcTree::updater(OcTreeNode* actualnode, double* gain){
  std::cout << "updater" << std::endl;
  float k = (float) ((*logitprob) * (*gain) /30.0);
  std::cout << "floated" << k << std::endl;
  if(actualnode==NULL)
    std::cout << "NUUUUUL" << std::endl;
  else
    std::cout << "FUUULLL" << actualnode << std::endl;
  actualnode->addValue((float) ((*logitprob) * (*gain) /30.0) );
  std::cout << "updated" << std::endl;
  
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
				 logodds(prob), false)) {
			//TODO exception
			

		}
                }
		current_origin = current_endpoint;
	}

	return true;
}

bool SonarOcTree::insertRealBeam(const base::samples::SonarBeam& beam,
		base::samples::RigidBodyState& sonar_state){
  
  double length = beam.speed_of_sound * beam.sampling_interval;
  
  for(int i = 0; i <= beam.beam.size(); i++)
    rbinlimits[i] =  i * length;
  
  
  #pragma omp parallel num_threads(6)
  {
    
  #pragma omp for schedule(dynamic)
  for(int i = 0; i < beam.beam.size(); i++){
//     if (beam.beam[i] == 0)
//       continue;
    
    float poweroffset = beam.beam[i]*80.0/255.0;
  
    if(poweroffset<value_trh)
      logitprob = &logit_empty_prob;
    else
      logitprob = &logit_full_prob;
    
    this->CreateBin( i, beam.bearing.rad, &octomap::SonarOcTree::updater, sonar_state );
    
  }
  
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


bool SonarOcTree::mergeTrees(SonarOcTree &tree2,point3d offset){
  tree2.expand();
  
  for (OcTree::leaf_iterator it = tree2.begin_leafs(), end = tree2.end_leafs(); it != end; ++it)
  {
	point3d tree2Point = it.getCoordinate();
	
	point3d root_point = offset + tree2Point; 

        // check occupancy prob:
        float tree2PointLogOdd = it->getLogOdds();
	this->updateNode(root_point,tree2PointLogOdd,false);
  }  
  
  return true;
  
}

double SonarOcTree::compareTrees(const SonarOcTree& tree, base::Vector3d& tree_position, base::Quaterniond& tree_orientatation)
{
    double kld_sum = 0;
    for (OcTree::leaf_iterator it = tree.begin_leafs(), end = tree.end_leafs(); it != end; ++it)
    {
	point3d tree_point = it.getCoordinate();
	
	pose6d tree_transfom = octomath::Pose6D(eigen2octoVector(tree_position),
							     eigen3octoQuaternion(tree_orientatation)); 
	
	
	point3d root2Tree = tree_transfom.transform(tree_point);

	OcTreeNode* root_node = this->search(root2Tree);
	if(!root_node)
	    continue;
	
	double root_point_prob = root_node->getOccupancy();
	double tree_point_prob = it->getOccupancy();
	
	double kld = 0;
	if (tree_point_prob < 0.0001)
	    kld =log((1-root_point_prob)/(1-tree_point_prob))*(1-root_point_prob);
	else if (tree_point_prob > 0.9999)
	    kld =log(root_point_prob/tree_point_prob)*root_point_prob;
	else
	    kld +=log(root_point_prob/tree_point_prob)*root_point_prob + log((1-root_point_prob)/(1-tree_point_prob))*(1-root_point_prob);

	kld_sum+=kld;
    }
    
    if (isnan(kld_sum))
	throw std::logic_error("found NaN when comparing octomaps");
    return kld_sum;
}

double octomap::SonarOcTree::evaluateSonarBeam(const base::samples::RigidBodyState& particle_pose, const base::samples::SonarBeam& sonar_beam)
{
  double weight = 0;
  sum_decay_occ = 0;
  sum_decay2 = 0;
  
           
  std::vector<uint8_t> projected_beam;
  
  
  for(int current_bin = 0; current_bin<sonar_beam.beam.size(); current_bin++)
  {
      //funcao do elael recebendo current_bin, sonar_beam.bearing, rbs, tlvz o valor do bin e o ponteiro para minha funcao 
  }    
  
  double projected_intensity = sum_decay_occ/sum_decay2;
  
  projected_intensity = floor(projected_intensity);
  
  projected_beam.push_back(projected_intensity);
  
  
 
  return weight;
}

void SonarOcTree::calculateIntensity(OcTreeNode* node, double* decay)
{
  double occupancy = node->getOccupancy();
  sum_decay_occ += *decay*occupancy;
  sum_decay2 += (*decay)*(*decay);

}


