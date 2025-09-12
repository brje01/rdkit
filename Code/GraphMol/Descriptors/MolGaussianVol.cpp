#include <GraphMol/RDKitBase.h>
#include <GraphMol/MolTransforms/MolTransforms.h>
#include <Geometry/point.h>

#include "MolGaussianVol.h"

#include <cmath>
#include <vector>
#include <map>
#include <set>
#include <cstring>
#include <iostream>
#include <algorithm>

namespace RDKit {
namespace Descriptors {
namespace {

double getD2CutOff(){
  return std::numeric_limits<double>::max();
}

double getA_ak( double alpha1, double alpha2 ){
  double a_ak = static_cast< double >( 0 );
  double sum = alpha1 + alpha2;
  if( sum > static_cast<double>( 0 ) ){
    a_ak = alpha1*alpha2 / sum;
  }
  return a_ak;
}

double getOverlap0( double alpha1, double alpha2 ){
  double a_pi = static_cast< double>( 0 );
  double sum = alpha1 + alpha2;
  if( sum > static_cast<double>( 0 ) ){
    a_pi = pi / sum;
  }
  double a_oe = _p2_ * a_pi * std::sqrt( a_pi );
  return a_oe;
}

void ComputeGaussianVolume(
                            const ROMol &mol, 
                            int confId,
                            double &o
						 )
{
    o = 0.0;
    const RDKit::Conformer &conf = mol.getConformer(confId);
    //get all atom positions
    std::vector<RDGeom::Point3D> all_positions;
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        RDGeom::Point3D pos = conf.getAtomPos(i);
        all_positions.push_back(pos);
    }
    // loop over all atoms
    for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
        RDGeom::Point3D pos1 = all_positions[i];
        for (unsigned int j = 0; j < i; ++j) {
            RDGeom::Point3D pos2 = all_positions[j];
            double d2 = (pos1.x - pos2.x)*(pos1.x - pos2.x) + 
                        (pos1.y - pos2.y)*(pos1.y - pos2.y) + 
                        (pos1.z - pos2.z)*(pos1.z - pos2.z);
            if ( d2 < getD2CutOff() ) {
                double alpha1 = alpha_fit_vector[i];
                double alpha2 = alpha_ref_vector[j];
                double a_ak = getA_ak( alpha1, alpha2 );
                double a_oe = getOverlap0( alpha1, alpha2 );
                double expMinus_ = exp( - d2 * a_ak );
                o += a_oe * expMinus_;
            }
    }
} // end of anonymous namespace
} // end of Descriptors namespace
} // end of RDKit namespace