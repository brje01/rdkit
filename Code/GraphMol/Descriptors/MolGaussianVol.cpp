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


const std::map<unsigned int, double> vdw_radii = {
        {1, 1.10},   // H
        {2, 1.40},   // He
        {3, 1.81},   // Li
        {4, 1.53},   // Be
        {5, 1.92},   // B
        {6, 1.70},   // C
        {7, 1.55},   // N
        {8, 1.52},   // O
        {9, 1.47},   // F
        {10, 1.54},  // Ne
        {11, 2.27},  // Na
        {12, 1.73},  // Mg
        {13, 1.84},  // Al
        {14, 2.10},  // Si
        {15, 1.80},  // P
        {16, 1.80},  // S
        {17, 1.75},  // Cl
        {18, 1.88},  // Ar
        {19, 2.75},  // K
        {20, 2.31},  // Ca
        {31, 1.87},  // Ga
        {32, 2.11},  // Ge
        {33, 1.85},  // As
        {34, 1.90},  // Se
        {35, 1.83},  // Br
        {36, 2.02},  // Kr
        {37, 3.03},  // Rb
        {38, 2.49},  // Sr
        {49, 1.93},  // In
        {50, 2.17},  // Sn
        {51, 2.06},  // Sb
        {52, 2.06},  // Te
        {53, 1.98},  // I
        {54, 2.16},  // Xe
        {55, 3.43},  // Cs
        {56, 2.68},  // Ba
        {81, 1.96},  // Tl
        {82, 2.02},  // Pb
        {83, 2.07},  // Bi
        {84, 1.97},  // Po
        {85, 2.02},  // At
        {86, 2.20},  // Rn
        {87, 3.48},  // Fr
        {88, 2.83},  // Ra
    };

namespace RDKit {
namespace Descriptors {
    double ComputeGaussianVolume(
                            const ROMol &mol, 
                            int confId
						 )
    {
        double o = 0.0;
        const RDKit::Conformer &conf = mol.getConformer(confId);
        //get all atom positions
        std::vector<RDGeom::Point3D> all_positions;
        std::vector<double> alpha_vector;
        for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
            RDGeom::Point3D pos = conf.getAtomPos(i);
            all_positions.push_back(pos);
            int atomic_num = mol.getAtomWithIdx(i)->getAtomicNum();
            alpha_vector.push_back(getAlpha(vdw_radii.at(atomic_num)));
        }
        // loop over all atoms
        for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
            RDGeom::Point3D pos1 = all_positions[i];
            for (unsigned int j = 0; j < mol.getNumAtoms(); ++j) {
                RDGeom::Point3D pos2 = all_positions[j];
                double d2 = (pos1.x - pos2.x)*(pos1.x - pos2.x) + 
                            (pos1.y - pos2.y)*(pos1.y - pos2.y) + 
                            (pos1.z - pos2.z)*(pos1.z - pos2.z);
                if ( d2 < getD2CutOff() ) {
                    double alpha1 = alpha_vector[i];
                    double alpha2 = alpha_vector[j];
                    double a_ak = getA_ak( alpha1, alpha2 );
                    double a_oe = getOverlap0( alpha1, alpha2 );
                    double expMinus_ = exp( - d2 * a_ak );
                    o += a_oe * expMinus_;
                }
            }
        }
        return o;
    }

} // end of Descriptors namespace
} // end of RDKit namespace

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

double getAlpha( double r ){
  double alpha( static_cast<double>( 0 ) );
  if( r > static_cast<double>( 0 ) ){
    alpha = kappa/( r*r );
  }
  return alpha;
}