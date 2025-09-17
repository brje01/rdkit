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
        {1, 1.00},   // H
        {2, 1.40},   // He
        {3, 1.81},   // Li
        {4, 1.53},   // Be
        {5, 1.92},   // B
        {6, 1.70},   // C
        {7, 1.65},   // N
        {8, 1.60},   // O
        {9, 1.47},   // F
        {10, 1.54},  // Ne
        {11, 2.27},  // Na
        {12, 1.73},  // Mg
        {13, 1.84},  // Al
        {14, 2.10},  // Si
        {15, 1.80},  // P
        {16, 1.90},  // S
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
        std::vector<std::vector<int>> neighbour_list;
        neighbour_list.resize(mol.getNumAtoms());
        double v = 0.0;
        double o_pairwise = 0.0;
        double o_triple = 0.0;
        double o_quadruple = 0.0;
        double o_pentuple = 0.0;
        double o_sextuple = 0.0;
        double epsilon = 1.0; // Gaussian cutoff parameter
        const RDKit::Conformer &conf = mol.getConformer(confId);
        //get all atom positions
        std::vector<RDGeom::Point3D> all_positions;
        std::vector<double> alpha_vector;
        for (unsigned int i = 0; i < mol.getNumAtoms(); ++i) {
            RDGeom::Point3D pos = conf.getAtomPos(i);
            all_positions.push_back(pos);
            int atomic_num = mol.getAtomWithIdx(i)->getAtomicNum();
            alpha_vector.push_back(getAlpha(vdw_radii.at(atomic_num)));
            v += _p_ * std::pow(pi/alpha_vector[i],1.5);
            //identify all neighbours
            for (unsigned int j = 1; j < mol.getNumAtoms(); ++j){
              RDGeom::Point3D pos_j = conf.getAtomPos(j);
              double d = std::pow((pos.x - pos_j.x)*(pos.x - pos_j.x) + 
                         (pos.y - pos_j.y)*(pos.y - pos_j.y) + 
                         (pos.z - pos_j.z)*(pos.z - pos_j.z), 0.5);
              int atomic_num_j = mol.getAtomWithIdx(j)->getAtomicNum();
              if (d <= vdw_radii.at(atomic_num_j)+vdw_radii.at(atomic_num)+epsilon){
                neighbour_list[i].push_back(j);
              }
            }
        }

        for(unsigned int i=0;i<mol.getNumAtoms();++i){
            RDGeom::Point3D pos1 = all_positions[i];
            for(unsigned int j=1;j<i;++j){
                RDGeom::Point3D pos2 = all_positions[j];
                // check if j is a neighbour of i
                if (std::find(neighbour_list[i].begin(), neighbour_list[i].end(), j) != neighbour_list[i].end()){ 
                  // compute pairwise overlap
                  double tmp1 = 0.0;
                  tmp1 = pairwiseOverlap(pos1, alpha_vector[i],
                                      pos2, alpha_vector[j]);
                  o_pairwise += tmp1;
                  for(unsigned int k=2;k<j;++k){
                      RDGeom::Point3D pos3 = all_positions[k];
                      if (std::find(neighbour_list[j].begin(), neighbour_list[j].end(), k) != neighbour_list[j].end()){
                        o_triple += tripleOverlap(pos1, alpha_vector[i],
                                                  pos2, alpha_vector[j],
                                                  pos3, alpha_vector[k]);
                        for(unsigned int l=3;l<k;++l){
                            RDGeom::Point3D pos4 = all_positions[l];
                            if (std::find(neighbour_list[k].begin(), neighbour_list[k].end(), l) != neighbour_list[k].end()){
                              o_quadruple += quadrupleOverlap(pos1, alpha_vector[i],
                                                  pos2, alpha_vector[j],
                                                  pos3, alpha_vector[k],
                                                  pos4, alpha_vector[l]);
                              for(unsigned int m=4;m<l;++m){
                                  RDGeom::Point3D pos5 = all_positions[m];
                                  if (std::find(neighbour_list[l].begin(), neighbour_list[l].end(), m) != neighbour_list[l].end()){
                                    o_pentuple += pentupleOverlap(pos1, alpha_vector[i],
                                                        pos2, alpha_vector[j],
                                                        pos3, alpha_vector[k],
                                                        pos4, alpha_vector[l],
                                                        pos5, alpha_vector[m]);
                                    for(unsigned int n=5;n<m;++n){
                                        RDGeom::Point3D pos6 = all_positions[n];
                                        if (std::find(neighbour_list[m].begin(), neighbour_list[m].end(), n) != neighbour_list[m].end()){
                                          o_sextuple += sextupleOverlap(pos1, alpha_vector[i],
                                                              pos2, alpha_vector[j],
                                                              pos3, alpha_vector[k],
                                                              pos4, alpha_vector[l],
                                                              pos5, alpha_vector[m],
                                                              pos6, alpha_vector[n]);    
                                          }
                                        }
                                  }   
                              }
                            }
                          }
                        }
                      }
                    }
                }
        }
        std::cout << "v: " << v << std::endl;
        std::cout << "o_pairwise: " << o_pairwise << std::endl;
        std::cout << "o_triple: " << o_triple << std::endl;
        std::cout << "o_quadruple: " << o_quadruple << std::endl;
        std::cout << "o_pentuple: " << o_pentuple << std::endl;
        std::cout << "o_sextuple: " << o_sextuple << std::endl;
        return v - o_pairwise + o_triple - o_quadruple + o_pentuple - o_sextuple;
        //return o;
    }

} // end of Descriptors namespace
} // end of RDKit namespace

double getD2CutOff(){
  // return std::numeric_limits<double>::max();
  return 9.0;
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

double pairwiseOverlap(const RDGeom::Point3D& Ri, double alphai,
                       const RDGeom::Point3D& Rj, double alphaj){
  double alpha_sum = alphai + alphaj;
  double prefactor = std::pow(pi/alpha_sum, 1.5);
  double Rij_squared = (Ri.x - Rj.x)*(Ri.x - Rj.x) +
             (Ri.y - Rj.y)*(Ri.y - Rj.y) +
             (Ri.z - Rj.z)*(Ri.z - Rj.z);
  double exponent = - ((alphai * alphaj * Rij_squared)/alpha_sum);
  return _p2_ * prefactor * exp(exponent);
}

double tripleOverlap(const RDGeom::Point3D& Ri, double alphai,
                       const RDGeom::Point3D& Rj, double alphaj,
                       const RDGeom::Point3D& Rk, double alphak) {
  double alpha_ij = alphai + alphaj;
  double alpha_ijk = alpha_ij + alphak;

  // Weighted centers
  RDGeom::Point3D Rij;
  Rij.x = (alphai * Ri.x + alphaj * Rj.x) / alpha_ij;
  Rij.y = (alphai * Ri.y + alphaj * Rj.y) / alpha_ij;
  Rij.z = (alphai * Ri.z + alphaj * Rj.z) / alpha_ij;
  RDGeom::Point3D Rijk;
  Rijk.x = (alpha_ij * Rij.x + alphak * Rk.x) / alpha_ijk;
  Rijk.y = (alpha_ij * Rij.y + alphak * Rk.y) / alpha_ijk;
  Rijk.z = (alpha_ij * Rij.z + alphak * Rk.z) / alpha_ijk;

  // Q terms
  double Qij = (alphai * alphaj / alpha_ij) * (Ri - Rj).lengthSq();
  double Qijk = (alpha_ij * alphak / alpha_ijk) * (Rij - Rk).lengthSq();
  double Q = Qij + Qijk;

  // Prefactor
  double prefactor = std::pow(pi / alpha_ijk, 1.5);

  // Triple overlap
  return _p2_ * _p_ * prefactor * std::exp(-Q);
}

double quadrupleOverlap(const RDGeom::Point3D& Ri, double alphai,
                       const RDGeom::Point3D& Rj, double alphaj,
                       const RDGeom::Point3D& Rk, double alphak,
                       const RDGeom::Point3D& Rl, double alphal) {
  double alpha_ij = alphai + alphaj;
  double alpha_ijk = alpha_ij + alphak;
  double alpha_ijkl = alpha_ijk + alphal;

  // Weighted centers
  RDGeom::Point3D Rij;
  Rij.x = (alphai * Ri.x + alphaj * Rj.x) / alpha_ij;
  Rij.y = (alphai * Ri.y + alphaj * Rj.y) / alpha_ij;
  Rij.z = (alphai * Ri.z + alphaj * Rj.z) / alpha_ij;
  RDGeom::Point3D Rijk;
  Rijk.x = (alpha_ij * Rij.x + alphak * Rk.x) / alpha_ijk;
  Rijk.y = (alpha_ij * Rij.y + alphak * Rk.y) / alpha_ijk;
  Rijk.z = (alpha_ij * Rij.z + alphak * Rk.z) / alpha_ijk;
  RDGeom::Point3D Rijkl;
  Rijkl.x = (alpha_ijk * Rijk.x + alphal * Rl.x) / alpha_ijkl;
  Rijkl.y = (alpha_ijk * Rijk.y + alphal * Rl.y) / alpha_ijkl;
  Rijkl.z = (alpha_ijk * Rijk.z + alphal * Rl.z) / alpha_ijkl;

  // Q terms
  double Qij = (alphai * alphaj / alpha_ij) * (Ri - Rj).lengthSq();
  double Qijk = (alpha_ij * alphak / alpha_ijk) * (Rij - Rk).lengthSq();
  double Qijkl = (alpha_ijk * alphal / alpha_ijkl) * (Rijk - Rl).lengthSq();
  double Q = Qij + Qijk + Qijkl;

  // Prefactor
  double prefactor = std::pow(pi / alpha_ijkl, 1.5);

  // Quadruple overlap
  return _p2_ * _p2_ * prefactor * std::exp(-Q);
}

double pentupleOverlap(const RDGeom::Point3D& Ri, double alphai,
                       const RDGeom::Point3D& Rj, double alphaj,
                       const RDGeom::Point3D& Rk, double alphak,
                       const RDGeom::Point3D& Rl, double alphal,
                       const RDGeom::Point3D& Rm, double alpham) {
  double alpha_ij = alphai + alphaj;
  double alpha_ijk = alpha_ij + alphak;
  double alpha_ijkl = alpha_ijk + alphal;
  double alpha_ijklm = alpha_ijkl + alpham;

  // Weighted centers
  RDGeom::Point3D Rij;
  Rij.x = (alphai * Ri.x + alphaj * Rj.x) / alpha_ij;
  Rij.y = (alphai * Ri.y + alphaj * Rj.y) / alpha_ij;
  Rij.z = (alphai * Ri.z + alphaj * Rj.z) / alpha_ij;
  RDGeom::Point3D Rijk;
  Rijk.x = (alpha_ij * Rij.x + alphak * Rk.x) / alpha_ijk;
  Rijk.y = (alpha_ij * Rij.y + alphak * Rk.y) / alpha_ijk;
  Rijk.z = (alpha_ij * Rij.z + alphak * Rk.z) / alpha_ijk;
  RDGeom::Point3D Rijkl;
  Rijkl.x = (alpha_ijk * Rijk.x + alphal * Rl.x) / alpha_ijkl;
  Rijkl.y = (alpha_ijk * Rijk.y + alphal * Rl.y) / alpha_ijkl;
  Rijkl.z = (alpha_ijk * Rijk.z + alphal * Rl.z) / alpha_ijkl;
  RDGeom::Point3D Rijklm;
  Rijklm.x = (alpha_ijkl * Rijkl.x + alpham * Rm.x) / alpha_ijklm;
  Rijklm.y = (alpha_ijkl * Rijkl.y + alpham * Rm.y) / alpha_ijklm;
  Rijklm.z = (alpha_ijkl * Rijkl.z + alpham * Rm.z) / alpha_ijklm;

  // Q terms
  double Qij = (alphai * alphaj / alpha_ij) * (Ri - Rj).lengthSq();
  double Qijk = (alpha_ij * alphak / alpha_ijk) * (Rij - Rk).lengthSq();
  double Qijkl = (alpha_ijk * alphal / alpha_ijkl) * (Rijk - Rl).lengthSq();
  double Qijklm = (alpha_ijkl * alpham / alpha_ijklm) * (Rijkl - Rm).lengthSq();
  double Q = Qij + Qijk + Qijkl + Qijklm; 

  // Prefactor
  double prefactor = std::pow(pi / alpha_ijklm, 1.5);

  // Pentuple overlap
  return _p2_ * _p2_ * _p_ * prefactor * std::exp(-Q);
}

double sextupleOverlap(const RDGeom::Point3D& Ri, double alphai,
                       const RDGeom::Point3D& Rj, double alphaj,
                       const RDGeom::Point3D& Rk, double alphak,
                       const RDGeom::Point3D& Rl, double alphal,
                       const RDGeom::Point3D& Rm, double alpham,
                       const RDGeom::Point3D& Rn, double alphan) {
  double alpha_ij = alphai + alphaj;
  double alpha_ijk = alpha_ij + alphak;
  double alpha_ijkl = alpha_ijk + alphal;
  double alpha_ijklm = alpha_ijkl + alpham;
  double alpha_ijklmn = alpha_ijklm + alphan;

  // Weighted centers
  RDGeom::Point3D Rij;
  Rij.x = (alphai * Ri.x + alphaj * Rj.x) / alpha_ij;
  Rij.y = (alphai * Ri.y + alphaj * Rj.y) / alpha_ij;
  Rij.z = (alphai * Ri.z + alphaj * Rj.z) / alpha_ij;
  RDGeom::Point3D Rijk;
  Rijk.x = (alpha_ij * Rij.x + alphak * Rk.x) / alpha_ijk;
  Rijk.y = (alpha_ij * Rij.y + alphak * Rk.y) / alpha_ijk;
  Rijk.z = (alpha_ij * Rij.z + alphak * Rk.z) / alpha_ijk;
  RDGeom::Point3D Rijkl;
  Rijkl.x = (alpha_ijk * Rijk.x + alphal * Rl.x) / alpha_ijkl;
  Rijkl.y = (alpha_ijk * Rijk.y + alphal * Rl.y) / alpha_ijkl;
  Rijkl.z = (alpha_ijk * Rijk.z + alphal * Rl.z) / alpha_ijkl;
  RDGeom::Point3D Rijklm;
  Rijklm.x = (alpha_ijkl * Rijkl.x + alpham * Rm.x) / alpha_ijklm;
  Rijklm.y = (alpha_ijkl * Rijkl.y + alpham * Rm.y) / alpha_ijklm;
  Rijklm.z = (alpha_ijkl * Rijkl.z + alpham * Rm.z) / alpha_ijklm;
  RDGeom::Point3D Rijklmn;
  Rijklmn.x = (alpha_ijklm * Rijklm.x + alphan * Rn.x) / alpha_ijklmn;
  Rijklmn.y = (alpha_ijklm * Rijklm.y + alphan * Rn.y) / alpha_ijklmn;
  Rijklmn.z = (alpha_ijklm * Rijklm.z + alphan * Rn.z) / alpha_ijklmn;
  // Q terms
  double Qij = (alphai * alphaj / alpha_ij) * (Ri - Rj).lengthSq();
  double Qijk = (alpha_ij * alphak / alpha_ijk) * (Rij - Rk).lengthSq();
  double Qijkl = (alpha_ijk * alphal / alpha_ijkl) * (Rijk - Rl).lengthSq();
  double Qijklm = (alpha_ijkl * alpham / alpha_ijklm) * (Rijkl - Rm).lengthSq();
  double Qijklmn = (alpha_ijklm * alphan / alpha_ijklmn) * (Rijklm - Rn).lengthSq();
  double Q = Qij + Qijk + Qijkl + Qijklm + Qijklmn; 
  // Prefactor
  double prefactor = std::pow(pi / alpha_ijklmn, 1.5);
  // Sextuple overlap
  return _p2_ * _p2_ * _p2_ * prefactor * std::exp(-Q);
}