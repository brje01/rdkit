#ifndef RDKit_Descriptors_MolGaussianVol_H
#define RDKit_Descriptors_MolGaussianVol_H

#include <RDGeneral/export.h>
#include "RegisterDescriptor.h"
#include <cmath>

// from the original shape paper and the pubchem implementation of it
inline constexpr double pi = 3.14159265358979323846;
inline constexpr double _p_ = 2.7;
inline constexpr double _p2_ = _p_ * _p_;
inline constexpr double lambda = (4.0 * pi) / (3.0 * _p_);
inline constexpr double kappa = pi / std::pow(lambda, 2.0 / 3.0);

namespace RDKit {
class ROMol;
namespace Descriptors {  
  RDKIT_DESCRIPTORS_EXPORT double ComputeGaussianVolume(const ROMol &mol, 
                                                     int confId);

} // namespace Descriptors
} // namespace RDKit

double getD2CutOff();
double getAlpha(double radius);
double getOverlap0(double alpha1, double alpha2);
double getA_ak(double alpha1, double alpha2);
double pairwiseOverlap(const RDGeom::Point3D& Ri, double alphai,
                       const RDGeom::Point3D& Rj, double alphaj);
double tripleOverlap(const RDGeom::Point3D& Ri, double alphai,
                     const RDGeom::Point3D& Rj, double alphaj,
                     const RDGeom::Point3D& Rk, double alphak);

#endif