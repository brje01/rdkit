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

void ComputeGaussianVolume(
                            const ROMol &mol, 
                            int confId,
                            double &o
						 )
{
    o = 0.0;
    double alpha = alpha_fit_vector[i];
	double a_ak = getA_ak( alpha, alpha );
	double a_oe = getOverlap0( alpha, alpha );
	double expMinus_ = exp( - d2 * a_ak );
	o += a_oe * expMinus_;
}

} // end of anonymous namespace
} // end of Descriptors namespace
} // end of RDKit namespace