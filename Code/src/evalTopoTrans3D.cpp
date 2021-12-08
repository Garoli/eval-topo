/**
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Lesser General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 **/

/**
 * @file images/exampleRigidtransformation2d.cpp
 * @ingroup Examples
 * @author Kacper Pluta (\c kacper.pluta@esiee.fr )
 * Laboratoire d'Informatique Gaspard-Monge - LIGM, France
 *
 * @date 2014/06/28
 *
 * An example file named rigidtransformation2d.
 *
 * This file is part of the DGtal library.
 */

/**
*  Example of 2D rigid transformation using forward and backward model.
   @see @ref moduleGeometricTransform
   \image html church_backward.jpg "Result for backward model" 
*  \example images/exampleRigidtransformation2d.cpp
**/


///////////////////////////////////////////////////////////////////////////////
#include <cmath>
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"
//! [include]
#include "Interpolation.h"
#include "AffineTransformation.h"
#include "ImageTransformation.h"
///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace functors;
using namespace Z3i;

///////////////////////////////////////////////////////////////////////////////

int main(int, char**){
    //image
    typedef unsigned char Value;
    typedef ImageContainerBySTLVector<Domain, Value> Image;

    //transformaton
    typedef AffineTransformation3D<Space> AffineTransformation;
    typedef DomainTransformer3D<Image::Domain, AffineTransformation> DomainTransformer;
    typedef NearestNeighbour<Space, Image> Interpolation;
    typedef ImageTransformation<Space, Image, AffineTransformation,
                                  DomainTransformer, Interpolation> ImageTransformation;

    //open image
    Image image = GenericReader<Image>::import("samples/Al100.vol");

    //instanciate
    Interpolation interpolation(&image);
    AffineTransformation transformation = AffineTransformation::rotation(RealPoint(50, 50, 50), RealPoint(1, 0, 0), M_PI_4);
    DomainTransformer domainTransformer(transformation);
    ImageTransformation imageTransformation(transformation, domainTransformer, interpolation);

    GenericWriter<Image>::exportFile("output/output.vol", imageTransformation(image));

    return 0;
}
//                                                                           //
///////////////////////////////////////////////////////////////////////////////
