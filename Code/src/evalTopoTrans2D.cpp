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
using namespace Z2i;

///////////////////////////////////////////////////////////////////////////////

int main(int, char**){
    //image
    typedef unsigned char Value;
    typedef ImageContainerBySTLVector<Domain, Value> Image;

    //transformaton
    typedef AffineTransformation2D<Space> AffineTransformation;
    typedef DomainTransformer2D<Image::Domain, AffineTransformation> DomainTransformer;
    typedef Bicubic<Space, Image> Interpolation;
    typedef ImageTransformation<Space, Image, AffineTransformation,
            DomainTransformer, Interpolation> ImageTransformation;

    //open image
    Image image = GenericReader<Image>::import("samples/church.pgm");

    //instanciate
    Interpolation interpolation(&image);
    AffineTransformation transformation = AffineTransformation::rotation(RealPoint(1, 0), M_PI_4);
    DomainTransformer domainTransformer(transformation);
    ImageTransformation imageTransformation(transformation, domainTransformer, interpolation);

    GenericWriter<Image>::exportFile("output/output.pgm", imageTransformation(image));

    return 0;
}