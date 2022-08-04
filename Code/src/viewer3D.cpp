/**
 * @file evalTopoTrans2D.cpp
 * @author Sarah Brood (\c sarahbrood\@ecole.ensicaen.fr)
 * @author Heithem Dridi (\c heithemdridi\@ecole.ensicaen.fr)
 *
 * @date 2022/01
 *
 * This file is part of the project : Topological evaluation of image transformations
 * The aim of this file is to give a tool to view 3D objects
 * It is inspired by https://dgtal-team.github.io/doctools-nightly/Doc3dVolBoundaryViewer.html
 */


#include <iostream>
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/io/viewers/Viewer3D.h"
#include "DGtal/images/imagesSetsUtils/SetFromImage.h"
#include "DGtal/topology/DigitalSurface.h"
#include "DGtal/topology/SetOfSurfels.h"

using namespace std;
using namespace DGtal;
using namespace Z3i;

int main(int argc, char** argv){
    typedef unsigned char Value;
    typedef ImageContainerBySTLVector<Domain,  Value> Image;
    typedef KSpace::SurfelSet SurfelSet;
    typedef functors::IntervalForegroundPredicate<Image> ThresholdImage;
    typedef SetOfSurfels<KSpace, SurfelSet> SetOfSurfels;
    typedef DigitalSurface<SetOfSurfels> DigitalSurface;
    typedef DigitalSurface::ConstIterator ConstIterator;
    typedef SurfelAdjacency<KSpace::dimension> SurfelAdjacency;

    QApplication window(argc,argv);

    trace.beginBlock( "Loading image into memory." );
        string file = argc >= 2 ? argv[1] : "samples/Al100.vol";
        Image image = GenericReader<Image>::import(file);
        trace.info() << "Image loaded: " << image << std::endl;
    trace.endBlock();

    trace.beginBlock( "Construct the Khalimsky space from the image domain." );
        const Domain& domain = image.domain();
        KSpace ks;
        if (!ks.init(domain.lowerBound(), domain.upperBound(), true)){
            trace.error() << "Error in the Khamisky space construction."<<std::endl;
            return 2;
        }
    trace.endBlock();

    trace.beginBlock("Wrapping a digital set around image. ");
        ThresholdImage thresholdImage(image, 0, 255);
    trace.endBlock();

    trace.beginBlock("Extracting boundary by scanning the space. ");
        SurfelAdjacency surfAdj(true); // interior in all directions.
        SetOfSurfels setOfSurfels(ks, surfAdj);
        Surfaces<KSpace>::sMakeBoundary(setOfSurfels.surfelSet(), ks, thresholdImage,
                                        domain.lowerBound(), domain.upperBound());
        DigitalSurface digSurf(setOfSurfels);
        trace.info() << "Digital surface has " << digSurf.size() << "surfels." << std::endl;
    trace.endBlock();

    trace.beginBlock( "Displaying everything. " );
        Viewer3D<Space,KSpace> viewer(ks);
        viewer.setWindowTitle("Simple boundary of volume Viewer");
        viewer.show();
        viewer << SetMode3D(ks.unsigns(*(digSurf.begin())).className(), "Basic");
        for(ConstIterator it = digSurf.begin(), itE = digSurf.end(); it != itE; ++it){
            viewer << ks.unsigns(*it);
        }
        viewer << Viewer3D<>::updateDisplay;
    trace.endBlock();

    return window.exec();
}
//