///////////////////////////////////////////////////////////////////////////////
// standard
#include <cmath>
#include <string>

// DGtal
#include "DGtal/io/writers/GenericWriter.h"
#include "DGtal/io/readers/GenericReader.h"
#include "DGtal/helpers/StdDefs.h"
#include "DGtal/images/ImageContainerBySTLVector.h"

// Gudhi
#include <gudhi/Bottleneck.h>

//! [include]
#include "AffineTransformation.h"
#include "Interpolation.h"
#include "ImageTransformation.h"
#include "Topology.h"
#include "Homology.h"
#include "utils/CLI11.hpp"
#include "utils/gnuplot-iostream.h"



/**
 * @file evalTopoTrans2D.cpp
 * @author Sarah Brood (\c sarahbrood\@ecole.ensicaen.fr)
 * @author Heithem Dridi (\c heithemdridi\@ecole.ensicaen.fr)
 *
 * @date 2022/01
 *
 * This file is part of the project : Topological evaluation of image transformations
 * The aim of this file is to show an example of use of our project for 2D
 */

///////////////////////////////////////////////////////////////////////////////

using namespace std;
using namespace DGtal;
using namespace functors;
using namespace Z2i;

///////////////////////////////////////////////////////////////////////////////


typedef unsigned char Value;
typedef ImageContainerBySTLVector<Domain, Value> Image_;

template<typename TInterpolation> 
int process(Image_ input, const string& CLI_output, bool homology){

    /// Transformation

    typedef AffineTransformation2D<Space> AffineTransformation;
    typedef DomainTransformer2D<Domain, AffineTransformation> DomainTransformer;
    typedef ImageTransformation<Space, Image_, AffineTransformation, DomainTransformer, TInterpolation> ImageTransformation;
    
    AffineTransformation transformation = AffineTransformation::rotation(RealPoint(1, 0), M_PI_4);
    DomainTransformer domainTransformer(transformation);
    TInterpolation interpolation(&input);
    ImageTransformation imageTransformation(transformation, domainTransformer, interpolation);
    Image_ output = imageTransformation(input);
    GenericWriter<Image_>::exportFile(CLI_output, output);

    // Topology
    
    typedef Topology2D<Image_, KSpace, DigitalSet, DT4_8 , Object4_8> Topology2D;
    typedef PersistentHomology<Image_, DT4_8, Topology2D> PersistentHomology;
    typedef typename PersistentHomology::Diagram Diagram;
    
    if(homology){
        PersistentHomology phIN(input, dt4_8, true);
        Diagram diagIN = phIN();
        PersistentHomology phOUT(output, dt4_8, true);
        Diagram diagOUT = phOUT();

        double b0 = Gudhi::persistence_diagram::bottleneck_distance(diagIN[0], diagOUT[0]);
        double b1 = Gudhi::persistence_diagram::bottleneck_distance(diagIN[1], diagOUT[1]);

        cout << "Bottleneck distance for b0 = " << b0 << endl;
        cout << "Bottleneck distance for b1 = " << b1 << endl;

        Gnuplot gpb0;
        gpb0 << "set xrange [-10:300]\n"
             << "set yrange [-10:300]\n"
             << "set key outside\n";
        gpb0 << "plot x,"
             << gpb0.file1d(diagIN[0]) << "with points pointtype 1 title 'b0 input',"
             << gpb0.file1d(diagOUT[0]) << "with points pointtype 4 title 'b0 output'"
             << endl;

        Gnuplot gpb1;
        gpb1 << "set xrange [-10:300]\n"
             << "set yrange [-10:300]\n"
             << "set key outside\n";
        gpb1 << "plot x,"
             << gpb1.file1d(diagIN[1]) << "with points pointtype 1 title 'b1 input',"
             << gpb1.file1d(diagOUT[1]) << "with points pointtype 4 title 'b1 output'"
             << endl;
    }else{
        Image_ inputExtended = PersistentHomology::extend(input);
        Image_ inputReversed = PersistentHomology::reverse(inputExtended);
        Topology2D inputTopo(&inputExtended, &inputReversed, dt4_8);
        for(int i : inputTopo(1, 255)){
            cout << i << " ";
        }
        cout << endl;
        Board2D inputBoard;
        inputTopo.toBoard(inputBoard);
        inputBoard.saveCairo("output/board_input.png");

        Image_ outputExtended = PersistentHomology::extend(output);
        Image_ outputReversed = PersistentHomology::reverse(outputExtended);
        Topology2D outputTopo(&outputExtended, &outputReversed, dt4_8);
        for(int i : outputTopo(1, 255)){
            cout << i << " ";
        }
        cout << endl;
        Board2D outputBoard;
        outputTopo.toBoard(outputBoard);
        outputBoard.saveCairo("output/board_output.png");
    }

    return 0;
}


int main(int argc, char** argv){

    // parse command line using CLI --------------------------------------------
    CLI::App app;
    app.description("This file is part of the project : Topological evaluation of image transformations\n"
                    "The aim of this file is to show an example of use of our project for 2D");
    string CLI_input;
    string CLI_output = "output/output.pgm";
    string CLI_interpolation = "NN";
    bool CLI_homology = true;

    app.add_option("-i,--input,1", CLI_input, "input file in .pgm / .ppm" )
            ->required()
            ->check(CLI::ExistingFile);

    app.add_option("-o,--output,2", CLI_output, "output file in .pgm / .ppm. Default to output/output.pgm" );

    app.add_option("-I,--interpolation, 3", CLI_interpolation,
                   "set interpolation: "
                    "NN: Nearest Neighbour, "
                    "BIL: Bilinear, "
                    "BIC: Bicubic")
            -> check(CLI::IsMember({"NN", "BIL", "BIC"}));

    app.add_option("-H,--homology,4", CLI_homology, "if true compute homology. Default to true" );


    app.get_formatter()->column_width(40);
    CLI11_PARSE(app, argc, argv);
    // END parse command line using CLI ----------------------------------------

    Image_ input = GenericReader<Image_>::import(CLI_input);

    typedef NearestNeighbour<Space, Image_> NN;
    typedef Bilinear<Space, Image_> BIL;
    typedef Bicubic<Space, Image_> BIC;

    if(CLI_interpolation == "BIL"){
        return process<BIL>(input, CLI_output, CLI_homology);
    }
    if(CLI_interpolation == "BIC"){
        return process<BIC>(input, CLI_output, CLI_homology);
    }

    return process<NN>(input, CLI_output, CLI_homology);
}