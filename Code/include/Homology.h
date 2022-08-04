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
 * @file Homology.h
 * @author Sarah Brood (\c sarahbrood\@ecole.ensicaen.fr)
 * @author Heithem Dridi (\c heithemdridi\@ecole.ensicaen.fr)
 *
 * @date 2022/01
 *
 * This file is part of the project : Topological evaluation of image transformations
 */
#pragma once

#if defined(Homology_RECURSES)
#error Recursive header files inclusion detected in Homology.h
#else
#define Homology_RECURSES

#if !defined Homology
#define Homology

namespace DGtal{
    /**
     * Description of class 'PersistentHomology' <p>
    * \brief Aim: compute persistent homology of an image. It can be 2D or 3D
    * @code

    typedef Topology2D<Image, KSpace, DigitalSet, DT4_8 , Object4_8> Topology2D;
    typedef PersistentHomology<Image, DT4_8, Topology2D> PersistentHomology;
    typedef PersistentHomology::Diagram Diagram;
    PersistentHomology phIN(input, dt4_8, true);
    Diagram diagIN = phIN();

     * @tparam TImage
     * @tparam TDigitalTopology Topology adjacency
     * @tparam TTopology Topology class 2D or 3D
     */
    template<typename TImage, typename TDigitalTopology, typename TTopology>
    class PersistentHomology{
    public:
        typedef std::pair<int, int> Pair;
        typedef std::vector<Pair> Persistence;
        typedef std::map<int, Persistence> Diagram;

        /**
         * @brief Constructor
         * @param image
         * @param dt digital topology to use (adjacency for foreground and background)
         * @param optimize if true, compute less value
         */
        PersistentHomology(TImage image, TDigitalTopology dt, bool optimize = false)
        : _image(extend(image)), _reverse(reverse(_image)), _topo(&_image, &_reverse, dt){
            typename TImage::ConstRange range = _image.constRange();
            int val;
            for(unsigned char value : range){
                val = (int) value;
                if(optimize){
                    int min = (val / 10) * 10;
                    int max = min + 10;
                    val = (val - min > max - val)? max : min;
                }

                if(val != 255){
                    _values.insert(val ? val : 1);
                }
            }
        }

        /**
         * @brief Compute topological characteristics in grayscale range
         * @return Diagram of birth/death of each Betti numsbers
         */
        Diagram operator()(){
            trace.beginBlock("Starting homology");
            Diagram diagram;
            std::vector<int> old;
            for(int d = 0; d < TImage::dimension; d++){
                diagram[d] = Persistence();
                old.push_back(0);
            }
            std::vector<int> current;
            Persistence::iterator younger;
            std::set<int>::iterator value = _values.begin();
            while(value != _values.end()){
                current = _topo(*value, 255);
                std::cout << *value << std::endl;
                for(int d = 0; d < TImage::dimension; d++){
                    if(current[d] > old[d]){
                        for(int i = 0; i < current[d] - old[d]; i++){
                            diagram[d].push_back(std::make_pair(*value, 255));
                        }
                    }else if(current[d] < old[d]){
                        for(int i = 0; i < old[d] - current[d]; i++){
                            younger = std::max_element(diagram[d].begin(), diagram[d].end());
                            diagram[d].push_back(std::make_pair(younger->first, *value));
                            diagram[d].erase(younger);
                        }
                    }
                }
                old = current;
                value++;
            }
            trace.endBlock();
            return diagram;
        }

        static TImage reverse(const TImage &image){
            typedef typename TImage::Value Value;
            TImage reversed(image);
            typename TImage::Iterator it;

            for (it = reversed.begin(); it != reversed.end(); ++it){
                *it = (Value) 255 - (int) *it;
            }

            return reversed;
        }

        static TImage extend(const TImage &image){
            typedef typename TImage::Domain::ConstIterator It;
            typename TImage::Domain domain(image.domain().lowerBound() - TImage::Point::diagonal(1),
                                           image.domain().upperBound() + TImage::Point::diagonal(1));
            TImage extended(domain);

            for (It it = domain.begin(); it != domain.end(); ++it){
                if(image.domain().isInside(*it)){
                    extended.setValue(*it, image(*it));
                }
            }

            return extended;
        }
    protected:
        TImage _image;
        TImage _reverse;
        TTopology _topo;
        std::set<int> _values;

    };

}

#endif

#undef Homology_RECURSES
#endif