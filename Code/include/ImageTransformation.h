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
 * @file images/ImageTransformation.h
 * @author Sarah Brood & Heithem Dridi

 * @date 2021/10/13
 *
 */
#include <utility>
#include "DGtal/base/Common.h"
#include "DGtal/kernel/BasicPointFunctors.h"
#include <DGtal/helpers/StdDefs.h>
#include <DGtal/kernel/domains/CDomain.h>
#include <DGtal/kernel/CSpace.h>

///////////////////////////////////////////////////////////////////////////////
namespace DGtal{
    namespace functors{
        template<typename TSpace, typename TImage, typename TAffineTransform, typename TDomainTransformer,
                typename TInterpolation>
        class ImageTransformation{
            BOOST_CONCEPT_ASSERT((concepts::CSpace<TSpace>));
        public:
            typedef typename TImage::Domain Domain;
            typedef typename TSpace::RealPoint RealPoint;
            typedef typename TDomainTransformer::Bounds Bounds;
        public:
            ImageTransformation(const TAffineTransform &transform, const TDomainTransformer &domainTransformer,
                                const TInterpolation &interpolation): _transform(transform),
                                _domainTransformer(domainTransformer), _interpolation(interpolation){}

            TImage operator()(const TImage &input) const{
                typedef typename Domain::ConstIterator It;

                trace.beginBlock("Starting process");
                Bounds bounds = _domainTransformer(input.domain());
                Domain transformedDomain(bounds.first, bounds.second);
                TImage output(transformedDomain);

                trace.info() << "Image: " << output << std::endl;
                RealPoint transformed;
                bool inside;
                for (It it = output.domain().begin(); it != output.domain().end(); ++it){
                    inside = true;
                    transformed = _transform(*it);

                    for(unsigned int i = 0; i < RealPoint::dimension; i++){
                        if(transformed[i] < input.domain().lowerBound()[i] ||
                           transformed[i] > input.domain().upperBound()[i]){
                            inside = false;
                            break;
                        }
                    }

                    if(inside){
                        output.setValue(*it, _interpolation(transformed));
                    }
                }

                trace.endBlock();
                return output;
            }

        protected:
            TAffineTransform _transform;
            TDomainTransformer _domainTransformer;
            TInterpolation _interpolation;
        };

    }
}
