#ifndef POST_CALCULATION_H
#define POST_CALCULATION_H

#include "components/Basic.h"
#include "utilities/numArray.h"
#include "utilities/VectorUtils.h"
#include "fields/GeometricField.h"

namespace SEM {
 
    namespace field
    {
        template<typename T>
        numArray<T> extrapolateInTime(const field::GeometricField<T> &field,unsigned int order)
        {
            std::vector<Scalar> coeff;
            switch(order)
            {
                case 0: 
                    coeff.push_back(1);
                    break;
                case 1:
                    coeff.push_back(2);
                    coeff.push_back(-1);
                    break;
                case 2:
                    coeff.push_back(3);
                    coeff.push_back(-3);
                    coeff.push_back(1);
                    break;
                default:
                {
                    using namespace SEM::iomanagment;
                    ErrorInFunction << "higher order then 2 in time extrapolation is not supported" <<endProgram;
                }
            }
            
            SEM::numArray<T> result(field.size());
            
            for(int i=0; i< coeff.size(); ++i)    
            {
                result += coeff[i] * field.cachedField(i);
            }
            
            return result;
        }
        
        template<typename T>
        numArray<T> extrapolateInTime(const field::GeometricField<T> &field, const std::vector<int> &mapping,unsigned int order)
        {
            std::vector<Scalar> coeff;
            switch(order)
            {
                case 0: 
                    coeff.push_back(1);
                    break;
                case 1:
                    coeff.push_back(2);
                    coeff.push_back(-1);
                    break;
                case 2:
                    coeff.push_back(3);
                    coeff.push_back(-3);
                    coeff.push_back(1);
                    break;
                default:
                {
                    using namespace SEM::iomanagment;
                    ErrorInFunction << "higher order then 2 in time extrapolation is not supported" <<endProgram;
                }
            }
            
            SEM::numArray<T> result(field.size());
            
            for(int i=0; i< coeff.size(); ++i)    
            {
                result += coeff[i] * field.cachedField(i).slice(mapping);
            }
            
            return result;
        }
    }
    
    
}


#endif // POST_CALCULATION_H
