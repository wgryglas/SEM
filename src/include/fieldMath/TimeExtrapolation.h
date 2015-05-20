
#ifndef TIMEEXTRAPOLATION_H
#define TIMEEXTRAPOLATION_H

#include "utilities/numArray.h"
#include "fields/GeometricField.h"

namespace SEM {

enum ExtrapolationScheme
{
    CurrentValue=0,
    First,
    Second,
    Third
};
    
template<typename T>
numArray<T> extInTime( const field::GeometricField<T> & field, ExtrapolationScheme sheme =Third )
{
//     if( sheme > field.oldFieldsNumber() )
//         sheme = (ExtrapolationScheme) field.oldFieldsNumber();
    
    switch(sheme)
    {
        case Third:
            return 3.*field.oldField(0) -3.*field.oldField(1) + field.oldField(2);
        case Second:
            return 2.*field.oldField(0) -field.oldField(1);
        case First:
            return field.oldField(0);
        default:
            return field;
    }
}

template<typename T>
numArray<T> extInTime( const field::GeometricField<T> & field, const std::vector<int> & map, ExtrapolationScheme sheme =Third )
{
//     if( sheme > field.oldFieldsNumber() )
//         sheme = (ExtrapolationScheme) field.oldFieldsNumber();
    
    switch(sheme)
    {
        case Third:
            return 3.*field.oldField(0).slice(map) - 3.*field.oldField(1).slice(map)  + field.oldField(2).slice(map);
        case Second:
            return 2.*field.oldField(0).slice(map) -field.oldField(1).slice(map);
        case First:
            return field.oldField(0).slice(map);
        default:
            return field.slice(map);
    }
}

} //SEM


#endif // TIMEEXTRAPOLATION_H
