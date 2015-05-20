#ifndef PROPERTY_H
#define PROPERTY_H

#include <string>
#include <map>
#include <stdlib.h>

#include "iomanagment/InfoStream.h"
#include "iomanagment/DictEntry.h"
#include "utilities/ET_Array.h"
#include "utilities/numArray.h"
#include "components/Basic.h"
#include <boost/container/flat_map.hpp>



namespace SEM { namespace materials {
    
    enum PropertyType
    {
        ScalarValue,
        ScalarField,
        VectorValue,
        VectorField,
        TensorValue,
        TensorField
    };
    
    
    
template<typename Derived>
struct PropertyUser
{
    template<typename PType>
    inline void takeProperty(const PType & property) const
    {
        static_cast<const Derived*>(this)->takeProperty(property);
    }
};

class Property
{
    //Scalar types
    Scalar m_scalarValue;
    numArray<Scalar> m_scalarField;
    
    //Vector types
    //TODO
    
    //Tensor tyeps
    //TODO
    
    PropertyType m_type;
    std::string m_name;
public:
    /// \brief Constructor which woudl runtime read it's property value and type from withing some file
    /// input passed by iomanagment::DictEntry
    Property(const std::string &name, const PropertyType& type, const iomanagment::DictEntry& valueFromFile);
    
    /// \brief copy constructor
    Property(const Property& other);
    
    /// \brief Constructor for creating property by passing any single-value property value.
    Property(const std::string &name,  const Scalar& valueToSet );
     
    /// \brief Constructor for passing scalar list like type-->by Expression Templates Array type, which 
    /// is base of all numerical list in this library
    template<typename DerivedExpr>
    Property(const std::string &name, const array::ET_Array<Scalar,DerivedExpr> & expr)
    :m_name(name), m_type(ScalarField), m_scalarField(expr)
    {
    }
    
    /// \TODO IMPLEMENT NEXT PROPETY TYPES, WHEN LIBRARY WOULD MANAGE TO
    /// FULLY SUPPORT VECTOR AND TENSOR TYPES.
    
    /// \brief Destructor
    virtual ~Property();
    
    /// \brief assigment operator
    Property& operator=(const Property& other);
    
    /// \brief insertProperty method for inserting
    /// runtime selected property type inside 
    /// class object which implements PropertyUser
    /// interface.
    /// \param user - object where shall be this property 
    /// applied
    /// \param expectedSize=-1, optional parameter, which
    /// shall be passed if we may expect some sort of array
    /// as property, which shall have desired size. Wrong
    /// size between provided to method and one read from
    /// file/provided by constructor would cause runtime 
    /// error.
    template<typename DerivedUser>
    inline void applyProperty(const PropertyUser<DerivedUser> & user, size_t expectedSize=0)
    {
        using namespace iomanagment;
        
        switch(m_type)
        {
            case ScalarValue:
                user.takeProperty(m_scalarValue);
                break;
            case ScalarField:
                using namespace iomanagment;
                if(expectedSize>0)
                    if(m_scalarField.size()==1)
                    {
                        m_scalarField.resize(expectedSize,m_scalarField[0]);
                    }
                    else
                    {
                        ErrorInFunction<<"Can't use provided in file size of ScalarField property type\n"
                                       <<"when program expects diffrent size of this field. Possible \n"
                                       <<"wrong definition of this property in appeared in user-file.\n"
                                       <<"Check apropriate file(or other source of property), for propetry\n"
                                       <<"named "<<m_name<<endProgram;
                    }
                
                user.takeProperty(m_scalarField);
                break;
            default:
            {
                ErrorInFunction<<"Runtime error,\n"
                                 "not supported property type "
                                 <<typeToString(m_type)<<" yet"<<endProgram;
            }
        }
    }
    
    typedef std::map<std::string,PropertyType> typesMap;
    static typesMap Type;
    
    
private:
    static std::string typeToString(const PropertyType &type);
    
};



}//materials
}//SEM

#endif // PROPERTY_H
