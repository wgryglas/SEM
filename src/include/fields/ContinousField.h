#ifndef CONTINOUSFIELD_H
#define CONTINOUSFIELD_H

#include <iostream>

#include "utilities/numArray.h"
#include "iomanagment/Dictionary2.h"
#include "iomanagment/DictEntry.h"
#include "iomanagment/RegistryObject.h"
#include "utilities/TypeDefs.h"

namespace SEM { namespace field {

///////////////////////////////////////////////////////////////////
/// \class ContinousField<T>
/// --------------------------------------------------------------
/// Extends numArray<T>
/// --------------------------------------------------------------
/// Allows to use advantage of expression templates
/// when calculating any kinde of expr.
/// Implements RegistryObject class to allow automatic
/// read and write.
/// --------------------------------------------------------------
/// operations like +-*/ are element-wise, so each operation is
/// done for each element. Operations thanks to ET are evaluated
/// in only one loop.
/// --------------------------------------------------------------
/// well known std::... functions are supported when creating
/// expression, eg. you can write:
/// field = ( filed1 + sin(field2)*abs(field3) )/5;
/// --------------------------------------------------------------
/// at first glance this class had to be derived from
/// NumericalListWrap, but when solution to expr. templates was
/// found, then there was no further need to use this class. It is
/// suggested that, you shall not use any more any ListWrap classes,
/// it shall be deprecated.
///////////////////////////////////////////////////////////////////
template <typename T>
class ContinousField : public numArray<T>, public iomanagment::RegistryObject
{
    REFERENCE_TYPE(ContinousField<T>)

protected:
    typedef numArray<T> baseListType;

public:
    using iomanagment::RegistryObject::setRegistryFile;

    /// \brief ContinousField - constructor with registration this object in file
    /// \param regFile - registry fiele, object which will call method read/write
    /// \param size - start size of this field, default 0
    ContinousField(iomanagment::RegistryFile::ref regFile, int size=0)
    : RegistryObject(regFile->fileName()),baseListType(size)
    {
        setRegistryFile(regFile);
    }

    /// \brief ContinousField - defualt constructor
    /// \param size -size of field, defualt 0
    ContinousField(int size=0, const std::string &objName="Field"): RegistryObject(objName), baseListType(size)
    {
    }

    /// \brief ContinousField copy-constructor
    ///            !!! - when filed is copied then field is not assigned to
    ///                regFile. It shall be constructed by other constructor,
    ///                and then by "=" operator asigned.
    /// \param other
    ContinousField(const ContinousField<T>& f)
    : RegistryObject("Field"), baseListType(f)
    {
    }
    
    /// \brief ContinousField copy-constructor
    ///            !!! - when filed is copied then field is not assigned to
    ///                regFile. It shall be constructed by other constructor,
    ///                and then by "=" operator asigned.
    /// \param other
    ContinousField(const baseListType& f)
    : RegistryObject("Field"), baseListType(f)
    {
    }
    

    /// \brief operator = -assigment operator (in fact,like numArray is unneccessery, but shuld be faster)
    /// \param f - other field
    /// \return - reference to this class object
    ContinousField<T> & operator = (const ContinousField<T> &f)
    {
        baseListType::operator =(f);
        return *this;
    }
    
    /// \brief operator = -assigment operator (in fact,like numArray is unneccessery, but shuld be faster)
    /// \param f - other field
    /// \return - reference to this class object
    ContinousField<T> & operator = (const baseListType &f)
    {
        baseListType::operator =(f);
        return *this;
    }

    /// \brief ContinousField -empty destructor
    virtual ~ContinousField()
    {
    }

    /// \brief read - implementation of RegistryObejct read method,
    ///               which is called when registryFile decide
    /// \param dict - dictionary read from file
    virtual void read(const iomanagment::Dictionary &dict)
    {
        // entry of dictionary under >> calls iomanagment::read
        dict.entry("value")>>*this;
    }

    /// \brief read - implementation of RegistryObejct read method,
    ///               which is called when registryFile decide
    /// \param dict - dictionary read from file
    virtual void write(iomanagment::Dictionary &dict)
    {
        iomanagment::DictEntry* entry = new iomanagment::DictEntry("value");
        // entry of dictionary under << calls iomanagment::write
        *entry<<*this;
        dict.add(entry);
    }


    // ----------------------------------------------------------------//
    //      EXPRESION TEMPLATES HANDLING
    // ----------------------------------------------------------------//
    /// copy constructor from expresion
    template<typename Exp>
    ContinousField(const array::ET_Array<T,Exp> & exp): RegistryObject("Field"),baseListType(exp)
    {
    }

    /// asigment operators from expressions
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE(T, SINGLE_ARG(ContinousField<T>),baseListType)

};


} //field
} //SEM

#endif // CONTINOUSFIELD_H
