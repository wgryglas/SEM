#ifndef REFERENCE_H
#define REFERENCE_H

#include "boost/shared_ptr.hpp"
#include "boost/weak_ptr.hpp"

/////////////////////////////////////////////////////////
/// \def macro for creating java-like memory
/// managment. Here it's done by boost::shared_ptr,
/// and due to this fact usage of object must be
/// similar to raw pointer.
/// --------------------------------------------------------
/// Usage of this type will be now more eye candy,
/// like SomeClassType::ref --> in contrast to
/// SomeClassType & is not very bad, and autocompleted.
/// --------------------------------------------------------
/// Note: best place to set this macro is at the
/// begining of the class, because by the default
/// class space is "private" and so this macro is
/// ending with "private" keyword.
/// --------------------------------------------------------
/// Note: ref is acctualy a boost::shared_ptr, so
/// when using ref you have to keep in mind that
/// circular references are not automaticaly solved.
/// Boost doc says:
/// "Because the implementation uses reference counting,
///  cycles of shared_ptr instances will not be reclaimed.
///  For example, if main() holds a shared_ptr to A, which
///  directly or indirectly holds a shared_ptr back to A,
///  A's use count will be 2. Destruction of the original
///  shared_ptr will leave A dangling with a use count of 1.
///  Use weak_ptr to "break cycles"."
/////////////////////////////////////////////////////////////
#define REFERENCE_TYPE(CLASS_NAME)                \
public:                                           \
    typedef boost::shared_ptr<CLASS_NAME> ref;    \
    typedef boost::weak_ptr<CLASS_NAME> weekRef;  \
private:                                          \
    
/////////////////////////////////////////////////////
/// \def COMMA - to use inside template CLASS_NAME
/// while user will specify REFERENCE_TYPE
/////////////////////////////////////////////////////
#ifndef COMMA
    #define COMMA ,
#endif
    
/////////////////////////////////////////////////////
/// \def SINGLE_ARG - to use inside template CLASS_NAME
/// while user will specify REFERENCE_TYPE
/////////////////////////////////////////////////////
#ifndef SINGLE_ARG
    #define SINGLE_ARG(...) __VA_ARGS__
#endif

#endif // REFERENCE_H
