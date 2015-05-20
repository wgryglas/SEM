#ifndef DERIVEDFACTORY_H
#define DERIVEDFACTORY_H

#include <map>
#include <string>

#include "iomanagment/InfoStream.h"

#include "boost/function.hpp"
#include "boost/functional/factory.hpp"

template<typename base, typename derived>
struct ImplementationRegistrar
{
    ImplementationRegistrar(const std::string & derivedName)
    {
        base::CSTORS_MAP()[derivedName]=boost::factory<derived*>();
    }
};

#define DECLARE_IMPLEMENTATION_FACTORY(BASE, ...)                                           \
    public:                                                                                 \
        typedef boost::function<BASE*( __VA_ARGS__)> CSTOR;                                 \
        static  CSTOR Impl(const std::string&);                                             \
    private:                                                                                \
        template<typename Inter, typename Impl>                                             \
        friend ImplementationRegistrar<Inter,Impl>::ImplementationRegistrar(const std::string &); \
        static std::map<std::string,CSTOR>& CSTORS_MAP();

#define DEFINE_IMPLEMENTATION_FACTORY(BASE, ...)                                        \
    std::map<std::string,typename BASE::CSTOR>& BASE::CSTORS_MAP()                      \
    {                                                                                   \
        static std::map<std::string,BASE::CSTOR> constructors;                          \
        return constructors;                                                            \
    }                                                                                   \
    BASE::CSTOR BASE::Impl(const std::string& name)                                     \
    {                                                                                   \
        if(BASE::CSTORS_MAP().find(name) == BASE::CSTORS_MAP().end())                   \
        {                                                                               \
            ErrorInFunction                                                             \
            <<"not known class #BASE implementation to create for specified type name"  \
            <<SEM::iomanagment::endProgram;                                             \
            return NULL;                                                                \
        }                                                                               \
        return BASE::CSTORS_MAP()[name];                                               \
    }


#define REGISTER_IMPLEMENATION(BASE,DERIVED,TYPE_NAME) static ImplementationRegistrar<BASE,DERIVED> registrar##DERIVED(TYPE_NAME);


#endif // DERIVEDFACTORY_H
