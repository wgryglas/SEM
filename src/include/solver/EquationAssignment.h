#ifndef EQUATION_ASSIGNMENT_H
#define EQUATION_ASSIGNMENT_H

namespace SEM { namespace las {
    
    
    /** ******************************************************
     * \class EquationAssigmentBase
     * base class for operator equation operator handling
     * Required for better memory managment -- avoid allocation
     * of any temporar equation matrix and vector.
     * EquationAssigmentBase role is to perform assigment.
     * Eny of funcator which implements some discretization
     * should not use direct assigment between matrix and
     * and it's coefficient, but shall use "assign" method,
     * which will provide apropriate operator defined in
     * full equation.
     **********************************************************/
    template<typename Derived>
    struct AssigmentBase
    {
        template<typename From, typename To>
        inline void operator()(const From & from, To & to) const
        {
            static_cast<const Derived&>(*this).assign(from,to);
        }
    };
    
    #define MAKE_EQ_ASSIGMENT_IMPLEMENTATION(OPERATOR, NAME) \
    struct NAME : public AssigmentBase<NAME>                 \
    {                                                        \
        template<typename From, typename To>                 \
        inline void assign(const From & from, To & to) const \
        {                                                    \
            to OPERATOR from;                                \
        }                                                    \
    };                                                       \

MAKE_EQ_ASSIGMENT_IMPLEMENTATION( =,EqAssigment)
MAKE_EQ_ASSIGMENT_IMPLEMENTATION(+=,AddEqAssigment)
MAKE_EQ_ASSIGMENT_IMPLEMENTATION(-=,SubstractEqAssigment)
MAKE_EQ_ASSIGMENT_IMPLEMENTATION(*=,MultiplyEqAssigment)
MAKE_EQ_ASSIGMENT_IMPLEMENTATION(/=,DivideEqAssigment)
    
    
}//las
}//SEM



#endif //EQUATION_ASSIGNMENT_H