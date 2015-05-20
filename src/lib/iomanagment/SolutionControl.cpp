#include "SolutionControl.h"

#include "iomanagment/InfoStream.h"
#include "iomanagment/DictEntry.h"

namespace SEM {

SolutionControl::SolutionControl(iomanagment::RegistryFile::ref regFile)
    : RegistryObject("SolutionControl"), m_solverType(BiCGSTAB), m_timeDiscretization(BDF2), m_maxIteration(100), m_tolerance(1e-5)
{
    setRegistryFile(regFile);
}

void SolutionControl::read(const iomanagment::Dictionary &dict)
{
    if(dict.hasEntry("solver"))
        m_solverType = str2SolverType(dict.entry("solver").value());
    
    if(dict.hasEntry("timeDiscretization"))
        m_timeDiscretization = str2TimeDiscretization(dict.entry("timeDiscretization").value() );
    
    if(dict.hasEntry("maxIteration"))
        dict.entry("maxIteration")>>m_maxIteration;
    
    if(dict.hasEntry("tolerance"))
        dict.entry("tolerance")>>m_tolerance;
    
}

void SolutionControl::write(iomanagment::Dictionary &dict)
{
    iomanagment::DictEntry *entry=new iomanagment::DictEntry("solver",solverType2Str(m_solverType));
    dict.add(entry);
    
    entry=new iomanagment::DictEntry("timeDiscretization",timeDiscretization2Str(m_timeDiscretization));
    dict.add(entry);
    
    entry=new iomanagment::DictEntry("maxIteration");
    *entry<<m_maxIteration;
    dict.add(entry);
    
    entry=new iomanagment::DictEntry("tolerance");
    *entry<<m_tolerance;
    dict.add(entry);
}


SolverType SolutionControl::str2SolverType(const std::string &s)
{
    if(s=="ConjugateGradient")
        return ConjugateGradient;
    
    if(s=="SimplicialLLT")
        return SimplicialLLT;
    
    if(s=="SimplicialLDLT")
        return SimplicialLDLT;
    
    if(s=="SimplicialCholesky")
        return SimplicialCholesky;
    
    if( s== "SEM_ConjugateGradient")
        return SEM_ConjugateGradient;
    
    
    return BiCGSTAB;
}

std::string SolutionControl::solverType2Str(const SolverType &type)
{
    switch(type)
    {
    case ConjugateGradient:
        return "ConjugateGradient";
    case SimplicialLLT:
        return "SimplicialLLT";
    case SimplicialLDLT:
        return "SimplicialLDLT";
    case SimplicialCholesky:
        return "SimplicialCholesky";
    case SEM_ConjugateGradient:
        return "SEM_ConjugateGradient";
    default:
        return "BiCGSTAB";
    }
}

TimeDiscretization SolutionControl:: str2TimeDiscretization(const std::string &s)
{
    if(s=="BDF2")
        return BDF2;
    
    if(s=="Euler")
        return Euler;
    
    return Steady;
}

std::string SolutionControl::timeDiscretization2Str(const TimeDiscretization &disc)
{
    if(disc==BDF2)
        return "BDF2";
    
    if(disc==Euler)
       return "Euler";
    
    return "Steady";
}

}//SEM
