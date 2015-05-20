#ifndef _SolutionControl_H_
#define _SolutionControl_H_

#include "iomanagment/RegistryFile.h"
#include "iomanagment/RegistryObject.h"
#include "iomanagment/Dictionary2.h"

namespace SEM {

    enum TimeDiscretization
    {
        Steady=0,
        Euler=1,
        BDF2=2
    };

    enum SolverType
    {
        ConjugateGradient,
        BiCGSTAB,
        SimplicialCholesky,
        SimplicialLLT,
        SimplicialLDLT,
        SEM_ConjugateGradient
    };

    class SolutionControl : public iomanagment::RegistryObject
    {
    public:
        SolutionControl(iomanagment::RegistryFile::ref regFile);

        void read(const iomanagment::Dictionary& dict);
        void write(iomanagment::Dictionary& dict);

        inline SolverType solverType() const {return m_solverType; }
        inline TimeDiscretization timeDiscretization() const { return m_timeDiscretization; }
        inline int maxIteration() const { return m_maxIteration;}
        inline double tolerance() const { return m_tolerance; }

    private:
        SolverType m_solverType;
        TimeDiscretization m_timeDiscretization;
        int m_maxIteration;
        double m_tolerance;

        static SolverType str2SolverType(const std::string& s);
        static std::string solverType2Str(const SolverType & type);

        static TimeDiscretization str2TimeDiscretization(const std::string &s);
        static std::string timeDiscretization2Str(const TimeDiscretization& disc);
    };

} //SEM

#endif //_SolutionControls_H_
