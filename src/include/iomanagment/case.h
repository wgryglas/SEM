#ifndef _Case_H_
#define _Case_H_

#include <iostream>
#include <string>

#include "boost/filesystem.hpp"

#include "SolutionControl.h"
#include "time/Time.h"
#include "materials/Material.h"

namespace SEM
{
	class Case
	{
        static std::string ENV_VARAIBLE_INSTALL_DIR_NAME; 
        
        static Time * m_time;
        static SolutionControl *m_solControl;
        static materials::Material *m_material;
        
        static boost::filesystem::path CASE_PATH;
        static boost::filesystem::path INSTALL_DIR;
	public:
        /// \brief setup -  sets path for case from main function arguments, creates time object
        /// \param size  -  main function number of input strings
        /// \param args  -  main function input strigns
        static void setup(int size, char* args[]);

        /// \brief setup - sets path for case, creates time object
        /// \param workingPath - file path to root case folder.
        static void setup( const boost::filesystem::path& workingPath, const boost::filesystem::path& installDir="");

        /// \brief path case path getter
        /// \return patch to root folder
        static boost::filesystem::path path();
        
        static boost::filesystem::path meshDir();
        
        static boost::filesystem::path meshPath();
        
        static boost::filesystem::path elementsPath();
        
        static boost::filesystem::path nodesPath();
        
        static boost::filesystem::path materialPath();
        
        static boost::filesystem::path materialDataBaseDir();
        
        static boost::filesystem::path postprocessingObjectsDir();
        
        static boost::filesystem::path probeLocationDefinition();

        /// \brief solutionControl
        ///        singelton provider of static variable type SolutionControl.
        ///        this variable shall be shard among solvers
        ///        If first time called method, then variable is initialized.
        /// \return the sam instance of SolutionControl-connected with file "control.sem"
        static SolutionControl& solutionControl();

        /// \brief time - getter for shared time object
        /// \return     - instance of case time-shared object
        static Time & time();
        
        
        static materials::Material &  material();
        
        static std::vector<double> timesInWorkingDir();

    private:
        /// Create static variable time - shall be shared by multiple objects
        static void setupTime();
        
        static char* getSEMInstallDirEnvVariable();
    };
    

    std::ostream& operator<<(std::ostream& o, Case );
}

#endif //_Case_H_
