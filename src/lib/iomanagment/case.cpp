#include "case.h"

#include <cstdlib>
#include <string>

#include "InfoStream.h"

namespace SEM
{

std::string Case::ENV_VARAIBLE_INSTALL_DIR_NAME = "SEM_INSTALL_DIR";
boost::filesystem::path Case::INSTALL_DIR="";
boost::filesystem::path Case::CASE_PATH="";
SolutionControl *Case::m_solControl = NULL;
Time *Case::m_time = NULL;
materials::Material *Case::m_material = NULL;

void Case::setup( int size,char* args[])
{
    using boost::filesystem::path;
    using boost::filesystem::initial_path;

    path wd;
    if(size==1)
    {
        wd = initial_path();
    }
    else
    {
        wd = path(args[1]);
        if(wd.is_relative())
        {
            wd = initial_path();
            wd.append<char*>(args[1], args[1]+sizeof(args[1]));
        }
    }
    
    Case::CASE_PATH = wd;
    Case::INSTALL_DIR = getSEMInstallDirEnvVariable();
    

    Case::setupTime();
}

void Case::setup( const boost::filesystem::path &workingPath, const boost::filesystem::path &installDir)
{
    if(installDir.empty())
    {
        Case::INSTALL_DIR = getSEMInstallDirEnvVariable();
    }
    else
    {
        Case::INSTALL_DIR = installDir;
    }
    
    Case::CASE_PATH=workingPath;

    Case::setupTime();
}


boost::filesystem::path Case::path()
{
    return Case::CASE_PATH;
}

boost::filesystem::path Case::meshDir()
{
    return Case::CASE_PATH/"mesh";
}

boost::filesystem::path Case::meshPath()
{
    return Case::meshDir()/"mesh.sem";
}

boost::filesystem::path Case::elementsPath()
{
    return Case::meshDir()/"elements.sem";;
}

boost::filesystem::path Case::nodesPath()
{
    return Case::meshDir()/"spectralNodes.sem";
}

boost::filesystem::path Case::materialPath()
{
    return Case::CASE_PATH/"material.sem";
}

boost::filesystem::path Case::materialDataBaseDir()
{
    return Case::INSTALL_DIR/"etc/materials";
}

boost::filesystem::path Case::postprocessingObjectsDir()
{
    return Case::CASE_PATH/"postprocessing.sem";
}

SolutionControl &Case::solutionControl()
{
    using namespace iomanagment;
    if(!Case::m_solControl)
    {
        Case::m_solControl = new SolutionControl(RegistryFile::ref(
                    new RegistryFile
                            (
                                Case::time(),
                                "solution",
                                boost::filesystem::path(),
                                READ_IF_MODIFIED,NO_WRITE
                            )
                ));
    }

    return *Case::m_solControl;
}

Time &Case::time()
{
    return *Case::m_time;
}

materials::Material& Case::material()
{
    if(!m_material)
    {
        using namespace iomanagment;
        RegistryFile::ref matFile
        (
            new RegistryFile
            (
                Case::time(), 
                "material", 
                boost::filesystem::path(),
                READ_IF_MODIFIED,NO_WRITE
            ) 
        );
        m_material = new materials::Material(matFile,Case::materialDataBaseDir());
    }
    
    return *m_material;
}

std::vector< double > Case::timesInWorkingDir() 
{
    std::vector<double> timesInDir;
    
    boost::filesystem::directory_iterator end;
    boost::filesystem::directory_iterator itr(CASE_PATH);
    
    for(;itr!=end; ++itr)
    {
        std::stringstream ss(itr->path().filename().string());
        double val = -1;
        
        if(ss >> val)
        {
            if(val >= 0. && boost::filesystem::is_directory(*itr))
            {
                timesInDir.push_back(val);
            }
        }
    }
    
    std::sort(timesInDir.begin(),timesInDir.end());
    
    return timesInDir;
    
}

void Case::setupTime()
{
    
    // Create new time object
    if(Case::m_time)
        delete Case::m_time;
    
    Case::m_time = new Time(); //<-- this constructor would manage to connect time object with "control.sem" file
    
    // Override start time if case path contains more then "0" directory
    std::vector<double> timeDirs = timesInWorkingDir();

//     boost::filesystem::directory_iterator end;
//     boost::filesystem::directory_iterator itr(CASE_PATH);    
//     double timeInDir = 0;
//     for(;itr!=end; ++itr)
//     {
//         std::stringstream ss(itr->path().filename().string());
//         double val = -1;
//         
//         if(ss >> val)
//         {
//             if(boost::filesystem::is_directory(*itr) && val > timeInDir)
//             {
//                 timeInDir  = val;
//             }
//         }
//     }

    if(timeDirs.size()>0)
    {
        time().setCurrentTime(timeDirs[timeDirs.size()-1]);
    }
}


char* Case::getSEMInstallDirEnvVariable()
{
    
    using namespace iomanagment;
    if( getenv( ENV_VARAIBLE_INSTALL_DIR_NAME.c_str() ) )
        return getenv( ENV_VARAIBLE_INSTALL_DIR_NAME.c_str() );
    else
        ErrorInFunction<<"There is no "<<ENV_VARAIBLE_INSTALL_DIR_NAME<<std::endl
        <<"envirtonment variable in your OS, but this variable"<<std::endl
        <<"is required for SEM programs, please add "<<ENV_VARAIBLE_INSTALL_DIR_NAME<<std::endl
        <<"to your system, and specifie there SEM install directory path"<<endProgram;
}


std::ostream& operator<<(std::ostream& o, Case c)
{
    o<<c.path();
    return o;
}
boost::filesystem::path Case::probeLocationDefinition() 
        {
            return Case::CASE_PATH/"probeLocation.sem";
        }


}


