#include "iomanagment/case.h"
#include "time/Time.h"
#include "mesh/Mesh.h"
#include "postprocessing/PostObject.h"
#include "iomanagment/InfoStream.h"

int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::mesh;
    
    //Setup case
    Case::setup(argc, args);
    
    //Create mesh
    Mesh mesh(Case::meshPath(),Case::elementsPath());

    
    iomanagment::Dictionary postDict;
    Case::postprocessingObjectsDir() >> postDict;
    
    //Create post objects defined in file "postprocessing.sem"
    std::vector<std::auto_ptr<PostObject> > postObjects;
    for(auto objDictPair : postDict.subDictionaries())
    {
        postObjects.push_back(std::auto_ptr<PostObject>
        (
            PostObject::Impl(objDictPair.second->entry("postCalculator").value())(*objDictPair.second, mesh, Case::time() )
        ));
    }


    // Get times which need to be evaluated
    std::vector<double> times =Case::timesInWorkingDir();
    if(postDict.hasEntry("time"))
    {
        std::string timeType = postDict.entry("time").value();
        if( timeType == "first" && times.size()>0)
        {
            times = {times[0]};
        }
        else if(timeType == "last" && times.size() > 0)
        {
            times = {times[times.size()-1]};
        }
        else if(timeType != "all")
        {
            ErrorInFunction << "not known option \"time "<<timeType<<"\"to do post calculation \n"
                            << "allowed: first,last,all"<<endProgram;
        }
    }
    else
    {
        WARNING(std::cout)<<"time for calculation is not set, calculation performed for all time results"<<endInfo;
    }
     
     
     
    //Iterate over times and perform calculation 
    for(double time: times)
    {
        Case::time().setCurrentTime(time);
        Info <<  "calculating postprocessing objects for time "<<time<<endInfo;
        
        for(auto & postObjPtr : postObjects)
        {
            postObjPtr->calculate();
        }
    }
}