#include "Time.h"

#include <cstdio>
#include <sstream>

#include "iomanagment/InfoStream.h"

namespace SEM
{

//----------------------------------------------------------------//
//          Implementations of WriteControlers                    //
//----------------------------------------------------------------//
    ////////////////////////////////////////////////
    ///     RuntimeWrite implementation
    ////////////////////////////////////////////////
    class Time::RuntimeWrite : public WriteControler
    {
    private:
        /// Dissalow copy/assigment -while copying time must be specified explicitly
        /// because RuntimeWrite is aimed to be used internaly for each instance of Time
        RuntimeWrite(const RuntimeWrite &):m_time(NULL){}
        RuntimeWrite& operator =(const RuntimeWrite &){}

        const Time * m_time;
        double m_interval;
    public:

        RuntimeWrite(const Time *time, const iomanagment::Dictionary& controlDict)
            : m_time(time)
        {
            controlDict.entry("writeInterval")>>m_interval;
        }

        RuntimeWrite(const Time *time, double interval)
            : m_time(time), m_interval(interval)
        {
        }

        ~RuntimeWrite(){}

        WriteControler* makeCopy(const Time *time)
        {
            return new RuntimeWrite(time,m_interval);
        }

        bool shallWrite()
        {
            return (m_time->lastWriteTime()+m_interval) < m_time->time() + m_time->timeStep();
        }

        void write(iomanagment::Dictionary& dict) const
        {
            using iomanagment::DictEntry;

            DictEntry * entry = new DictEntry("writeType");
            entry->setValue("runTime");
            dict.add(entry);

            entry = new DictEntry("writeInterval");
            *entry<<m_interval;
            dict.add(entry);
        }
    };

    ////////////////////////////////////////////////
    ///     TimeStepWrite implementation
    ////////////////////////////////////////////////
    class Time::TimeStepWrite : public WriteControler
    {
    private:
        /// Dissalow copy/assigment -while copying time must be specified explicitly
        /// because RuntimeWrite is aimed to be used internaly for each instance of Time
        TimeStepWrite(const TimeStepWrite &):m_time(NULL){}
        TimeStepWrite& operator =(const TimeStepWrite &){}

        const Time* m_time;
        int m_interval;
    public:

        TimeStepWrite(const Time *time, const iomanagment::Dictionary& controlDict)
            : m_time(time)
        {
            controlDict.entry("writeInterval")>>m_interval;
        }

        TimeStepWrite(const Time *time, int interval)
            : m_time(time), m_interval(interval)
        {
        }

        ~TimeStepWrite(){}

        WriteControler* makeCopy(const Time *time)
        {
            return new TimeStepWrite(time, m_interval);
        }

        bool shallWrite()
        {
            return m_time->lastWriteStep()+m_interval <= m_time->step();
        }

        void write(iomanagment::Dictionary& dict) const
        {
            using iomanagment::DictEntry;

            DictEntry * entry = new DictEntry("writeType");
            entry->setValue("timeStep");
            dict.add(entry);

            entry = new DictEntry("writeInterval");
            *entry<<m_interval;
            dict.add(entry);
        }
    };


//----------------------------------------------------------------//
//                  Time method definitions                       //
//----------------------------------------------------------------//

    Time::Time(double start, double step, double end, int writeInterval )
        : RegistryObject("Time"),
          m_time(start),
          m_step(0),
          m_lastWriteTime(start),
          m_lastWriteStep(0),
          m_timeStep(step),
          m_timeEnd(end),
          m_writeControler(new TimeStepWrite(this,writeInterval)),
          m_localPath(new TimeLocal(*this))
    {
    }

    Time::Time()
    :RegistryObject("Time"),
    m_time(0),
    m_step(0),
    m_lastWriteTime(0),
    m_lastWriteStep(0),
    m_timeStep(1),
    m_timeEnd(1),
    m_localPath(new TimeLocal(*this))
    {
        using namespace iomanagment;

        RegistryFile::ref regFile
                (
                    new RegistryFile
                    (
                        *this,
                        "control",
                        boost::filesystem::path(),
                        READ_IF_MODIFIED,
                        NO_WRITE
                     )
                 );
        setRegistryFile(regFile);
    }

    Time::Time(const Time &other)
        : RegistryObject("Time"),
          m_time(other.m_time),
          m_lastWriteTime(other.m_lastWriteTime),
          m_timeStep(other.m_timeStep),
          m_timeEnd(other.m_timeEnd),
          m_step(other.m_step),
          m_lastWriteStep(other.m_lastWriteStep),
          m_writeControler(other.m_writeControler->makeCopy(this)),
          m_localPath(new TimeLocal(*this))
    {
    }

    Time &Time::operator =(const Time &other)
    {
        m_time           = other.m_time;
        m_lastWriteTime  = other.m_lastWriteTime;
        m_timeStep       = other.m_timeStep;
        m_timeEnd        = other.m_timeEnd;
        m_step           = other.m_step;
        m_lastWriteStep  = other.m_lastWriteStep;
        m_writeControler = other.m_writeControler->makeCopy(this);

        return *this;
    }

    Time::~Time()
    {
        delete m_writeControler;
    }
    
    void Time::setCurrentTime(double cTime) 
    {
//         if(cTime < m_time )
//         {
//             using namespace SEM::iomanagment;
//             Warning << "trying to set time smaller than is already, time change aborted"<<std::endl;
//             return;
//             
//         } else 
            
        if(cTime == m_time)
        {
            return;
        }
        
        if(cTime < m_timeEnd)
        {
            m_step = (cTime-m_time)/m_timeStep +1;
            m_time = cTime;
        }
        else
        {
            m_time = m_timeEnd;
            m_step = (m_timeEnd - m_time) / m_timeStep +1;
        }
        
        fireReading();
        
        callTimeListners();
        
        using namespace iomanagment;
        Info << "Time::"<<timeName()<<endInfo;
    }
    
    void Time::setCurrentStep(size_t new_step) 
    {
        setCurrentTime(time() + (new_step - step() )*timeStep());
    }
    
    Time& Time::operator++()
    {
        m_time+=m_timeStep;
        ++m_step;

        fireReading();
        
        callTimeListners();
        
        using namespace iomanagment;
        Info << "Time::"<<timeName()<<endInfo;

        return *this;
    }

    Time Time::operator++(int unused)
    {
        Time prevTime(*this);
        operator ++();
        return prevTime;
    }

    void Time::fireWriting() 
    {
        if(m_writeControler->shallWrite())
        {
            iomanagment::Register::fireWriting();
            m_lastWriteTime = m_time;
            m_lastWriteStep = m_step;
        }
    }
    
    void Time::registerTimeListner(time::TimeChangeListner *listner)
    {
        timeListners.push_back(listner);
    }

    void Time::unregisterTimeListner(time::TimeChangeListner *listner)
    {
        std::vector<time::TimeChangeListner*>::iterator itr = std::find(timeListners.begin(), timeListners.end(), listner);
        if(itr!=timeListners.end())
            timeListners.erase(itr);
    }

    void Time::read(const iomanagment::Dictionary &dict)
    {
        if(m_time == 0) //change time only if shure that time is not yet running
        {
            dict.entry("timeStart")>>m_time;
            m_lastWriteTime = m_time; //set last write time to current--> avoid revriting intial setup
        }
        
        dict.entry("timeStep")>>m_timeStep;
        dict.entry("timeEnd")>>m_timeEnd;

        iomanagment::DictEntry& writeType = dict.entry("writeType");

        if(writeType=="runTime")
        {
            m_writeControler = new RuntimeWrite(this,dict);
        }
        else if(writeType=="timeStep")
        {
            m_writeControler = new TimeStepWrite(this,dict);
        }
        else
        {
            ErrorInFunction<<"Not known writing type "<<writeType.value()<<std::endl
                           <<"Allowed: runTime | timeStep"<<iomanagment::endProgram;
        }
    }

    void Time::write(iomanagment::Dictionary &dict)
    {
        using namespace iomanagment;

        m_writeControler->write(dict);

        DictEntry * entry = new DictEntry("timeStart");
        *entry<<m_time;
        dict.add(entry);

        entry=new DictEntry("timeStep");
        *entry<<m_timeStep;
        dict.add(entry);

        entry=new DictEntry("v");
        *entry<<m_timeEnd;
        dict.add(entry);
    }

    std::string Time::timeName() const
    {
        std::ostringstream stream;
        stream<<m_time;
        return stream.str();
    }

    boost::shared_ptr<iomanagment::LocalPath> Time::localPath() const
    {
        return m_localPath;
    }

    double Time::time() const
    {
        return m_time;
    }

    double Time::lastWriteTime() const
    {
        return m_lastWriteTime;
    }

    int Time::step() const
    {
        return m_step;
    }

    int Time::lastWriteStep() const
    {
        return m_lastWriteStep;
    }


    bool Time::end() const
    {
        return m_time >= m_timeEnd;
    }

    double Time::timeStep() const
    {
        return m_timeStep;
    }

    std::ostream& operator<<(std::ostream& o,Time& t)
    {
        return o<<t.timeName();
    }
    

}
