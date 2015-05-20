
#ifndef _Time_H_
#define _Time_H_

#include <iostream>
#include <string>

#include "boost/shared_ptr.hpp"

#include "iomanagment/Register.h"
#include "iomanagment/Dictionary2.h"
#include "TimeLocal.h"
#include "TimeChangeListner.h"

namespace SEM
{
    //////////////////////////////////////////////
    /// \class WriteManager
    /// -----------------------------------------
    /// Runtime selection of writing type
    /// This class will control writing registered
    /// objects in Time
    //////////////////////////////////////////////
    class Time;

    class WriteControler
    {
    public:
        virtual bool shallWrite()=0;
        virtual WriteControler* makeCopy(const Time *time)=0;
        virtual void write(iomanagment::Dictionary& dict) const=0;
    };

    //////////////////////////////////////////////
    /// \class Time
    /// -----------------------------------------
    /// Global controler for start/end calculation
    /// Extends Register- then this class is able
    /// to control write/read operations on
    /// registered objects.
    //////////////////////////////////////////////
    class Time : public iomanagment::Register, public iomanagment::RegistryObject
	{
        /// Privet implementations of
        /// WriteControler - one of those types
        /// is chosen at program runtime
        class RuntimeWrite;
        class TimeStepWrite;

        double m_time;
        double m_timeStep;
        double m_timeEnd;
        int m_step;

        double m_lastWriteTime;
        int m_lastWriteStep;

        WriteControler* m_writeControler;

        boost::shared_ptr<iomanagment::LocalPath> m_localPath;

        std::vector<time::TimeChangeListner*> timeListners;

    public:
        Time(double start, double step, double end, int writeInterval=1 );
        Time();
        Time(const Time& other);
        Time& operator =(const Time& other);
        virtual ~Time();

        std::string timeName() const;
        boost::shared_ptr<iomanagment::LocalPath> localPath() const;

        void setCurrentStep(size_t new_step);
        
        void setCurrentTime(double cTime);
        
        bool end() const;
        double timeStep() const;
        double time() const;
        double lastWriteTime() const;

        int step() const;
        int lastWriteStep() const;

        Time& operator++();
        Time operator++(int unused);

        void registerTimeListner(time::TimeChangeListner* listner);
        void unregisterTimeListner(time::TimeChangeListner* listner);

        
        void fireWriting();
        
        /// \brief RegistryObject interface
        void read(const iomanagment::Dictionary &dict);

        /// \brief RegistryObject interface
        void write(iomanagment::Dictionary &dict);


    private:
        inline void callTimeListners();

	};

    std::ostream& operator<<(std::ostream& o,Time& t);


    //---------------- inline definitions -----------------//
    void Time::callTimeListners()
    {
        for(time::TimeChangeListner* listner: timeListners)
        {
            listner->timeChanged(time());
        }
    }

}

#endif //_Time_H_
