#ifndef LINEARINTERPOLATOR_H
#define LINEARINTERPOLATOR_H

#include <vector>
#include <array>
#include <sstream>

#include "components/Basic.h"
#include "components/CmpTraits.h"
#include "iomanagment/InfoStream.h"

namespace SEM {
    

template<typename P, typename T>
struct InterpolatorHelper;

template<typename T>
struct InterpolatorHelper<Scalar,T>
{
    static std::pair<size_t, size_t> findPair(const Scalar &pos, const std::vector<std::pair<Scalar,T> > & table)
    {
        for(size_t i=0; i<table.size()-1; ++i)
        {
            if(pos >= table[i].first && pos < table[i+1].first)
            {
                return std::make_pair(i, i+1);
            }
        }
        
        if(table[0].first < pos)
        {
            return std::make_pair(table.size()-1,table.size()-1);
        }
        else
        {
            return std::make_pair(0, 0);
        }
    }

    static T doInterpolation(const Scalar & pos, const std::pair<size_t,size_t> & bounds, const std::vector<std::pair<Scalar,T> > & table)
    {
        Scalar c = (pos-table[bounds.first].first)/(table[bounds.second].first-table[bounds.first].first);
        return (table[bounds.second].second - table[bounds.first].second)*c + table[bounds.first].second;
    }
};

template<typename T>
struct InterpolatorHelper<Vector,T>
{
    static std::pair<size_t, size_t> findPair(const Vector &pos, const std::vector<std::pair<Vector,T> > & table)
    {
        Vector tmp1, tmp2, dir;
        
        for(size_t i=0; i<table.size()-1; ++i)
        {
            for(size_t d=0; d<pos.size(); ++d)
            {
                tmp1[d] = table[i].first[d] - pos[d];
                tmp2[d] = table[i+1].first[d] - pos[d];
                dir[d]  = table[i+1].first[d] - table[i].first[d];
            }
            
            if( dotProd(tmp1, dir) * dotProd(tmp2,dir) <= 0)
                return std::make_pair(i, i+1);
        }
        
        for(size_t d=0; d<pos.size(); ++d)
        {
            tmp1[d] = table[1].first[d] - pos[d];
            dir[d]  = table[1].first[d] - table[0].first[d];
        }
        
        if(dotProd(tmp1,dir) < 0)
            return std::make_pair(0,0);
        else
            return std::make_pair(table.size()-1, table.size()-1);
        
    }
    
    static T doInterpolation(const Vector & pos, const std::pair<size_t,size_t> & bounds, const std::vector<std::pair<Vector,T> > & table)
    {
        Vector tmp1, dir;
        
        for(size_t d=0; d<pos.size(); ++d)
        {
            tmp1[d] = pos[d]-table[bounds.first].first[d];
            dir[d]  = table[bounds.second].first[d]-table[bounds.first].first[d];
        }
        dir.normalize();
        
        Scalar c = dotProd(tmp1,dir);
        
        return (table[bounds.second].second - table[bounds.first].second)*c + table[bounds.first].second;
    }
    
};


    
    
template<typename Pos, typename T>
struct LinearInterpolator
{
   LinearInterpolator()
   {
   }
    
   LinearInterpolator(const std::vector<std::pair<Pos,T> > & values)
   : m_values(values)
   {
   }
   
   LinearInterpolator(const LinearInterpolator<Pos, T> & other)
   : m_values(other.m_values)
   {
   }
   
   void update(const std::vector<std::pair<Pos,T> > & values)
   {
       using namespace SEM::iomanagment;
       
       if(values.size()<2)
           ErrorInFunction<<"table for interpolation need at least 2 entries"<<endProgram;
       
       m_values = values;
   }
   
   void update(const iomanagment::DictEntry & valueEntry)
   {
       std::vector<std::array<Scalar,CmpTraits<Pos>::DIM_SIZE + CmpTraits<T>::DIM_SIZE> > tableInput;
       valueEntry>>tableInput;
       
       size_t pDim = CmpTraits<Pos>::DIM_SIZE;
       size_t vDim = CmpTraits<T>::DIM_SIZE;
       
       m_values.clear();
       m_values.reserve(tableInput.size());
       for(size_t i=0; i<tableInput.size(); ++i)
       {
           T val;
           Pos pos;
           
           for(size_t d=0; d<pDim; ++d)
                CmpTraits<Pos>::component(pos,d)=tableInput[i][d];
           
           for(size_t d=0; d< vDim; ++d)
               CmpTraits<T>::component(val,d)=tableInput[i][d+pDim];
           
           m_values.push_back(std::pair<Pos,T>(pos,val));
       }
   }
   
   void writeValues(iomanagment::DictEntry & valueEntry) const
   {
       std::stringstream out;
       out << "(\n";
       for(unsigned int i=0; i<m_values.size(); ++i)
       {
           
            out<< "( ";
            
            for(unsigned int d=0; d< CmpTraits<Pos>::dim(); ++d)
            {
                out << ' ' << CmpTraits<Pos>::component(m_values[i].first,d);
            }
            
            for(unsigned int d=0; d< CmpTraits<T>::dim(); ++d)
            {
                out << ' ' << CmpTraits<T>::component(m_values[i].second,d);
            }
           
           out << ")\n";
       }
       out << ")\n";
       
       valueEntry.setStrValue(out.str());
   }
    
    T operator()(const Pos & current) const
    {
//         unsigned int begin =-1;
//         for(unsigned int i=0; i<m_values.size()-1; ++i)
//         {
//             if(current >= m_values[i].first && current < m_values[i+1].first)
//             {
//                 begin = i;
//             }
//         }
//         
//         if(begin ==-1)
//         {
//             if(m_values[0].first < current)
//                 return m_values.back().second;
//             else
//                 return m_values[0].second;
//         }

        
//         Pos c = (current-m_values[begin].first)/(m_values[begin+1].first-m_values[begin].first);
//         return (m_values[begin+1].second -m_values[begin].second)*c + m_values[begin].second;

        std::pair<size_t,size_t> bounds = InterpolatorHelper<Pos,T>::findPair(current,m_values);
        
        if(bounds.first == bounds.second )
            return m_values[bounds.first].second;

        return InterpolatorHelper<Pos,T>::doInterpolation(current,bounds, m_values);
    }
  
  private:
      std::vector<std::pair<Pos,T> > m_values;
};  


} //SEM

#endif // LINEARINTERPOLATOR_H
