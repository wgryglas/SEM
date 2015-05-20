#ifndef READ_WRITE_STREAM_H
#define READ_WRITE_STREAM_H

#include <iostream>
#include <string>
#include <list>
#include <vector>
#include <valarray>

#include "boost/array.hpp"
#include "iomanagment/InfoStream.h"
/////////////////////////////////////////////////////
/// \file ReadWriteStream.h
/// General notes for read/write functions:
/// --------------------------------------------------
/// When user would like to use read/write method on
/// own type, then he should support operators '<<'
/// and '>>' in own class.
/// --------------------------------------------------
/// When user type is some sort of List-container
/// class then there is no need to write reading method.
/// Below templates shall do the work, until container
/// supports operator [], and "size" (and "resize" for
/// dynamicly changing size containers).
/// --------------------------------------------------
/// It is assumed that containers with <T,size_t> template
/// argumetns are fixed size, and conatiners wit <T> arg
/// are dynamicly changed.
/// --------------------------------------------------
/// Inner container type and container type itself is
/// resolved by read(write)Helper struct specialization
/// --------------------------------------------------
/// Standard containers(like std::vector, std::list,
/// std::valaray, boost::array) shall work fine
/// with read/write functions
/// --------------------------------------------------
/// Additionaly there is introdueced other type, which
/// aim is to support reading into Lists (list wrappers)
/// derived(wraped) from complex template type container.
/// It's done here by ContainerType provider template
/// argument. This prvider shall be simple structure
/// which makes typedef "type" depending on one template
/// parameters, eg. Make List with 1 template Type arg.
/// from SpecialList<T,A,B,C>:
/// -
/// SpecialList<T,a,b,c>--> desired container to be derived
/// VectorProvider would be like for it like below:
/// -
/// template<T> struct VectorProvider
/// {typedef SpecialList<T,a,b,c> type;};
/// -
/// In fact it may be a little be complex, because we
/// could ommit Provider and just dircetly derive into
/// single template arg. class and use standard readers.
/// described above this note.
/// But in this approach, with provider, our new List
/// is not limitted to one SpecialList instance, but
/// to as many as we would make providers for them
///////////////////////////////////////////////////////

namespace SEM { namespace iomanagment {

//////////////////////////////////////////////////////////////////////////////////////////////////
// WRITING INTO STREAM
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////
/// function to write basic data into stream
/// including containers and special string form
//////////////////////////////////////////////////
template<typename T>
std::ostream& write(const T& val, std::ostream& stream);

//////////////////////////////////////////////////
/// \struct writeHelper
/// ----------------------------------------------
/// workaround to allow writing list/single value
/// by calling one template function
/////////////////////////////////////////////////
template<typename T>
struct writeHelper
{
    static void writeValue(const T & val, std::ostream & s)
	{
		s<<val;
	}
};

//////////////////////////////////////////////////
/// \struct writeHelper<std::string>
/// ----------------------------------------------
/// specialization for writing string in quotted
/// form
/////////////////////////////////////////////////
template<>
struct writeHelper<std::string>
{
    static void writeValue(const std::string& val, std::ostream& s)
    {
        s<<'"'<<val<<'"';
    }
};


////////////////////////////////////////////
/// container write/read workaround
/// --------------------------------
/// workaround for writing into any list
/// Inner type is resolved by readHelper
/// sepcialization (done inside macro)
/// ----------------------------------
/// note: Container type require access
/// operator and size method
////////////////////////////////////////////
template<typename Inner, typename Cont>
struct listWriter
{
    static void wrtieContainter(const Cont & list, std::ostream& stream, bool writeSize=true)
    {
        if(list.size()==0)
        {
            stream<<"0()";
            return;
        }

        // Check if container is uniform
        // (lower file, and if is not uniform,
        // then ususaly it's very fast found,
        // and time spend on checking is gained
        // on writing time, so this solution is
        // quite good
        bool uniform=true;
        const Inner& first=list[0];
        for(int i=1;i<list.size(); ++i)
        {
            if(list[i]!=first)
            {
                uniform=false;
                break;
            }
        }

        //Write size
        if(writeSize)
        {
            stream<<list.size()<<' ';
        }

        //Write in appropriate form
        if(uniform)
        {
            stream<<'[';
            write<Inner>(first, stream);
            stream<<']';
        }
        else
        {
            stream<<'('<<'\n';
            for(int i=0;i<list.size();i++)
            {
                write<Inner>(list[i], stream);
                stream<<'\n';
            }
            stream<<')';
        }

    }
};

//////////////////////////////////////////////////////
/// \struct Templete spacialization for writeHelper
/// for fixed sieze containers
/// --------------------------------------------------
/// It resloves template types in compilation time
/// after just calling write<AnyCompoundType>(in,out)
//////////////////////////////////////////////////////
template<typename T, size_t N, template<class,size_t> class C>
struct writeHelper<C<T,N> >
{
    static void writeValue(const C<T,N>& list, std::ostream& stream)
    {
        listWriter<T,C<T,N> >::wrtieContainter(list,stream,false);
    }
};

//////////////////////////////////////////////////////
/// \struct Templete spacialization for writeHelper
/// for dynamic change size containers (one Template
/// argument, eg. std::valarray)
/// --------------------------------------------------
/// It resloves template types in compilation time
/// after just calling write<AnyCompoundType>(in,out)
//////////////////////////////////////////////////////
template<typename T,template<class> class C>
struct writeHelper<C<T> >
{
    static void writeValue(const C<T>& list, std::ostream& stream)
    {
        listWriter<T,C<T> >::wrtieContainter(list,stream);
    }
};

//////////////////////////////////////////////////////
/// \struct Templete spacialization for writeHelper
/// for dynamic change size containers (two Template
/// argument, eg. std::vector, std::list)
/// --------------------------------------------------
/// It resloves template types in compilation time
/// after just calling write<AnyCompoundType>(in,out)
//////////////////////////////////////////////////////
template<typename T, typename Alloc,template<class, class> class C>
struct writeHelper<C<T,Alloc> >
{
    static void writeValue(const C<T,Alloc>& list, std::ostream& stream)
    {
        listWriter<T,C<T,Alloc> >::wrtieContainter(list,stream);
    }
};


//////////////////////////////////////////////////////
/// \struct Templete spacialization for writeHelper
/// for wrappers of dynamic change size containers.
/// Those wrappers can use default containers or any more
/// complex template container as wrapped data.
/// --------------------------------------------------
/// It resloves template types in compilation time
/// after just calling write<AnyCompoundType>(in,out)
/// --------------------------------------------------
/// Wrapped container type is provided by C_Provider type,
/// which must define typedef C_Provider::type-->see
/// utilities/ListWrapp.h-->this specialization is
/// aimed for those types.
/// --------------------------------------------------
/// Wrapper shall support operator [] and "size"
/// method to work correctyly with this writer.
/// --------------------------------------------------
/// This writer can also be used for types that derives
/// (not only wrappes) from complex container type.
/// That only neccessery thing is template structure:
/// *
/// templat<class T, template<class> class W>
/// class DerivedContainer
/// {
///  ...
/// };
/// *
/// Where:
/// (T-element type, W-Container type provider)
/// *
/// (Container provider introduced for keeping any list like
/// continer as general as it is possible - Provider choses
/// from complex template only those instances which will
/// look like list)
/// --------------------------------------------------
/// For derived classes from complex template, which
/// use only one(or 2 with size) template type
/// (eg.:template<T> class{..};) compilator will chose
/// to use fixed/dynamic size readers described above
//////////////////////////////////////////////////////
template<typename T,template<class> class C_Provider, template<class,template<class> class > class W>
struct writeHelper<W<T,C_Provider> >
{
    static void writeValue(const W<T,C_Provider >& list, std::ostream& stream)
    {
        listWriter<T,W<T,C_Provider> >::wrtieContainter(list,stream);
    }
};

//////////////////////////////////////////////////
/// write function definition
//////////////////////////////////////////////////
template<typename T>
std::ostream& write(const T& val, std::ostream& stream)
{
    writeHelper<T>::writeValue(val, stream);
    return stream;
}


//////////////////////////////////////////////////////////////////////////////////////////////////
// READING FROM STREAM
//////////////////////////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////
/// function to read basic data from stream
/// including containers and special string form
//////////////////////////////////////////////////
template<typename T>
std::istream & read(std::istream& stream, T& val );

//////////////////////////////////////////////////
/// \struct readHelper
/// wrokaround to allow reading in the same manner
/// single/container data and special string form
//////////////////////////////////////////////////
template<typename T> struct readHelper
{
    static void readValue(std::istream & in, T& out)
    {
        in>>out;
    }
};

//////////////////////////////////////////////////
/// \struct readHelper<std::string>
/// specialization for string quotted form
//////////////////////////////////////////////////
template<> struct readHelper<std::string>
{
    static void readValue(std::istream & in, std::string & out)
    {
        using namespace SEM::iomanagment;

        char c;
        if(!in.get(c))
        {
            ErrorInFunction<<"Can't read string value - stream is empty"<<endProgram;
        }
        
        while(std::iswspace(c))
        {
            if(!in.get(c))
            {
                ErrorInFunction<<"Can't read string value - stream is empty"<<endProgram;
            }
        }
        

        if(c!='"')
        {
            ErrorInFunction<<"Wrong string definition - expected \" char for string begining"<<endProgram;
        }

        do
        {
            if(!in.get(c))
            {
                ErrorInFunction<<"While readign string value found unexpected end of stream before \" char was reached"<<endProgram;
            }

            if(c=='"')
            {
                break;
            }
            else
            {
                out.push_back(c);
            }
        }
        while(true);
    }
};

//////////////////////////////////////////////////
/// \struct fixedListReadHelper
/// ---------------------------------------------
/// workaround for reading into fixe size list.
/// Inner and Container (note: its any container,
/// with any number of templ. arguments) type is
/// resolved by readHelper.
/// sepcialization is done inside below struct
/// ---------------------------------------------
/// fixedListReadHelper - struct allows easy
/// reading into any container, which suppor basic
/// operation -access operator [] and method "size"
//////////////////////////////////////////////////
template<typename Inner, typename Cont> struct fixedListReader
{
    static void readContainer( std::istream & in, Cont & out )
    {
        using namespace std;

        in>>skipws;
        char c;
        in>>c;

        if(std::isdigit(c)) //Unneccessary token for fixed list. Supported only when size is equal
        {
            in.unget();//move back to read full number
            int size;
            in>>size;

            if(size!=out.size())
            {
                ErrorInFunction<<"Wrong size specified in stream for reading fixed container. \n"
                               <<"Fixed container know its size before runtime read, so possible \n"
                               <<"you should remove size token from place before list or change it \n"
                               <<"to match expected size(size token is supported only for consistency \n"
                               <<"with dynamic size list) \n"
                               <<"If you are sure that here shall be placed size, then consider changing \n"
                               <<"in source code container to the one that can dynamicly change size \n"
                               <<endProgram;
            }

            //Go to '(' or '[' char
            in>>c;
        }


        if(c!='(' && c!='[')
            ErrorInFunction<<"Unexpected char \'"<<c<<"\' found while reading list, expected '(' or ']' for uniform list"<<endProgram;

        readFixedSize_StartTokensChecked(in, out, c);
    }

    /// extracted part of reading, because this method will be used
    /// by dynamicListReadHelper in case when size is known
    /// --> one method call more, but kept basic rule not to duplicate code.
    /// \param c - opening list char
    static void readFixedSize_StartTokensChecked( std::istream & in, Cont & out, char c )
    {
        if(c=='[')//uniform list
        {
            if(out.size()>0)
                read<Inner>(in,out[0]);
            
            in>>c;
            if(c!=']')
                ErrorInFunction<<"Expected ']' at end of uniform list, found "<<c<<" token"<<endProgram;

            //Write list
            for(int i=1; i<out.size(); ++i)
                out[i]=out[0];

        }
        else //nonuniform list
        {
            for(int i=0; i<out.size(); ++i)
            {
//                 Inner val;
                read<Inner>(in,out[i]);
//                 out[i]=val;
            }

            in>>c;
            if(c!=')')
                ErrorInFunction<<"Expected ')' at end of nonuniform list, found "<<c<<" token"<<endProgram;
        }
    }
};

//////////////////////////////////////////////////
/// \struct Template specialization for fixed
/// size contianers reading.
/// ---------------------------------------------
/// After resolving approprite types T,N,C struct
/// delegates operation into general container type
/// reader-fixedListReadHelper
/// ---------------------------------------------
/// Container type must support access operator []
/// and "size" method.
//////////////////////////////////////////////////
template<typename T, size_t N, template<class,size_t> class C>
struct readHelper<C<T, N> >
{
    static void readValue(std::istream & in, C<T, N>& out)
    {
        fixedListReader<T,C<T,N> >::readContainer(in,out);
    }
};

//////////////////////////////////////////////////
/// \struct dynamicListReadHelper
/// ---------------------------------------------
/// workaround for reading into dynamicly
/// allcated list.
/// Inner and Container (note: its any container,
/// with any number of templ. arguments) type is
/// resolved by readHelper.
/// sepcialization is done inside below struct
/// macro
/// ---------------------------------------------
/// dynamicListReadHelper - struct allows easy
/// reading into any container, which suppor basic
/// operation -access operator "[]",
/// methods "size" and "resize".
/// ---------------------------------------------
/// note: if container size is set, then
/// program will try to read only this
/// size of input (unless stream specifies
/// diffrent size)
/// If size is equal to "0", and in stream
/// there is no size specifying digit, then
/// memory allocation is not performed
//////////////////////////////////////////////////
template<typename Inner, typename Cont>
struct dynamicListReader
{
    static void readContainer(std::istream & in, Cont & out )
    {
        using namespace std;

        in>>skipws;
        char c;
        in>>c;

        if(std::isdigit(c))
        {
            in.unget();//move back to read full number
            int size;
            in>>size;
            //Allocate memory
            out.resize(size);

            //Go to '(' char
            in>>c;
        }

        if(c!='(' && c!='[')
            ErrorInFunction<<"Unexpected char \'"<<c<<"\' found while reading list, expected '(' or '[' for uniform list"<<endProgram;

        if(out.size()==0)//Dynamic allocation basing on chars (...)
        {
            if(c!='(')
                ErrorInFunction<<"Wrong token found, expected '(' for list without specified size. \n"
                               <<"'[' is reserevd for uniform list, so size must be know befor parsing. \n"
                               <<"Write down size before '[' token in file, or change source code  \n"
                               <<"to one that will allocate memeory in list before stream is parsed"<<endProgram;

            do
            {
                Inner val;
                read<Inner>(in,val);

                out.resize(out.size()+1);
                out[out.size()-1]=val;

                //out.push_back(val);

                in>>c;//check next char if it's not end of list

                if(!in.good())
                    ErrorInFunction<<"Unexpected end of stream before ')' token was reachd \n"
                                   <<"while reading list from stream with size auto deducated"<<endProgram;

                if(c!=')')
                    in.unget();
                else
                    break;
            }
            while(true);
        }
        else //Memory is allocated read basing on list size (like fixed size container)
        {
            fixedListReader<Inner,Cont>::readFixedSize_StartTokensChecked(in, out, c);
//            if(c=='{')//uniform list
//            {
//                Inner val;
//                read<Inner>(in,val);
//                in>>c;
//                if(c!='}')
//                    ErrorInFunction<<"Expected '}' at end of uniform list, found "<<c<<" token"<<endProgram;

//                //Write list
//                for(int i=0; i<out.size(); ++i)
//                    out[i]=val;

//            }
//            else //nonuniform list
//            {
//                Inner val;
//                for(int i=0; i<out.size(); ++i)
//                {
//                    read<Inner>(in,val);
//                    out[i]=val;
//                }

//                in>>c;
//                if(c!=')')
//                    ErrorInFunction<<"Expected ')' at end of nonuniform list, found "<<c<<" token"<<endProgram;
//            }
        }
    }
};

//////////////////////////////////////////////////
/// \struct Template specialization for dynamic
/// size contianers reading (one template arg,
/// eg. std::valarray).
/// ---------------------------------------------
/// After resolving approprite types T,C struct
/// delegates operation into general container type
/// reader-dynamicListReadHelper
/// ---------------------------------------------
/// Container type must support access operator [],
/// "size" and "resize" methods.
//////////////////////////////////////////////////
template<typename T,template<class> class C>
struct readHelper<C<T> >
{
    static void readValue(std::istream & in, C<T>& out)
    {
        dynamicListReader<T,C<T> >::readContainer(in,out);
    }
};

//////////////////////////////////////////////////
/// \struct Template specialization for dynamic
/// size contianers reading (two template arg,
/// eg. std::vector, std::list).
/// ---------------------------------------------
/// After resolving approprite types T,C struct
/// delegates operation into general container type
/// reader-dynamicListReadHelper
/// ---------------------------------------------
/// Container type must support access operator [],
/// "size" and "resize" methods.
//////////////////////////////////////////////////
template<typename T, typename Alloc, template<class,class> class C>
struct readHelper<C<T,Alloc> >
{
    static void readValue(std::istream & in, C<T,Alloc>& out)
    {
        dynamicListReader<T,C<T,Alloc> >::readContainer(in,out);
    }
};

//////////////////////////////////////////////////////
/// \struct Templete spacialization for readHelper
/// for wrappers of dynamic change size containers.
/// Those wrappers can use default containers or any more
/// complex template container as wrapped data.
/// --------------------------------------------------
/// It resloves template types in compilation time
/// after just calling read<AnyCompoundType>(in,out)
/// --------------------------------------------------
/// Wrapped container type is provided by C_Provider type,
/// which must define typedef C_Provider::type-->see
/// utilities/ListWrapp.h-->this specialization is
/// aimed for those types.
/// --------------------------------------------------
/// Wrapper shall support operator [], "size" and
/// "resize" methods to work correctyly with this
/// reader.
/// --------------------------------------------------
/// This reader can also be used for types that derives
/// (not only wrapps) from complex container type.
/// That only neccessery thing is template structure:
/// *
/// templat<class T, template<class> class W>
/// class DerivedContainer
/// {
///  ...
/// };
/// *
/// Where:
/// (T-element type, W-Container type provider)
/// *
/// (Container provider introduced for keeping any list like
/// continer as general as it is possible - Provider choses
/// from complex template only those instances which will
/// look like list)
/// --------------------------------------------------
/// For derived classes from complex template, which
/// use only one(or 2 with size) template type
/// (eg.:template<T> class{..};) compilator will chose
/// to use fixed/dynamic size readers described above
//////////////////////////////////////////////////////
template<typename T,template<class> class C_Provider, template<class,template<class> class> class W>
struct readHelper< W<T, C_Provider> >
{
    static void readValue(std::istream & in, W<T,C_Provider>& out)
    {
        dynamicListReader<T,W<T,C_Provider> >::readContainer(in,out);
    }
};

//////////////////////////////////////////////////
/// read function definition
//////////////////////////////////////////////////
template<typename T>
std::istream & read(std::istream& stream, T& val )
{
    readHelper<T>::readValue(stream, val);
}

}//iomanagment
}//SEM




#endif// READ_WRITE_STREAM_H

