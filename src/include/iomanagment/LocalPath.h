#ifndef LOCALPATH_H
#define LOCALPATH_H

#include "boost/filesystem/path.hpp"

#include "utilities/Reference.h"

namespace SEM { namespace iomanagment {

class LocalPath
{
    REFERENCE_TYPE(LocalPath)

protected:
    boost::filesystem::path m_path;

public:
    LocalPath();
    LocalPath(const std::string & path);
    LocalPath(const boost::filesystem::path& path);

    virtual boost::filesystem::path path() const;
};

}//iomanagment
}//SEM

#endif // LOCALPATH_H
