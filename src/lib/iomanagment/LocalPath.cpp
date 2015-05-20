#include "LocalPath.h"

namespace SEM { namespace iomanagment {

LocalPath::LocalPath()
{
}

LocalPath::LocalPath(const std::string &path): m_path(path)
{
}

LocalPath::LocalPath(const boost::filesystem::path &path): m_path(path)
{
}

boost::filesystem::path LocalPath::path() const
{
    return m_path;
}




}//iomanagment
}//SEM
