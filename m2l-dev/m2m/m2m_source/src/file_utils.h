/******************************************************************************
* Functions to work with files.
*
* Author: Vitaliy Ogarko, vogarko@gmail.com
*******************************************************************************/

#ifndef SRC_FILEUTILS_H_
#define SRC_FILEUTILS_H_

#if defined(_WIN32)
    #include <direct.h>
#else
    #include <sys/stat.h>
    #include <sys/types.h>
#endif

#include <fstream>

namespace FileUtils {

// Create a directory.
inline void CreateDirectory(const char* path)
{
#if defined(_WIN32)
    _mkdir(path);
#else
    mkdir(path, 0777); // notice that 777 is different than 0777
#endif
}

// Makes a copy of a file.
void CopyFile(const char* source, const char* destination)
{
    std::ifstream file_src(source, std::ios::binary);
    std::ofstream file_dst(destination, std::ios::binary);

    file_dst << file_src.rdbuf();
}
}

#endif
