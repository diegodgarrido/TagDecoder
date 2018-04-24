#include "fileutils.h"
#include <stdarg.h>
#include <algorithm>
#include <set>


//All of the hard work
#if 0
void lsfiles(string folder, vector<string> &files)
{
    //vector<string> files;  // Will be added to List
    vector<string> folders; // Will be added to List
    char search_path[200];
    sprintf_s(search_path, "%s*.jpg", folder.c_str());
    WIN32_FIND_DATA fd;
    HANDLE hFind = ::FindFirstFile(search_path, &fd);
    if (hFind != INVALID_HANDLE_VALUE)
    {
        do
        {
            // read all (real) files in current folder, delete '!' read other 2 default folder . and ..
            if (!(fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
            {

                files.push_back(fd.cFileName);
            }
            else //Put folders into vector
            {
                folders.push_back(fd.cFileName);
            }
        }
        while (::FindNextFile(hFind, &fd));
        ::FindClose(hFind);
    }

    //return (files);
}
#endif
bool removeStringFromString(std::string& str, const std::string& from, const std::string& to)
{
    size_t start_pos = str.find(from);
    if (start_pos == std::string::npos)
        return false;
    str.replace(start_pos, from.length(), to);
    return true;
}



int getFrameSuffixes(const char *dirname, const string suffix_name, vector<string> &suffix_names)
{

    DIR *dir;
    std::set<std::string> sortedItems;
    struct dirent *ent;
    if ((dir = opendir(dirname)) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL)
        {
            char buffer[200];
#ifdef __GNUC__
            strcpy(buffer, (char *)ent->d_name);
#else
            strcpy_s(buffer, (char *)ent->d_name);
#endif
            string full_name(buffer);
            size_t found = full_name.find(suffix_name);
            if (found != std::string::npos)
            {
                sortedItems.insert(full_name);
            }
        }
        // Iterate through all the elements in a set and store the value.
        for (std::set<std::string>::iterator it=sortedItems.begin(); it!=sortedItems.end(); ++it)
            suffix_names.push_back(*it);
        closedir(dir);
    }
    else
    {
        /* could not open directory */
        perror("");
        return EXIT_FAILURE;
    }
    return(1);
}

int getFramePrefixes(const char *dirname, const string prefix_name,vector<string> &prefix_names)
{
    DIR *dir;
    std::set<std::string> sortedItems;
    struct dirent *ent;
    if ((dir = opendir(dirname)) != NULL)
    {
        /* print all the files and directories within directory */
        while ((ent = readdir(dir)) != NULL)
        {
            char buffer[200];
#ifdef __GNUC__
            strcpy(buffer, (char *)ent->d_name);
#else
            strcpy_s(buffer, (char *)ent->d_name);
#endif
            string full_name(buffer);
            size_t found = full_name.find(prefix_name);
            if (found != std::string::npos)
            {
                sortedItems.insert(full_name);
            }
            //printf("%s\n", buffer);
        }
        // Iterate through all the elements in a set and store the value.
        for (std::set<std::string>::iterator it = sortedItems.begin(); it != sortedItems.end(); ++it)
            prefix_names.push_back(*it);
        closedir(dir);

    }
    else
    {
        /* could not open directory */
        perror("");
        return EXIT_FAILURE;
    }
    return(1);
}

bool getFrameNames(const string prefix_name, const vector<string> &prefix_names,string index_name_begin,string index_name_end,vector<string> &frame_names)
{
    int idx0 = atoi(index_name_begin.c_str());
    int idx1 = atoi(index_name_end.c_str());
    size_t count = 0;
    for (int idx = idx0; idx <= idx1; idx++)
    {
        string idx_name(prefix_name);
        idx_name += to_string(idx) + '-';
        size_t found = prefix_names[count].find(prefix_name);
        if (found != std::string::npos)
        {
            frame_names.push_back(prefix_names[count]);
        }
        else
        {
            return(false);
        }
        count++;
    }
    return(true);
}



