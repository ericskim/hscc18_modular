#ifndef FILEHANDLER_HH_
#define FILEHANDLER_HH_

#include <iostream>
#include <sstream>
#include <fstream>
#include <string>
#include <vector>
#include <limits>

#include "Maps.hh"

#define SCOTS_FH_VERSION    1.3
#define SCOTS_FH_SYMBOL     "#"
#define SCOTS_FH_SEPERATOR  ";"

#define SCOTS_FH_KEY        "#SCOTSFILE_FH_VERSION :"
#define SCOTS_FH_TYPE       "#TYPE :"
#define SCOTS_FH_MEMBER     "#MEMBER :"
#define SCOTS_FH_DETMAP     "#DETERMINISTICMAP: "
#define SCOTS_FH_NDETMAP    "#NON-DETERMINISTICMAP: "
#define SCOTS_FH_VECTOR     "#VECTOR :"
#define SCOTS_FH_ARRAY      "#ARRAY :"
#define SCOTS_FH_BEGIN      "#BEGIN :"
#define SCOTS_FH_END        "#END"

namespace scots {

class FileHandler{

public:
    FileHandler();
    FileHandler(std::string);
    FileHandler(const FileHandler&);
    ~FileHandler();

    std::string getFileName();

protected:
    std::string m_filename;
};

/*!
 * \brief The FileWriter class is used to write information into files
 */
class FileWriter :public FileHandler
{
public:
    FileWriter();
    FileWriter(std::string);
    FileWriter(const FileWriter&);
    ~FileWriter();

    bool create(std::string = "");  //!< tries to open a file erasing its previous content
    bool open();                    //!< tries to open the current file keeping its previous content
    void close();                   //!< closes the current file

    bool add_VERSION(void);
    bool add_TYPE(std::string);

    template<class T>
    bool add_MEMBER(std::string,T&);
    template<class T>
    bool add_VECTOR(std::string,std::vector<T>&);
    template<class T>
    bool add_ARRAY(std::string,T&,size_t);

    bool add_deterministicMap(deterministicMap&,std::string = "");
    bool add_nonDeterministicMap(nonDeterministicMap&,std::string = "");

private:
    std::ofstream   m_file;
};

/*!
 * \brief The FileReader class is used to read information from files
 */
class FileReader :public FileHandler
{
public:
    FileReader();
    FileReader(std::string);
    FileReader(const FileReader&);
    ~FileReader();

    bool open();    //!< tries to open the file
    void close();   //!< closes the file

    size_t get_VERSION(double&, size_t = 0);
    size_t get_TYPE(std::string &, size_t = 0);

    template<class T>
    size_t get_MEMBER(std::string,T&,size_t = 0);
    template<class T>
    size_t get_VECTOR(std::string,std::vector<T>&,size_t = 0);
    template<class T>
    size_t get_ARRAY(std::string,T&,size_t,size_t = 0);

    size_t get_deterministicMap(deterministicMap&,std::string ="",size_t = 0);
    size_t get_nonDeterministicMap(nonDeterministicMap&,std::string ="",size_t = 0);

private:
    std::ifstream       m_file;
    std::string         m_line;

    bool skipOffset(size_t&);
    void backToFirstLine();
};

FileHandler::FileHandler() :m_filename("")
{
}

FileHandler::FileHandler(std::string filename) : m_filename(filename)
{
}

FileHandler::FileHandler(const FileHandler &other) :m_filename(other.m_filename)
{
}

FileHandler::~FileHandler()
{
}

std::string FileHandler::getFileName()
{
    return m_filename;
}


FileWriter::FileWriter() : FileHandler::FileHandler()
{
}

FileWriter::FileWriter(std::string filename) :FileHandler::FileHandler(filename)
{

}

FileWriter::FileWriter(const FileWriter& other) :FileHandler::FileHandler(other)
{
}

FileWriter::~FileWriter()
{
}

bool FileWriter::create(std::string filename)
{
    m_file.close();

    if(filename != "")
    {
        m_filename = filename;
    }

    m_file.open(m_filename,std::fstream::out);

    return m_file.is_open();
}

bool FileWriter::open()
{
    m_file.close();

    m_file.open(m_filename,std::fstream::app);

    return m_file.is_open();
}

void FileWriter::close()
{
    m_file.close();
}

bool FileWriter::add_VERSION()
{
    if(m_file.is_open())
    {
        m_file << SCOTS_FH_KEY << SCOTS_FH_VERSION <<std::endl;

        return true;
    }
    else
    {
        return false;
    }
}

bool FileWriter::add_TYPE(std::string type)
{
    if(m_file.is_open())
    {
        m_file << SCOTS_FH_TYPE << type << std::endl;
        return true;
    }
    else
    {
        return false;
    }
}

bool FileWriter::add_deterministicMap(deterministicMap& map, std::string name)
{
    if(m_file.is_open())
    {
        m_file << SCOTS_FH_DETMAP << name << std::endl;
        m_file << map.NrOfCells() << SCOTS_FH_SEPERATOR << map.NrOfInputIndices() << SCOTS_FH_SEPERATOR << map.NrOfValidPairs() << std::endl;
        m_file << SCOTS_FH_BEGIN << std::endl;

        for(size_t cell = 0; cell < map.NrOfCells(); cell ++)
        {
            if(map[cell] != std::numeric_limits<abs_type>::max())
            {
                m_file << cell << SCOTS_FH_SEPERATOR << map[cell] << std::endl;
            }
        }
        m_file << SCOTS_FH_END << std::endl;

        return true;
    }
    else
    {
        return false;
    }

}

bool FileWriter::add_nonDeterministicMap(nonDeterministicMap& map, std::string name)
{
    if(m_file.is_open())
    {
        m_file << SCOTS_FH_NDETMAP << name << std::endl;
        m_file << map.NrOfCells() << SCOTS_FH_SEPERATOR << map.NrOfInputIndices() << SCOTS_FH_SEPERATOR << map.NrOfValidPairs() << std::endl;
        m_file << SCOTS_FH_BEGIN << std::endl;

        booleanMatrix minimalValidPairs = map.extractMinimalValidPairs();

        for(size_t cell = 0; cell < map.NrOfCells(); cell ++)
        {
            if(minimalValidPairs.get(cell,0))
            {
                m_file << cell;

                for(size_t index = 0; index < map.NrOfInputIndices();index++)
                {
                    if(map.get(cell,index))
                    {
                        m_file << SCOTS_FH_SEPERATOR << index;
                    }
                }

                m_file << std::endl;
            }
        }

        m_file << SCOTS_FH_END << std::endl;

        return true;
    }
    else
    {
        return false;
    }
}
template<class T>
bool FileWriter::add_MEMBER(std::string name ,T &member)
{
    if(m_file.is_open())
    {
        m_file << SCOTS_FH_MEMBER << name << std::endl;
        m_file << member << std::endl;

        return true;
    }
    else
    {
        return false;
    }
}
template<class T>
bool FileWriter::add_VECTOR(std::string name, std::vector<T> & vector)
{
    if(m_file.is_open())
    {
        m_file << SCOTS_FH_VECTOR << name << std::endl;
        m_file << SCOTS_FH_BEGIN << vector.size() << std::endl;
        for(size_t index = 0; index < vector.size(); index++)
        {
            m_file << vector.at(index) << std::endl;
        }
        m_file << SCOTS_FH_END << std::endl;

        return true;

    }
    else
    {
        return false;
    }
}
template<class T>
bool FileWriter::add_ARRAY(std::string name, T& array, size_t arraySize)
{
    if(m_file.is_open())
    {
        m_file << SCOTS_FH_ARRAY << name << std::endl;
        m_file << SCOTS_FH_BEGIN << arraySize << std::endl;
        for(size_t index = 0; index < arraySize; index++)
        {
            m_file << array[index] << std::endl;
        }
        m_file << SCOTS_FH_END << std::endl;

        return true;

    }
    else
    {
        return false;
    }
}

FileReader::FileReader() :FileHandler::FileHandler()
{
}

FileReader::FileReader(std::string filename) :FileHandler::FileHandler(filename)
{
}

FileReader::FileReader(const FileReader &other) :FileHandler::FileHandler(other)
{
}

FileReader::~FileReader()
{
}

bool FileReader::open()
{
    m_file.open(m_filename);

    return m_file.good();
}

void FileReader::close()
{
    m_file.close();
}

size_t FileReader::get_VERSION(double& version, size_t offset)
{
    backToFirstLine();
    if(skipOffset(offset))
    {
        size_t counter =0;

        while(std::getline(m_file,m_line))
        {
            counter++;
            if(m_line.find(SCOTS_FH_KEY)!=std::string::npos)
            {
                std::istringstream stream(m_line.substr(m_line.find(":")+1));
                stream >> version;

                return counter;
            }
        }

        return 0;

    }
    else
    {
        return 0;
    }
}

size_t FileReader::get_TYPE(std::string& string, size_t offset)
{
    backToFirstLine();
    if(skipOffset(offset))
    {
        string.clear();
        size_t counter =0;

        while(std::getline(m_file,m_line))
        {
            counter++;
            if(m_line.find(SCOTS_FH_TYPE)!=std::string::npos)
            {
                std::istringstream stream(m_line.substr(m_line.find(":")+1));
                stream >> string;

                return counter;
            }
        }

        return 0;

    }
    else
    {
        return 0;
    }
}

size_t FileReader::get_deterministicMap(deterministicMap& map, std::string name, size_t offset)
{
    backToFirstLine();
    if(skipOffset(offset))
    {
        size_t counter = 0;
        std::string match = SCOTS_FH_DETMAP + name;

        while(std::getline(m_file,m_line))
        {
            counter++;
            if(m_line.find(match)!=std::string::npos)
            {
                if(!std::getline(m_file,m_line))
                {
                    return 0;
                }

                counter++;

                size_t first  = m_line.find(SCOTS_FH_SEPERATOR);
                size_t second = m_line.find(SCOTS_FH_SEPERATOR,first+1);

                if((first == std::string::npos)||(second == std::string::npos))
                {
                    return 0;
                }

                std::istringstream stream_1(m_line.substr(0,first));
                std::istringstream stream_2(m_line.substr(first+1,(second-first)));

                abs_type cells;
                abs_type indices;

                stream_1 >> cells;
                stream_2 >> indices;

                map.resize(cells,indices);
                map.initialize();

                break;
            }
        }
        if(!std::getline(m_file,m_line))
        {
            map.resize(0,0);
            return 0;
        }
        counter++;
        if(m_line.find(SCOTS_FH_BEGIN)== std::string::npos)
        {
            map.resize(0,0);
            return 0;
        }
        while(std::getline(m_file,m_line))
        {
            counter++;
            abs_type cell;
            abs_type index;

            if(m_line.find(SCOTS_FH_END) != std::string::npos)
            {
                map.checkNrOfValidPairs();
                return  counter;
            }
            else
            {
                size_t placeofcomma = m_line.find(SCOTS_FH_SEPERATOR);
                if(placeofcomma == std::string::npos)
                {
                    map.resize(0,0);
                    return 0;
                }
                std::istringstream stream_1(m_line.substr(0,placeofcomma));
                std::istringstream stream_2(m_line.substr(placeofcomma+1));

                stream_1 >> cell;
                stream_2 >> index;

                map[cell] = index;
            }
        }

        map.resize(0,0);
        return 0;
    }
    else
    {
        map.resize(0,0);
        return 0;
    }
}

size_t FileReader::get_nonDeterministicMap(nonDeterministicMap& map, std::string name, size_t offset)
{
    backToFirstLine();
    if(skipOffset(offset))
    {
        size_t counter = 0;
        std::string match = SCOTS_FH_NDETMAP + name;

        while(std::getline(m_file,m_line))
        {
            counter++;
            if(m_line.find(match)!=std::string::npos)
            {
                if(!std::getline(m_file,m_line))
                {
                    return 0;
                }

                counter++;

                size_t first  = m_line.find(SCOTS_FH_SEPERATOR);
                size_t second = m_line.find(SCOTS_FH_SEPERATOR,first+1);

                if((first == std::string::npos)||(second == std::string::npos))
                {
                    return 0;
                }

                std::istringstream stream_1(m_line.substr(0,first));
                std::istringstream stream_2(m_line.substr(first+1,(second-first)));

                abs_type cells;
                abs_type indices;

                stream_1 >> cells;
                stream_2 >> indices;

                map.resize(cells,indices);

                break;
            }
        }
        if(!std::getline(m_file,m_line))
        {
            map.resize(0,0);
            return 0;
        }
        counter++;
        if(m_line.find(SCOTS_FH_BEGIN)== std::string::npos)
        {
            map.resize(0,0);
            return 0;
        }

        while(std::getline(m_file,m_line))
        {
            counter++;
            abs_type cell;
            abs_type index;

            if(m_line.find(SCOTS_FH_END) != std::string::npos)
            {
                map.checkNrOfValidPairs();
                return  counter;
            }
            else
            {
                size_t last = 0;
                size_t next = 0;

                next = m_line.find(";",last);
                std::istringstream stream_cell(m_line.substr(last,next-last));
                stream_cell >> cell;

                last = next;
                next = m_line.find(";",last+1);

                while(next != std::string::npos)
                {
                    std::istringstream stream_index(m_line.substr(last+1,next-last));
                    stream_index >> index;

                    map.set(cell,index,true);

                    last = next;
                    next = m_line.find(";",last+1);
                }
                std::istringstream stream_last_index(m_line.substr(last+1));
                stream_last_index >> index;

                map.set(cell,index,true);
            }

        }

        map.resize(0,0);
        return 0;
    }
    else
    {
        map.resize(0,0);
        return 0;
    }
}

bool FileReader::skipOffset(size_t& offset)
{
    for(size_t index = 0; index<offset; index++)
    {
        if(!(m_file.ignore(std::numeric_limits<std::streamsize>::max(),'\n')))
        {
            return false;
        }
    }

    return true;
}

void FileReader::backToFirstLine()
{
    m_file.clear();
    m_file.seekg(0,std::ios::beg);
}

template<class T>
size_t FileReader::get_MEMBER(std::string memberName,T& member, size_t offset)
{
    backToFirstLine();
    if(skipOffset(offset))
    {
        size_t counter =0;
        std::string match = SCOTS_FH_MEMBER + memberName;

        while(std::getline(m_file,m_line))
        {
            counter++;
            if(m_line.find(match)!=std::string::npos)
            {
                if(std::getline(m_file,m_line))
                {
                    std::istringstream stream(m_line);
                    stream >> member;
                    counter++;
                    return counter;
                }
                else
                {
                    return 0;
                }
            }
        }

        return 0;
    }
    else
    {
        return 0;
    }
}
template<class T>
size_t FileReader::get_VECTOR(std::string vectorName, std::vector<T>& vector, size_t offset)
{
    backToFirstLine();
    if(skipOffset(offset))
    {
        size_t counter  = 0;
        size_t size     = 0;
        vector.clear();
        std::string match = SCOTS_FH_VECTOR + vectorName;

        while(std::getline(m_file,m_line))
        {
            counter++;
            if(m_line.find(match)!=std::string::npos)
            {
                if(std::getline(m_file,m_line))
                {
                    if(m_line.find(SCOTS_FH_BEGIN)!=std::string::npos)
                    {
                        std::istringstream stream(m_line.substr(m_line.find(":")+1));
                        stream >> size;
                        counter++;

                        break;
                    }
                    else
                    {
                        return 0;
                    }
                }
                else
                {
                    return 0;
                }
            }
        }
        if(m_file.eof())
        {
            return 0;
        }
        else
        {
            vector.resize(size);
            for(size_t index = 0; index < size; index++)
            {
                if(std::getline(m_file,m_line))
                {
                    if(m_line.find(SCOTS_FH_SYMBOL)!=std::string::npos)
                    {
                        vector.clear();
                        return 0;
                    }
                    else
                    {
                        counter++;
                        std::istringstream stream(m_line);
                        stream >> vector[index];
                    }
                }
                else
                {
                    vector.clear();
                    return 0;
                }
            }
            if(std::getline(m_file,m_line))
            {
                if(m_line.find(SCOTS_FH_END)!=std::string::npos)
                {
                    counter++;
                    return counter;
                }
                else
                {
                    return 0;
                }
            }
            else
            {
                return 0;
            }
        }
    }
    else
    {
        return 0;
    }
}
template<class T>
size_t FileReader::get_ARRAY(std::string arrayName,T& array, size_t arraySize, size_t offset)
{
    backToFirstLine();
    if(skipOffset(offset))
    {
        size_t counter  = 0;
        size_t size     = 0;
        std::string match = SCOTS_FH_ARRAY + arrayName;

        while(std::getline(m_file,m_line))
        {
            counter++;
            if(m_line.find(match)!=std::string::npos)
            {
                if(std::getline(m_file,m_line))
                {
                    if(m_line.find(SCOTS_FH_BEGIN)!=std::string::npos)
                    {
                        std::istringstream stream(m_line.substr(m_line.find(":")+1));
                        stream >> size;
                        if(size == arraySize)
                        {
                            counter++;
                            break;
                        }
                        else
                        {
                            return 0;
                        }

                    }
                    else
                    {
                        return 0;
                    }
                }
                else
                {
                    return 0;
                }
            }
        }
        if(m_file.eof())
        {
            return 0;
        }
        else
        {
            for(size_t index = 0; index < size; index++)
            {
                if(std::getline(m_file,m_line))
                {
                    if(m_line.find(SCOTS_FH_SYMBOL)!=std::string::npos)
                    {
                        return 0;
                    }
                    else
                    {
                        counter++;
                        std::istringstream stream(m_line);
                        stream >> array[index];
                    }
                }
                else
                {
                    return 0;
                }
            }
            if(std::getline(m_file,m_line))
            {
                if(m_line.find(SCOTS_FH_END)!=std::string::npos)
                {
                    counter++;
                    return counter;
                }
                else
                {
                    return 0;
                }
            }
            else
            {
                return 0;
            }
        }
    }
    else
    {
        return 0;
    }
}

} /*end of namepace scots*/

#endif /* FILEHANDLER_HH_ */
