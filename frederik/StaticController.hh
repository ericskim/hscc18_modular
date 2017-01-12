#ifndef STATICCONTROLLER_HH_
#define STATICCONTROLLER_HH_

#include <vector>
#include <array>
#include <stdint.h>
#include <ctime>

#include "UniformGrid_FK.hh"
#include "Maps.hh"

#define  SCOTS_SC_DEFAULT_NAME  "StaticController"
#define  SCOTS_DSC_TYPE         "STATICCONTROLLER_DETERMINISTIC"
#define  SCOTS_NDSC_TYPE        "STATICCONTROLLER_NONDETERMINISTIC"

namespace scots{


bool writeStaticControllerToFile(UniformGrid* const, UniformGrid* const, deterministicMap* const,std::string = "");
bool writeStaticControllerToFile(UniformGrid* const, UniformGrid* const, nonDeterministicMap* const,std::string = "");

/*!
 * \brief The nonDeterministic_staticController class
 */
class nonDeterministic_staticController
{
public:
    nonDeterministic_staticController();
    nonDeterministic_staticController(const UniformGrid&, const UniformGrid&, const nonDeterministicMap&, const std::string = "");
    nonDeterministic_staticController(const nonDeterministic_staticController&);
    ~nonDeterministic_staticController();

    nonDeterministic_staticController& operator=(const nonDeterministic_staticController&); //!< copy asignment operator

    virtual bool writeToFile(std::string = "");
    virtual bool readFromFile(FileReader&);

    bool validation();
    bool IsValid();


    template< class grid_point_t_in, class grid_point_t_state>
    bool control(const grid_point_t_state,std::vector<grid_point_t_in>&);

    static std::string createDefaultFilename();

protected:
    UniformGrid m_inputGrid;
    UniformGrid m_stateGrid;
    std::string m_fileName;
    bool        m_validController;

private:
    nonDeterministicMap m_map;

};

class deterministic_staticController :public nonDeterministic_staticController
{
public:
    deterministic_staticController();
    deterministic_staticController(const UniformGrid&, const UniformGrid&, const deterministicMap&, const std::string = "");
    deterministic_staticController(const deterministic_staticController&);
    ~deterministic_staticController();

    deterministic_staticController& operator=(const deterministic_staticController&); //!< copy asignment operator

    virtual bool writeToFile(std::string = "");
    virtual bool readFromFile(FileReader&);

    template< class grid_point_t_in, class grid_point_t_state>
    bool control(const grid_point_t_state,grid_point_t_in&);

private:
    deterministicMap m_map;
};

bool writeStaticControllerToFile(UniformGrid* const inputGrid, UniformGrid* const stateGrid, nonDeterministicMap* const map, std::string filename)
{
    if(filename == "")
    {
       filename = nonDeterministic_staticController::createDefaultFilename();
    }
    scots::FileWriter writer(filename);

    if(!writer.create())
    {
        return false;
    }
    if(!writer.add_TYPE(SCOTS_DSC_TYPE))
    {
        return false;
    }
    if(!inputGrid->addGridToFile(writer))
    {
        return false;
    }
    if(!stateGrid->addGridToFile(writer))
    {
        return false;
    }
    if(!writer.open())
    {
        return false;
    }
    if(!writer.add_nonDeterministicMap(*map))
    {
        return false;
    }

    writer.close();

    return true;
}

bool writeStaticControllerToFile(UniformGrid* const inputGrid, UniformGrid* const stateGrid, deterministicMap* const map, std::string filename)
{
    if(filename == "")
    {
       filename = deterministic_staticController::createDefaultFilename();
    }
    scots::FileWriter writer(filename);

    if(!writer.create())
    {
        return false;
    }
    if(!writer.add_TYPE(SCOTS_NDSC_TYPE))
    {
        return false;
    }
    if(!inputGrid->addGridToFile(writer))
    {
        return false;
    }
    if(!stateGrid->addGridToFile(writer))
    {
        return false;
    }
    if(!writer.open())
    {
        return false;
    }
    if(!writer.add_deterministicMap(*map))
    {
        return false;
    }

    writer.close();

    return true;
}

nonDeterministic_staticController::nonDeterministic_staticController()
{
    m_validController   = false;
}

nonDeterministic_staticController::nonDeterministic_staticController(const UniformGrid& inputGrid, const UniformGrid& stateGrid, const nonDeterministicMap& domain, const std::string filename) :nonDeterministic_staticController()
{
    m_inputGrid = inputGrid;
    m_stateGrid = stateGrid;
    m_map       = domain;
    m_fileName  = filename;

    m_validController = validation();
}

nonDeterministic_staticController &nonDeterministic_staticController::operator=(const nonDeterministic_staticController &other)
{
    m_inputGrid = other.m_inputGrid;
    m_stateGrid = other.m_stateGrid;
    m_map       = other.m_map;
    m_fileName  = other.m_fileName;

    m_validController   = validation();

    return *this;
}

nonDeterministic_staticController::nonDeterministic_staticController(const nonDeterministic_staticController &model) :nonDeterministic_staticController()
{
    *this = model;
}

scots::nonDeterministic_staticController::~nonDeterministic_staticController()
{
}

bool nonDeterministic_staticController::writeToFile(const std::string filename)
{
    return writeStaticControllerToFile(&m_inputGrid,&m_stateGrid,&m_map,filename);
}

bool nonDeterministic_staticController::readFromFile(FileReader& reader)
{
    size_t offsetOfController = 0;
    size_t offsetOfInputGrid  = 0;
    size_t offsetOfStateGrid  = 0;
    std::string type;
    double version = 0;

    m_validController = false;

    if(!reader.open())
    {
        return false;
    }
    reader.get_VERSION(version,0);
    if(version != SCOTS_FH_VERSION)
    {
        //WARNING
    }
    offsetOfController = reader.get_TYPE(type,0);
    if(!offsetOfController || (type != SCOTS_DSC_TYPE))
    {
        return false;
    }
    offsetOfInputGrid += reader.get_TYPE(type,offsetOfController);
    if(!offsetOfInputGrid || (type != SCOTS_UG_TYPE))
    {
        return false;
    }
    reader.close();
    offsetOfInputGrid += offsetOfController;
    if(!m_inputGrid.readFromGridFile(reader,offsetOfInputGrid))
    {
        return false;
    }
    if(!reader.open())
    {
        return false;
    }
    offsetOfStateGrid = reader.get_TYPE(type,offsetOfInputGrid);
    if(!offsetOfStateGrid || (type != SCOTS_UG_TYPE))
    {
        return false;
    }
    reader.close();
    offsetOfStateGrid += offsetOfInputGrid;
    if(!m_stateGrid.readFromGridFile(reader,offsetOfStateGrid))
    {
        return false;
    }
    if(!reader.open())
    {
        return false;
    }
    if(!reader.get_nonDeterministicMap(m_map,"",offsetOfStateGrid))
    {
        return false;
    }
    reader.close();

    m_validController   = validation();

    return true;
}

bool nonDeterministic_staticController::validation()
{
    uint8_t flag = 0;
    if((m_inputGrid.getTotalNrOfGridPoints() !=0 )& (m_stateGrid.getTotalNrOfGridPoints() != 0) & (m_map.NrOfCells() !=0) & (m_map.NrOfInputIndices() != 0))
    {
        flag++;
    }
    if(m_inputGrid.getTotalNrOfGridPoints() == m_map.NrOfInputIndices())
    {
        flag++;
    }
    if(m_stateGrid.getTotalNrOfGridPoints() == m_map.NrOfCells())
    {
        flag++;
    }
    if(flag == 3)
    {
        m_validController = true;
    }
    else
    {
        m_validController = false;
    }

    return m_validController;
}

bool nonDeterministic_staticController::IsValid()
{
    return m_validController;
}

std::string scots::nonDeterministic_staticController::createDefaultFilename()
{
    std::string filename = SCOTS_SC_DEFAULT_NAME;

    time_t now  = time(0);
    char* dt    = ctime(&now);

    filename.append(dt);

    return filename;
}

template< class grid_point_t_in, class grid_point_t_state>
bool nonDeterministic_staticController::control(const grid_point_t_state feed,std::vector<grid_point_t_in>& input)
{
    abs_type cellIndex;
    std::vector<abs_type> inputIndices;

    input.clear();

    m_stateGrid.xtoi<grid_point_t_state>(cellIndex,feed);

    inputIndices = m_map.getIndices(cellIndex);

    if(!inputIndices.size())
    {
        return false;
    }
    else
    {
        input.resize(inputIndices.size());

        for(size_t index = 0; index < inputIndices.size();index++)
        {
            m_inputGrid.itox<grid_point_t_in>(inputIndices[index],input[index]);
        }

        return true;
    }
}

deterministic_staticController::deterministic_staticController()
{
    m_validController = false;
}

deterministic_staticController::deterministic_staticController(const UniformGrid& inputGrid, const UniformGrid& stateGrid, const deterministicMap& domain, const std::string filename) :deterministic_staticController()
{
    m_inputGrid = inputGrid;
    m_stateGrid = stateGrid;
    m_map       = domain;
    m_fileName  = filename;

    m_validController = validation();
}

deterministic_staticController::deterministic_staticController(const deterministic_staticController& model) :deterministic_staticController()
{
    *this = model;
}

deterministic_staticController::~deterministic_staticController()
{
}

deterministic_staticController &deterministic_staticController::operator=(const deterministic_staticController& other)
{
    m_inputGrid = other.m_inputGrid;
    m_stateGrid = other.m_stateGrid;
    m_map       = other.m_map;
    m_fileName  = other.m_fileName;

    m_validController   = validation();

    return *this;
}

bool deterministic_staticController::writeToFile(std::string filename)
{
    return writeStaticControllerToFile(&m_inputGrid,&m_stateGrid,&m_map,filename);
}

bool deterministic_staticController::readFromFile(FileReader& reader)
{
    size_t offsetOfController = 0;
    size_t offsetOfInputGrid  = 0;
    size_t offsetOfStateGrid  = 0;
    std::string type;
    double version = 0;

    m_validController = false;

    if(!reader.open())
    {
        return false;
    }
    reader.get_VERSION(version,0);
    if(version != SCOTS_FH_VERSION)
    {
        //WARNING
    }
    offsetOfController = reader.get_TYPE(type,0);
    if(!offsetOfController || (type != SCOTS_NDSC_TYPE))
    {
        return false;
    }
    offsetOfInputGrid += reader.get_TYPE(type,offsetOfController);
    if(!offsetOfInputGrid || (type != SCOTS_UG_TYPE))
    {
        return false;
    }
    reader.close();
    offsetOfInputGrid += offsetOfController;
    if(!m_inputGrid.readFromGridFile(reader,offsetOfInputGrid))
    {
        return false;
    }
    if(!reader.open())
    {
        return false;
    }
    offsetOfStateGrid = reader.get_TYPE(type,offsetOfInputGrid);
    if(!offsetOfStateGrid || (type != SCOTS_UG_TYPE))
    {
        return false;
    }
    reader.close();
    offsetOfStateGrid += offsetOfInputGrid;
    if(!m_stateGrid.readFromGridFile(reader,offsetOfStateGrid))
    {
        return false;
    }
    if(!reader.open())
    {
        return false;
    }
    if(!reader.get_deterministicMap(m_map,"",offsetOfStateGrid))
    {
        return false;
    }
    reader.close();

    m_validController   = validation();

    return true;
}

template< class grid_point_t_in, class grid_point_t_state>
bool deterministic_staticController::control(const grid_point_t_state feed, grid_point_t_in& input)
{
    abs_type cellIndex;
    abs_type inputIndex;

    m_stateGrid.xtoi<grid_point_t_state>(cellIndex,feed);

    inputIndex = m_map[cellIndex];

    if(inputIndex == Scots_Domain_InvalidIndex)
    {
        return false;
    }
    else
    {
        m_inputGrid.itox<grid_point_t_in>(inputIndex,input);
        return true;
    }
}

} /* end of namespace scots */

#endif /* STATICCONTROLLER_HH_ */
