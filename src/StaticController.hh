/*
 * StaticController.hh
 *
 *  created: Dec 2016
 *   author: Frederik Kunik
 *           Matthias Rungger
 *           
 */

/** @file **/

#ifndef STATICCONTROLLER_HH_
#define STATICCONTROLLER_HH_

#include <vector>
#include <array>
#include <exception>
#include <chrono>

#include "UniformGrid.hh"
#include "WinningDomain.hh"

#define  SCOTS_SC_DEFAULT_NAME  "StaticController"
#define  SCOTS_SC_TYPE          "STATICCONTROLLER"

/** @namespace scots **/ 
namespace scots{

//bool writeStaticControllerToFile(const UniformGrid&, const UniformGrid&, const StaticController& ,std::string = "");

/**
 * @brief StaticController class to simualte the closed loop
 **/
class StaticController {
public:
  ///** @cond  EXCLUDE from doxygen **/
  ///* default constructor */
  //StaticController()=default;                      
  ///* destructor */
  //~StaticController()=default;
  ///* copy constructor deleted (cannot be copied) */
  //StaticController(const StaticController&);
  ///* move constructor */
  //StaticController(StaticController&&);
  ///* copy assignment operator */
  //StaticController& operator=(const StaticController&)=delete; 
  ///* move assignment operator */
  //StaticController& operator=(StaticController&&)=default;
  ///* @endcond */

  //StaticController& operator=(const StaticController&); 
  //StaticController& operator=(StaticController&&); 

  /* @brief controller constructor */

  StaticController(const UniformGrid& state_grid,
									 const UniformGrid& input_grid,
                   WinningDomain&& wining_domain,
                   const std::string filename = "") {
    m_state_grid = state_grid;
    m_input_grid = input_grid;
		m_winning_domain = std::move(wining_domain);
    m_file_name  = filename;
  }

  template<class state_type, class input_type>
  std::vector<input_type> get_control(const state_type &x) {
    /* abstract state index */
    abs_type i;
    m_state_grid.xtoi(i,x);
    std::vector<abs_type> abs_inputs = m_winning_domain.get_inputs(i);

    if(!m_winning_domain.is_winning(i)) {
      std::ostringstream os;
      os << "\nscots::StaticController: state ";
      for(int i=0; i<m_state_grid.getDimension(); i++) {
        os << x[i] << " ";
      }
      os << "is out of winning domain: no progress possible.";

      throw std::runtime_error(os.str().c_str());
    }
    std::vector<input_type> inputs(abs_inputs.size());

    input_type u;
    for(abs_type i=0; i<abs_inputs.size();i++) {
      m_input_grid.itox(abs_inputs[i],u);
      inputs[i]=u;
    }
    return inputs;
  }

  //bool writeToFile(std::string = "");
  //bool readFromFile(FileReader&);

  static std::string createDefaultFilename() {
		std::string filename = SCOTS_SC_DEFAULT_NAME;
		/* create time stamp */
    auto time_point = std::chrono::system_clock::now(); 
    std::time_t now_c = std::chrono::system_clock::to_time_t(time_point);
    std::ostringstream os;
    os << "_" << now_c;
		filename.append(os.str());

		return filename;
	}

protected:
  UniformGrid m_input_grid;
  UniformGrid m_state_grid;
  std::string m_file_name;

private:
  WinningDomain m_winning_domain;

};

//bool writeStaticControllerToFile(const UniformGrid& state_alphabet, const UniformGrid& input_alphabet, const StaticController& controller, std::string filename) {
//  if(filename == "") {
//    filename = nonDeterministic_staticController::createDefaultFilename();
//  }
//  scots::FileWriter writer(filename);
//
//  if(!writer.create()) {
//    return false;
//  }
//  if(!writer.add_TYPE(SCOTS_SC_TYPE)) {
//    return false;
//  }
//  if(!input_alphabet.addGridToFile(writer)) {
//    return false;
//  }
//  if(!state_alphabet->addGridToFile(writer)) {
//    return false;
//  }
//  if(!writer.open()) {
//    return false;
//  }
//  if(!writer.add_nonDeterministicMap(*map)) {
//    return false;
//  }
//  writer.close();
//
//  return true;
//}



//bool readFromFile(FileReader& reader)
//{
//  size_t offsetOfController = 0;
//  size_t offsetOfInputGrid  = 0;
//  size_t offsetOfStateGrid  = 0;
//  std::string type;
//  double version = 0;
//
//  m_validController = false;
//
//  if(!reader.open())
//  {
//    return false;
//  }
//  reader.get_VERSION(version,0);
//  if(version != SCOTS_FH_VERSION)
//  {
//    //WARNING
//  }
//  offsetOfController = reader.get_TYPE(type,0);
//  if(!offsetOfController || (type != SCOTS_DSC_TYPE))
//  {
//    return false;
//  }
//  offsetOfInputGrid += reader.get_TYPE(type,offsetOfController);
//  if(!offsetOfInputGrid || (type != SCOTS_UG_TYPE))
//  {
//    return false;
//  }
//  reader.close();
//  offsetOfInputGrid += offsetOfController;
//  if(!m_inputGrid.readFromGridFile(reader,offsetOfInputGrid))
//  {
//    return false;
//  }
//  if(!reader.open())
//  {
//    return false;
//  }
//  offsetOfStateGrid = reader.get_TYPE(type,offsetOfInputGrid);
//  if(!offsetOfStateGrid || (type != SCOTS_UG_TYPE))
//  {
//    return false;
//  }
//  reader.close();
//  offsetOfStateGrid += offsetOfInputGrid;
//  if(!m_stateGrid.readFromGridFile(reader,offsetOfStateGrid))
//  {
//    return false;
//  }
//  if(!reader.open())
//  {
//    return false;
//  }
//  if(!reader.get_nonDeterministicMap(m_map,"",offsetOfStateGrid))
//  {
//    return false;
//  }
//  reader.close();
//
//  m_validController   = validation();
//
//  return true;
//}

//bool nonDeterministic_staticController::validation()
//{
//  uint8_t flag = 0;
//  if((m_inputGrid.getTotalNrOfGridPoints() !=0 )& (m_stateGrid.getTotalNrOfGridPoints() != 0) & (m_map.NrOfCells() !=0) & (m_map.NrOfInputIndices() != 0))
//  {
//    flag++;
//  }
//  if(m_inputGrid.getTotalNrOfGridPoints() == m_map.NrOfInputIndices())
//  {
//    flag++;
//  }
//  if(m_stateGrid.getTotalNrOfGridPoints() == m_map.NrOfCells())
//  {
//    flag++;
//  }
//  if(flag == 3)
//  {
//    m_validController = true;
//  }
//  else
//  {
//    m_validController = false;
//  }
//
//  return m_validController;
//}
//bool deterministic_staticController::writeToFile(std::string filename)
//{
//  return writeStaticControllerToFile(&m_inputGrid,&m_stateGrid,&m_map,filename);
//}

//bool deterministic_staticController::readFromFile(FileReader& reader)
//{
//  size_t offsetOfController = 0;
//  size_t offsetOfInputGrid  = 0;
//  size_t offsetOfStateGrid  = 0;
//  std::string type;
//  double version = 0;
//
//  m_validController = false;
//
//  if(!reader.open())
//  {
//    return false;
//  }
//  reader.get_VERSION(version,0);
//  if(version != SCOTS_FH_VERSION)
//  {
//    //WARNING
//  }
//  offsetOfController = reader.get_TYPE(type,0);
//  if(!offsetOfController || (type != SCOTS_NDSC_TYPE))
//  {
//    return false;
//  }
//  offsetOfInputGrid += reader.get_TYPE(type,offsetOfController);
//  if(!offsetOfInputGrid || (type != SCOTS_UG_TYPE))
//  {
//    return false;
//  }
//  reader.close();
//  offsetOfInputGrid += offsetOfController;
//  if(!m_inputGrid.readFromGridFile(reader,offsetOfInputGrid))
//  {
//    return false;
//  }
//  if(!reader.open())
//  {
//    return false;
//  }
//  offsetOfStateGrid = reader.get_TYPE(type,offsetOfInputGrid);
//  if(!offsetOfStateGrid || (type != SCOTS_UG_TYPE))
//  {
//    return false;
//  }
//  reader.close();
//  offsetOfStateGrid += offsetOfInputGrid;
//  if(!m_stateGrid.readFromGridFile(reader,offsetOfStateGrid))
//  {
//    return false;
//  }
//  if(!reader.open())
//  {
//    return false;
//  }
//  if(!reader.get_deterministicMap(m_map,"",offsetOfStateGrid))
//  {
//    return false;
//  }
//  reader.close();
//
//  m_validController   = validation();
//
//  return true;
//}
//

} /* end of namespace scots */

#endif /* STATICCONTROLLER_HH_ */
