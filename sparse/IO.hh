#ifndef IO_HH_
#define IO_HH_

#include <iostream>
#include <limits>
#include "UniformGrid.hh"
#include "StaticController.hh"
#include "TransitionSystem.hh"
#include "AbstractionGB.hh"
#include "ReachabilityGame.hh"
namespace scots {
/*
 * class: IO
 *
 *
 */

// short changed to unsigned short and for pointer state_size_t changed to size_t
/* class IO */
class IO {
public:

/* constructor */
IO(){
}

/* destructor */
~IO(){
}

/* function: createUniformGridParameterString */
template<class U>
static std::string createUniformGridParameterString(const UniformGrid<U>* gs, const char* key="scots") {
    /* produce string with parameters to be written into file */
    std::stringstream dim;
    std::stringstream z;
    std::stringstream eta;
    std::stringstream first;
    std::stringstream nofGridPoints;

    dim <<"#"<<key<<" dimension: " << gs->dim_ << std::endl;
    z <<"#"<<key<<" measurement error bound: ";
    eta <<"#"<<key<<" grid parameter eta: ";
    first <<"#"<<key<<" coordinate of first grid point: ";
    nofGridPoints <<"#"<<key<<" number of grid points (per dim): ";

    for(size_t i=0; i<gs->dim_; i++) {
        z << gs->z_[i] << " ";
        eta << gs->eta_[i] << " ";
        first << gs->firstGridPoint_[i] << " ";
        nofGridPoints << gs->nofGridPoints_[i] << " ";
    }
    z << std::endl;
    eta << std::endl;
    first << std::endl;
    nofGridPoints << std::endl;

    std::stringstream paramss;
    paramss << dim.str() \
        << eta.str() \
        << z.str() \
        << first.str() \
        << nofGridPoints.str();
    return paramss.str();
}

/* function:  readDimensionFromFile
 * a function to read the dimension of the uniform grid from file */
template<class T>
//static void readDimensionFromFile(UniformGrid<T> &ss, const char& filename, const char* key="scots")
static void readDimensionFromFile(UniformGrid<T> *ss, const char *filename, const char* key="scots") {
    ss->dim_=0;
    /* open file */
    std::ifstream file(filename);
    if(!file.good()) {
        std::ostringstream os;
        os << "IO: UniformGrid: Error: Unable to open file:" << filename << "'.";
        throw std::runtime_error(os.str().c_str());
    }
    /* create key string: parameters are only read from lines that start with
     * '#key' */
    std::string mykey("#");
    mykey.append(key);

    /* read dimension from file */
    std::string line;
    while(!file.eof()) {
        std::getline(file,line);
        if (line.substr(0,mykey.length())==mykey) {
            if(line.find("dimension")!=std::string::npos) {
                std::istringstream sline(line.substr(line.find(":")+1));
                sline >> ss->dim_;
            }
        }
    }
    if(ss->dim_==0) {
        std::ostringstream os;
        os << "IO: UniformGrid: Error: Could not read dimension from file: " << filename << ". ";
        os << "Was " << filename << " created with scots::UniformGrid::writeToFile?";
        throw std::runtime_error(os.str().c_str());
    }
}

/* function:  readMembersFromFile
 * a function to read the uniform grid data from file */
template<class T>
// static void readMembersFromFile(UniformGrid<T> &ss, const char &filename, const char* key="scots")
static void readMembersFromFile(UniformGrid<T> *ss, const char *filename, const char* key="scots") {

    std::ifstream file(filename);
    if(!file.good()) {
        std::ostringstream os;
        os << "IO: UniformGrid: Error: Unable to open file:" << filename << "'.";
        throw std::runtime_error(os.str().c_str());
    }
    /* create key string: parameters are only read from lines that start with
     * '#key' */
    std::string mykey("#");
    mykey.append(key);

    /* read eta/first/last/no of grid points */
    std::string line;
    int check=0;
    while(!file.eof()) {
        std::getline(file,line);
        if (line.substr(0,mykey.length())==mykey) {
            /* read eta */
            if(line.find("eta")!=std::string::npos) {
                check++;
                std::istringstream sline(line.substr(line.find(":")+1));
                for(size_t i=0; i<ss->dim_; i++)
                    sline >> ss->eta_[i];
            }
            /* read z */
            if(line.find("measurement")!=std::string::npos) {
                check++;
                std::istringstream sline(line.substr(line.find(":")+1));
                for(size_t i=0; i<ss->dim_; i++)
                    sline >> ss->z_[i];
            }
            /* read first grid point*/
            if(line.find("first")!=std::string::npos) {
                check++;
                std::istringstream sline(line.substr(line.find(":")+1));
                for(size_t i=0; i<ss->dim_; i++)
                    sline >> ss->firstGridPoint_[i];
            }
            /* read no of grid points */
            if(line.find("number")!=std::string::npos) {
                check++;
                std::istringstream sline(line.substr(line.find(":")+1));
                for(size_t i=0; i<ss->dim_; i++) {
                    sline >> ss->nofGridPoints_[i];
                }
            }
        }
        if(check==4)
            break;
    }
    if(check<4) {
        std::ostringstream os;
        os << "IO: UniformGrid: Error: Could not read all parameters from file: " << filename << ". ";
        os << "Was " << filename << " created with scots::UniformGrid::writeToFile?";
        throw std::runtime_error(os.str().c_str());
    }
    /* read the abstract set indices */
    file.clear();
    file.seekg(0, std::ios::beg);
    int now=0;
    size_t idx;
    while(!file.eof()) {
        std::getline(file,line);
        if(now==1) {
            try {
//    idx = std::stoull(line);
                idx =(size_t)strtol(line.c_str(),NULL,10);
                ss->abstractSet_.insert(idx);
            }
            catch(...) {
            }
        }
        if(line.find("#AbstractSet")!=std::string::npos)
            now=1;
    }
}
/* function: writeToFile()
 * gs: UniformGrid
 * filename: name of file */
template<class U>
static void writeToFile(const UniformGrid<U>* gs, const char* filename, const char* key="scots") {
    /* write string to file */
    std::ofstream file(filename);
    if(!file.good()) {
        std::ostringstream os;
        os << "IO: UniformGrid: Error: Unable to open file:" << filename << "'.";
        throw std::runtime_error(os.str().c_str());
    }
    file << "################################################################################" << std::endl;
    file << "########### SCOTS: uniform grid parameter data      ############" << std::endl;
    file << "################################################################################" << std::endl;
    /* append the parameter data of the uniform grid */
    file << createUniformGridParameterString(gs,key) << std::endl;
    file << "#AbstractSet" << std::endl;

    std::set<abs_type>::iterator it;
    for (it=gs->abstractSet_.begin(); it != gs->abstractSet_.end(); ++it)
        file << *it << std::endl;
    file.close();
    std::cout << "Uniform grid saved to file: "<< filename << std::endl;
}
/* function: writeToFile()
 * ts: TransitionSystem
 * filename: name of file */
static void writeToFile(const TransitionSystem* ts, const char* filename) {
    /* open file */
    std::ofstream file(filename);
    if(!file.good()) {
        std::ostringstream os;
        os << "IO: TransitionSystem: Error: Unable to open file:" << filename << "'.";
        throw std::runtime_error(os.str().c_str());
    }
    /* write to file */
    file << "################################################################################" << std::endl;
    file << "###########  SCOTS: TransitionSystem information    ##########"<< std::endl;
    file << "################################################################################" << std::endl;
    file << std::endl;
    file << "#Transition System" << std::endl;
    file << "#no of states in the transition system: " <<ts->N_ << std::endl;
    file << "#no of labels in the transition system: " << ts->M_ << std::endl;
    file << "#no of transitions in the transition system: " << ts->T_ << std::endl;

    if(ts->pre_!=NULL) {
        file << "#transition relation:" << std::endl;
        file << "#\"idx of source state\"  \"idx of input\"  \"idx of target state\" " << std::endl;
        for(size_t k=0; k<ts->N_; k++) {
            for(size_t j=0; j<ts->M_; j++) {
                if(ts->noPre_[k] && ts->noPre_[k*ts->M_+j]!=0) {
                    for(size_t no=0; no<(size_t)ts->noPre_[k*ts->M_+j]; no++) {
                        size_t pos=ts->prePointer_[k*ts->M_+j]+no;
                        size_t i=ts->pre_[pos];
                        file << i << " " << j << " " << k << std::endl;
                    }
                }
            }
        }
    }
    else if(ts->post_!=NULL) {
        file << "#transition relation:" << std::endl;
        file << "#\"idx of source state\"  \"idx of input\"  \"idx of target state\" " << std::endl;
        for(size_t i=0; i<ts->N_; i++) {
            for(size_t j=0; j<ts->M_; j++) {
                if(ts->noPost_[i] && ts->noPost_[i*ts->M_+j]!=0) {
                    for(size_t no=0; no<(size_t)ts->noPost_[i*ts->M_+j]; no++) {
                        size_t pos=ts->postPointer_[i*ts->M_+j]+no;
                        size_t k=ts->post_[pos];
                        file << i << " " << j << " " << k << std::endl;
                    }
                }
            }
        }
    }
    else {
        std::ostringstream os;
        os << "IO: TransitionSystem: Error: Unable to write to file. Transition relation is empty, i.e.,  pre_ and post_ are NULL.";
        throw std::runtime_error(os.str().c_str());
    }
    file.close();
    std::cout << "Transition System saved to file: "<< filename << std::endl;
}

/* function: writeControllerToFile()
* only save the state-input pairs */
//const ReachabilityGame* reach,
static void writeControllerToFile(const ReachabilityGame* reach, const char* filename) {
    /* open file */
    std::ofstream file(filename);
    if(!file.good()) {
        std::ostringstream os;
        os << "IO: Game: Error: Unable to open file:" << filename << "'.";
        throw std::runtime_error(os.str().c_str());
    }
    /* write to file */
    file << "################################################################################" << std::endl;
    file << "#############    SCOTS: controller information    ############" << std::endl;
    file << "################################################################################" << std::endl;
    file << std::endl;
    file << "#Controller" << std::endl;
    file << "#no of states in the transition system: " << reach->N_ << std::endl;
    file << "#no of labels in the transition system: " << reach->M_ << std::endl;
    file << "#\"idx of state\"  \"value of value function\"  \"idx of input\" " << std::endl;

    bool m=false;
    for(size_t i=0; i<reach->N_; i++) {
        for(size_t j=0; j<reach->M_; j++) {
            if(reach->inputs_[i]!=-1) {
                if(!m) {
                    file << i << " " << reach->val_[i] <<" ";
                    m=true;
                }
                file << j << " ";
            }
            m=false;
            file << std::endl;
        }
    }
    file.close();
    std::cout << "controller saved to file: "<< filename << std::endl;
}


/* function: writeControllerToFile()  (save the graph and grid points)
 * ss,is: UniformGrid
 * con: StaticController
 * filename: name of file */
// template<class T=std::array<double,1>, class U=std::array<double,1>>
//const ReachabilityGame* reach,
template <class T, class U>
static void writeControllerToFile(const ReachabilityGame* reach, const char* filename, const UniformGrid<T>* ss, const UniformGrid<U>* is) {
  /* open file */
  std::ofstream file(filename);
  if(!file.good()) {
    std::ostringstream os;
    os << "IO: Game: Error: Unable to open file:" << filename << "'.";
    throw std::runtime_error(os.str().c_str());
  }
  /* write to file */
  file << "################################################################################" << std::endl;
  file << "#############    SCOTS: controller information    ############" << std::endl;
  file << "################################################################################" << std::endl;
  file << std::endl;
  file << "#StateSpace Parameters " << std::endl;
  file << createUniformGridParameterString(ss,"state space");
  file << std::endl;
  file << "#InputSpace Parameters" << std::endl;
  file << createUniformGridParameterString(is,"input space");
  file << std::endl;
  file << "#Controller" << std::endl;
  file << "#no of states in the transition system: " << reach->N_ << std::endl;
  file << "#no of labels in the transition system: " << reach->M_ << std::endl;
  file << "#\"idx of state\"  \"value of value function\"  \"idx of input\" " << std::endl;

  for(size_t i=0; i<reach->N_; i++) {
    if(reach->inputs_[i]!=-1 && reach->val_[i] < std::numeric_limits<double>::infinity()) {
      file << i << " " << reach->val_[i] <<" ";
      for(size_t j=0; j<reach->M_; j++)
        if(reach->inputs_[i]!=-1 && reach->inputs_[i]==j)
          file << j << " ";
      file << std::endl;
    }
  }

  file.close();
  std::cout << "controller saved to file: "<< filename << std::endl;
}


/* function: readFromFile()
* ts: TransitionSystem
* filename: name of file */
/* close function readFromFile */


}; /* close class def */
} /* close namespace */

#endif /* IO_HH_ */
