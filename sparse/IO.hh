#ifndef IO_HH_
#define IO_HH_

#include <iostream>
#include <limits>
#include "UniformGrid.hh"
#include "StaticController.hh"
#include "TransitionSystem.hh"
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
static std::string createUniformGridParameterString(const UniformGrid<U> &gs, const char* key="scots") {
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
static void readDimensionFromFile(UniformGrid<T> &ss, const char& filename, const char* key="scots") {
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
static void readMembersFromFile(UniformGrid<T> &ss, const char &filename, const char* key="scots") {

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
//          idx = std::stoull(line);
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
        file << "########### SCOTS: uniform grid parameter data                      ############" << std::endl;
        file << "################################################################################" << std::endl;
        /* append the parameter data of the uniform grid */
        file << createUniformGridParameterString(gs,key) << std::endl;
        file << "#AbstractSet" << std::endl;

        std::set<size_t>::iterator it;
        for (it=gs->abstractSet_.begin(); it != gs->abstractSet_.end(); ++it)
                file << *it << std::endl;
        file.close();
        std::cout << "Uniform grid saved to file: "<< filename << std::endl;
}

/* function: writeToFile()
 * ts: TransitionSystem
 * filename: name of file */
static void writeToFile(const TransitionSystem& ts, const char* filename) {
        /* open file */
        std::ofstream file(filename);
        if(!file.good()) {
                std::ostringstream os;
                os << "IO: TransitionSystem: Error: Unable to open file:" << filename << "'.";
                throw std::runtime_error(os.str().c_str());
        }
        /* write to file */
        file << "################################################################################" << std::endl;
        file << "###########    SCOTS: TransitionSystem information              ##########"<< std::endl;
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
                                if(ts->noPre_[k] && ts->noPre_[k][j]!=0) {
                                        for(size_t no=0; no<(size_t)ts->noPre_[k][j]; no++) {
                                                size_t pos=ts->prePointer_[k][j]+no;
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
                                if(ts->noPost_[i] && ts->noPost_[i][j]!=0) {
                                        for(size_t no=0; no<(size_t)ts->noPost_[i][j]; no++) {
                                                size_t pos=ts->postPointer_[i][j]+no;
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
static void writeControllerToFile(const StaticController &con, const char* filename) {
  /* open file */
  std::ofstream file(filename);
  if(!file.good()) {
    std::ostringstream os;
    os << "IO: Game: Error: Unable to open file:" << filename << "'.";
    throw std::runtime_error(os.str().c_str());
  }
  /* write to file */
  file << "################################################################################" << std::endl;
  file << "#############            SCOTS: controller information              ############" << std::endl;
  file << "################################################################################" << std::endl;
  file << std::endl;
  file << "#Controller" << std::endl;
  file << "#no of states in the transition system: " << con.N_ << std::endl;
  file << "#no of labels in the transition system: " << con.M_ << std::endl;
  file << "#\"idx of state\"  \"value of value function\"  \"idx of input\" " << std::endl;

  bool m=false;
  for(size_t i=0; i<con.N_; i++) {
    for(size_t j=0; j<con.M_; j++) {
      if(con.domain_[i*con.M_ + j]) {
        if(!m) {
          file << i << " " << con.val_[i] <<" ";
          m=true;
        }
        file << j << " ";
    }
    m=false;
    file << std::endl;
  }

  file.close();
  std::cout << "controller saved to file: "<< filename << std::endl;
}

/* function: writeControllerToFile()  (save the graph and grid points)
 * ss,is: UniformGrid
 * gg: Game
 * filename: name of file */
//  template<class T=std::array<double,1>, class U=std::array<double,1>>
template<class T, class U>
static void writeControllerToFile(const StaticController &con, const UniformGrid<T> &ss, const UniformGrid<U> &is, const char* filename) {
        /* open file */
        std::ofstream file(filename);
        if(!file.good()) {
                std::ostringstream os;
                os << "IO: Game: Error: Unable to open file:" << filename << "'.";
                throw std::runtime_error(os.str().c_str());
        }
        /* write to file */
        file << "################################################################################" << std::endl;
        file << "#############            SCOTS: controller information              ############" << std::endl;
        file << "################################################################################" << std::endl;
        file << std::endl;
        file << "#StateSpace Parameters " << std::endl;
        file << createUniformGridParameterString(ss,"state space");
        file << std::endl;
        file << "#InputSpace Parameters" << std::endl;
        file << createUniformGridParameterString(is,"input space");
        file << std::endl;
        file << "#Controller" << std::endl;
        file << "#no of states in the transition system: " << gg->N_ << std::endl;
        file << "#no of labels in the transition system: " << gg->M_ << std::endl;
        file << "#\"idx of state\"  \"value of value function\"  \"idx of input\" " << std::endl;

        for(size_t i=0; i<gg->N_; i++) {
                if(gg->domain_[i] && gg->val_[i] < std::numeric_limits<double>::infinity()) {
                        file << i << " " << gg->val_[i] <<" ";
                        for(size_t j=0; j<gg->M_; j++)
                                if(gg->domain_[i][j])
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
static void readFromFile(TransitionSystem& ts, const char* filename) {
        /* open file */
        std::ifstream file(filename);
        if(!file.good()) {
                std::ostringstream os;
                os << "IO: TransitionSystem: Error: Unable to open file:" << filename << "'.";
                throw std::runtime_error(os.str().c_str());
        }
        /* check if file contains safety controller */
        std::string line;
        /* should be in the second line */
        std::getline(file,line);
        std::getline(file,line);
        /* check if second line contains safety */
        if(line.find("TransitionSystem")==std::string::npos) {
                std::ostringstream os;
                os << "IO: TransitionSystem: Error: Could not read file: " << filename << ". ";
                os << "Was " << filename << " created with scots::IO::writeToFile?";
                throw std::runtime_error(os.str().c_str());
        }
        /* read N_ ,M_ and T_ parameters */
        int check=3;
        while(!file.eof()) {
                std::getline(file,line);
                if(line.find("#no of states")!=std::string::npos) {
                        check--;
                        try {
                                std::string subline(line.substr(line.find(":")+1));
                                //N_ = std::stoull(subline);
                                ts->N_=(size_t)strtol(subline.c_str(),NULL,10);
                        } catch(...) { check++; }
                }
                if(line.find("#no of labels")!=std::string::npos) {
                        check--;
                        try {
                                std::string subline(line.substr(line.find(":")+1));
                                //M_ = std::stoull(subline);
                                ts->M_=(size_t)strtol(subline.c_str(),NULL,10);
                        } catch(...) { check++; }
                }
                if(line.find("#no of transitions")!=std::string::npos) {
                        check--;
                        try {
                                std::string subline(line.substr(line.find(":")+1));
                                //N_ = std::stoull(subline);
                                ts->T_=(size_t)strtol(subline.c_str(),NULL,10);
                        } catch(...) { check++; }
                }
                if(!check)
                        break;
        }
        if(check) {
                std::ostringstream os;
                os << "IO: TransitionSystem: Error: Could not read parameters N_, M_ and T_ from file: " << filename << ". ";
                os << "Was " << filename << " created with scots::IO::writeToFile?";
                throw std::runtime_error(os.str().c_str());
        }
        if(ts->pre_!=NULL || ts->post_!=NULL) {
                free(ts->pre_);
                ts->pre_=NULL;
                for(size_t i=0; i<ts->N_; i++) {
                        free(ts->prePointer_[i]);
                        free(ts->postPointer_[i]);
                        free(ts->noPre_[i]);
                        free(ts->noPost_[i]);
                }
                free(ts->prePointer_);
                free(ts->postPointer_);
                free(ts->noPre_);
                free(ts->noPost_);
        }
        //ts->prePointer_=(state_size_t**)calloc(ts->N_,sizeof(state_size_t*));
        //ts->postPointer_=(state_size_t**)calloc(ts->N_,sizeof(state_size_t*));

        ts->prePointer_=(size_t**)calloc(ts->N_,sizeof(size_t*));
        ts->postPointer_=(size_t**)calloc(ts->N_,sizeof(size_t*));


        ts->noPre_=(unsigned short**)calloc(ts->N_,sizeof(unsigned short*));
        ts->noPost_=(unsigned short**)calloc(ts->N_,sizeof(unsigned short*));
        unsigned short** tempPointer=(unsigned short**)calloc(ts->N_,sizeof(unsigned short*));
        for(size_t i=0; i<ts->N_; i++) {
                ts->prePointer_[i]=NULL;
                tempPointer[i]=NULL;
                ts->noPre_[i]=NULL;
                ts->noPost_[i]=NULL;
                ts->postPointer_[i]=NULL;
        }
        /* first loop */
        while(!file.eof()) {
                std::getline(file,line);
                /* read transition relation */
                if(line.find("#transition relation")!=std::string::npos) {
                        std::getline(file,line);
                        for(size_t v=0; v<ts->T_; v++) {
                                size_t i,j,k;
                                std::getline(file,line);
                                std::istringstream sline(line);
                                sline >> i;
                                sline >> j;
                                sline >> k;
                                if(!ts->noPre_[k]) {
                                        ts->noPre_[k]=(unsigned short*)calloc(ts->M_,sizeof(unsigned short));

                                        //ts->prePointer_[k]=(state_size_t*)calloc(ts->M_,sizeof(state_size_t));
                                        ts->prePointer_[k]=(size_t*)calloc(ts->M_,sizeof(size_t));

                                        tempPointer[k]=(unsigned short*)calloc(ts->M_,sizeof(unsigned short));
                                }
                                ts->noPre_[k][j]++;

                                if(!ts->noPost_[i]) {
                                        ts->noPost_[i]=(unsigned short*)calloc(ts->M_,sizeof(unsigned short));

                                        //ts->postPointer_[i]=(state_size_t*)calloc(ts->M_,sizeof(state_size_t));
                                        ts->postPointer_[i]=(size_t*)calloc(ts->M_,sizeof(size_t));
                                }
                                ts->noPost_[i][j]++;
                        }
                }
        }
        file.close();
        /* integrate noPre and noPost */
        size_t sumPre=0;
        size_t sumPost=0;
        for(size_t i=0; i<ts->N_; i++) {
                for(size_t j=0; j<ts->M_; j++) {
                        if(ts->noPre_[i] && ts->noPre_[i][j]!=0) {
                                ts->prePointer_[i][j]=sumPre;
                                sumPre=sumPre+ts->noPre_[i][j];
                        }
                        if(ts->noPost_[i] && ts->noPost_[i][j]!=0) {
                                ts->postPointer_[i][j]=sumPost;
                                sumPost=sumPost+ts->noPost_[i][j];
                        }
                }
        }
        ts->pre_=(state_size_t*)malloc(sizeof(state_size_t)*ts->T_);
        file.open(filename,std::ifstream::in);
        if(!file.good()) {
                std::ostringstream os;
                os << "IO: TransitionSystem: Error: Unable to open file:" << filename << "'.";
                throw std::runtime_error(os.str().c_str());
        }
        while(!file.eof()) {
                std::getline(file,line);
                /* read transition relation */
                if(line.find("#transition relation")!=std::string::npos) {
                        std::getline(file,line);
                        for(size_t v=0; v<ts->T_; v++) {
                                size_t i,j,k;
                                std::getline(file,line);
                                std::istringstream sline(line);
                                sline >> i;
                                sline >> j;
                                sline >> k;
                                size_t pos=ts->prePointer_[k][j]+tempPointer[k][j];
                                tempPointer[k][j]++;
                                ts->pre_[pos]=i;
                        }
                }
        }
        for(size_t i=0; i<ts->N_; i++)
                free(tempPointer[i]);
        free(tempPointer);
}   /* close function readFromFile */

template<class T>
static void readFromFile(UniformGrid<T> &ss, const char* filename, const char *key="scots") {
        /* read the dimension  from file */
        readDimensionFromFile(ss,filename,key);
        ss->nofGridPoints_.resize(ss->dim_);

        /* read the UnfiormGrid members from file */
        readMembersFromFile(ss,filename,key);
        /* total number of total grid points */
        ss->N_=1;
        ss->NN_.resize(ss->dim_);
        for(size_t i=0; i<ss->dim_; i++) {
                ss->NN_[i]=ss->N_;
                ss->N_*=ss->nofGridPoints_[i];
        }
}

static void readFromFile(UniformGrid<std::vector<double> > *ss, const char* filename, const char *key="scots") {
        /* read the dimension  from file */
        readDimensionFromFile(ss,filename,key);
        ss->nofGridPoints_.resize(ss->dim_);

        /* set size of vectors */
        ss->eta_.resize(ss->dim_);;
        ss->z_.resize(ss->dim_);;
        ss->firstGridPoint_.resize(ss->dim_);

        /* read the UnfiormGrid members from file */
        readMembersFromFile(ss,filename,key);
        /* total number of total grid points */
        ss->N_=1;
        ss->NN_.resize(ss->dim_);
        for(size_t i=0; i<ss->dim_; i++) {
                ss->NN_[i]=ss->N_;
                ss->N_*=ss->nofGridPoints_[i];
        }
}

static void readControllerFromFile(StaticController con, const char* filename) {
        gg->ts_=NULL;
        /* open file */
        std::ifstream file(filename);
        if(!file.good()) {
                std::ostringstream os;
                os << "IO: Game: Error: Unable to open file:" << filename << "'.";
                throw std::runtime_error(os.str().c_str());
        }
        /* check if file contains safety controller */
        std::string line;
        /* should be in the second line */
        std::getline(file,line);
        std::getline(file,line);
        /* check if second line contains safety */
        if(line.find("controller")==std::string::npos) {
                std::ostringstream os;
                os << "IO: Game: Error: Could not read file: " << filename << ". ";
                os << "Was " << filename << " created with scots::IO::writeControllerToFile?";
                throw std::runtime_error(os.str().c_str());
        }
        /* read N_ and M_ parameters */
        int check=2;
        while(!file.eof()) {
                std::getline(file,line);
                if(line.find("#no of states")!=std::string::npos) {
                        check--;
                        try {
                                std::string subline(line.substr(line.find(":")+1));
                                //            N_ = std::stoull(subline);
                                gg->N_=(size_t)strtol(subline.c_str(),NULL,10);
                        } catch(...) { check++; }
                }
                if(line.find("#no of labels")!=std::string::npos) {
                        check--;
                        try {
                                std::string subline(line.substr(line.find(":")+1));
                                //            M_ = std::stoull(subline);
                                gg->M_=(size_t)strtol(subline.c_str(),NULL,10);
                        } catch(...) { check++; }
                }
                if(!check)
                        break;
        }
        if(check) {
                std::ostringstream os;
                os << "IO: Game: Error: Could not read parameters N_ and M_ from file: " << filename << ". ";
                os << "Was " << filename << " created with scots::IO::writeControllerToFile?";
                throw std::runtime_error(os.str().c_str());
        }
        /* init and read domain_ from file */
        gg->domain_ = (bool**)malloc(gg->N_*sizeof(bool*));
        gg->val_ = (double*)malloc(gg->N_*sizeof(double));
        for(size_t i=0; i<gg->N_; i++) {
                gg->domain_[i]=NULL;
                gg->val_[i]=std::numeric_limits<double>::infinity();
        }
        /* read domain from file */
        size_t state=0;
        size_t input=0;
        size_t value=0;
        std::getline(file,line);
        while(!file.eof()) {
                std::getline(file,line);
                std::istringstream sline(line);
                sline >> state;
                gg->domain_[state]=(bool*)calloc(gg->M_,sizeof(bool));
                sline >> value;
                gg->val_[state]=value;
                while(1) {
                        sline >> input;
                        if(!sline)
                                break;
                        gg->domain_[state][input]=1;
                }
        }
}

}; /* close class def */
} /* close namespace */

#endif /* IO_HH_ */
