#ifndef SAFETYGAME_HH_
#define SAFETYGAME_HH_

#include <iostream>
#include <limits>
#include "TransitionSystem.hh"
#include <queue>
#include "Game.hh"
namespace scots {
/*
 * class: SafetyGame
 *
 *
 */


/* class SafetyGame */
class SafetyGame : public Game {

public:
/* constructor */
SafetyGame(const TransitionSystem *ts) {
        ts_=ts;
        N_=ts->N_;
        M_=ts->M_;
        /* allocate memory for the controller domain */
        domain_ = (bool**)malloc(N_*sizeof(bool*));
        val_ = (double*)malloc(N_*sizeof(double));
        for(size_t i=0; i<N_; i++) {
                domain_[i]=NULL;
                val_[i]=std::numeric_limits<double>::infinity();
        }
}

/* constructor */
SafetyGame(){
}

/* destructor */
~SafetyGame() {
        if(domain_) {
                for(size_t i=0; i<N_; i++) {
                        if(domain_[i]!=NULL)
                                free(domain_[i]);
                }
                free(domain_);
        }
}

/* function: solve the safety game with respect to the set defined by
 * specification(idx)
 *
 * if specification(idx)==true -> grid point with index idx is in safe set
 * if specification(idx)==false -> grid point with index idx is not in safe set
 *
 */
template<class F>
void solve(F &specification) {
        if(!ts_) {
                std::ostringstream os;
                os << "SafetyGame.hh: Error: There is no transition system initialized.";
                throw std::runtime_error(os.str().c_str());
        }

        if(ts_->pre_==NULL) {
                std::ostringstream os;
                os << "SafetyGame.hh: Error: There is no pre list initialized. Please call firstly the function ComputeTransitionRelation.";
                throw std::runtime_error(os.str().c_str());
        }

        std::deque<size_t> heap;
        /* noLabel: keep track of the number of labels leading to safe states */
        size_t* noLabel=(size_t*)calloc(N_,sizeof(size_t));
        /* initialization */
        for(size_t i=0; i<N_; i++) {
                if(specification(i)) {
                        if(!ts_->noPost_[i]) {
                                domain_[i]=NULL;
                                heap.push_back(i);
                        }
                        else {
                                domain_[i]=(bool*)calloc(M_,sizeof(bool));
                                for(size_t j=0; j<M_; j++) {
                                        if(ts_->noPost_[i][j]!=0) {
                                                domain_[i][j]=1;
                                                noLabel[i]++;
                                        }
                                }
                        }
                }
                else {
                        domain_[i]=NULL;
                        heap.push_back(i);
                }
        }

        while(!heap.empty()) {
                size_t length=heap.size();
                /* loop over all the unsafe states found in the previous loop */
                for(size_t count=0; count<length; count++) {
                        size_t k=heap.front();
                        heap.pop_front();
                        /* loop over all the labels */
                        for(size_t j=0; j<M_; j++) {
                                size_t numOfPre=0;
                                if(ts_->noPre_[k])
                                        numOfPre=ts_->noPre_[k][j];
                                /* loop over all the source states */
                                for(size_t v=0; v<numOfPre; v++) {
                                        size_t pos=ts_->prePointer_[k][j]+v;
                                        size_t i=ts_->pre_[pos];
                                        /* check if the source states is safe */
                                        if(noLabel[i]!=0) {
                                                /* set source states i with label j as unsafe pair */
                                                if(domain_[i][j]!=0) {
                                                        domain_[i][j]=0;
                                                        noLabel[i]--;
                                                }
                                                /* add to unsafe set if the source state has no label leading to safe states */
                                                if(noLabel[i]==0) {
                                                        heap.push_back(i);
                                                        free(domain_[i]);
                                                        domain_[i]=NULL;
                                                }
                                        }
                                }
                        }
                }
        }
        for(size_t i=0; i<N_; i++) {
                if(noLabel[i]!=0)
                        val_[i]=noLabel[i];
        }
        free(noLabel);

}     /* end of solve function*/

}; /* close class def */
} /* close namespace */

#endif /* SAFETYGAME_HH_ */
