/*
 * ReachabilityGame.hh
 *
 *  created on: 08.01.2016
 *      author: rungger, rigas
 */
#ifndef REACHABILITYGAME_HH_
#define REACHABILITYGAME_HH_

#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <cassert>
#include <vector>
#include <queue>
#include <functional>

#include "TransitionSystem.hh"
#include "Heap.hh"
#include "Game.hh"


namespace scots{

/*
 * class: ReachabilityGame
 *
 * 
 * Implementation of the Dijkstra algorithm for hyper-graphs according to 
 *
 * Gallo G, Longo G, Pallottino S, Nguyen S. Directed hypergraphs and
 * applications. Discrete applied mathematics. 1993 Apr 27;42(2):177-201.
 *
 * to solve a reachability problem 
 *
 */
class ReachabilityGame : public Game {

  public: 
    /* constructor */
    ReachabilityGame(const TransitionSystem *ts) {
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
    ReachabilityGame(){}

    ~ReachabilityGame() {
      if(domain_) {
        for(size_t i=0; i<N_; i++) {
          if(domain_[i])
            free(domain_[i]);
        } 
        free(domain_);
        free(val_);
      }
    }

    /* function: solve the reachability problem with respect to target set 
     * 
     * if target(idx)==true -> grid point with index idx is in target set
     * if target(idx)==false -> grid point with index idx is not in target set
     *
     */
    template<class F>
    void solve(F &target) {
      if(!ts_) {
        std::ostringstream os;
        os << "ReachabilityGame.hh: Error: There is no transition system initialized.";
        throw std::runtime_error(os.str().c_str());
      }


      //Heap<double> heap(val_,N_);
      std::deque<size_t> heap;
      //std::vector<size_t> heap;

      /* initialize heap with indices of the target set and auxiliary matrices */
      size_t** k=(size_t**)malloc(N_*sizeof(size_t*));
      double** edgevalues=(double**)malloc(N_*sizeof(double*));
      size_t* hidx=(size_t*)malloc(N_*sizeof(size_t));
      for(size_t i=0; i<N_; i++) {
        k[i]=NULL;
        edgevalues[i]=NULL;
        if(target(i)) {
          /* nodes in the target set have value zero */
          val_[i]=0;
          domain_[i]=(bool*)calloc(M_,sizeof(bool));
          domain_[i][0]=1;
          /* put node in heap */
          //heap.push(i);
          //hidx[i]=heap.size();
          hidx[i]=1;
          heap.push_back(i);
        }
        else 
          hidx[i]=0;
      }
      /* done heap initilization */

      /* while heap non-empty */
      while(!heap.empty()) { 
        /* pop idx of domain_ with minimum val from heap */
        //size_t idx = heap.pop();
        //size_t idx = heap.back();
        //heap.pop_back();
        size_t idx = heap.front();
        heap.pop_front();
        hidx[idx]=0;
        /* is the node in the avoid set ? */
        if(!ts_->noPre_[idx] && !ts_->noPost_[idx])
          continue;
        /* loop over each label */
        for(size_t j=0; j<M_; j++) {
          /* loop over pre's associated wiht this label */
          size_t numOfPre=0;
          if(ts_->noPre_[idx])
            numOfPre=ts_->noPre_[idx][j];
          for(size_t v=0;v<numOfPre;v++) {
            size_t pos=ts_->prePointer_[idx][j]+v;
            size_t i=ts_->pre_[pos];
            /* check if the array to keep track of number of processed post's is
             * initialized */
            if(!k[i]) {
              k[i]=(size_t*)malloc(M_*sizeof(size_t));
              for(size_t l=0; l<M_; l++)
                k[i][l]=0;
            }
            /* update the number of processed posts */
            k[i][j]++;
            /* check if the edgevalues to keep track of max value in the post's is
             * initialized */
            if(!edgevalues[i]) {
              edgevalues[i]=(double*)malloc(M_*sizeof(double));
              for(size_t l=0; l<M_; l++)
                edgevalues[i][l]=0;
            }
            /* update the max value of processed posts */
            edgevalues[i][j]=(edgevalues[i][j]>=1+val_[idx] ? edgevalues[i][j] :  1+val_[idx]);
            /* check if for node i and label j all posts are processed */
            if(ts_->noPost_[i] && k[i][j]==ts_->noPost_[i][j]) {
              /* leads this label to a better value ? */
              if(val_[i]>edgevalues[i][j]) {
                //domain_[i].val=domain_[i].edgevalues[j];
                //domain_[i].idx=j;
                /* if not in heap put in heap */
                if(hidx[i]==0) {
                  //heap.push(i);
                  //hidx[i]=heap.size();
                  hidx[i]=1;
                  heap.push_back(i);
                  if(val_[i]<std::numeric_limits<double>::infinity()) {
                    size_t numOfPre2=0;
                    if(ts_->noPre_[i])
                      numOfPre2=ts_->noPre_[i][j];
                    for(size_t u=0;u<numOfPre2;u++) {
                      size_t pos2=ts_->prePointer_[i][j]+u;
                      size_t pre2=ts_->pre_[pos2];
                      k[pre2][j]--;
                    }
                  }
                } else { /* update heap */
                  //heap.heapifyup(hidx[i]-1);
                }
                val_[i]=edgevalues[i][j];
                domain_[i]=(bool*)calloc(M_,sizeof(bool));
                domain_[i][j]=1;
              }
            } /* end update of heap */
          }/* end loop over all pres of node i  under label j */
        } /* end loop over all label j */
      } /* heap is empty */
      if(edgevalues)
        for(size_t i=0;i<N_;i++) {
          if(edgevalues[i]!=NULL)
            free(edgevalues[i]);
        }
      if(k)
        for(size_t i=0;i<N_;i++) {
          if(k[i]!=NULL)
            free(k[i]);
        }
      free(k);
      free(edgevalues);
      free(hidx);
    }/* end solve() */

}; /* close class def */
} /* close namespace */


#endif
