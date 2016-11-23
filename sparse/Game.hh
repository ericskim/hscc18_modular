#ifndef GAME_HH_
#define GAME_HH_

#include <iostream>
#include <vector>
#include <string>
#include "TransitionSystem.hh"

namespace scots {

class Game {
friend class IO;
protected:
/* var: N_
 * number of states in the transition system */
size_t N_;
/* var: M_
 * number of labels in the transition system */
size_t M_;
/* var: ts_
 * pointer to the transition system */
const TransitionSystem *ts_;
/* var: domain_
 * contains the controller domain */
bool** domain_;
/* var: val_ */
double* val_;

public:
/* constructor */
Game() {
}

/* destructor */
~Game() {
}

/* function: getLabel
 * return the index of the controller input */
std::vector<size_t> getLabel(size_t idx) const {
        std::vector<size_t> label;
        label.clear();
        if (domain_[idx] && val_[idx] < std::numeric_limits<double>::infinity()) {
                for (size_t j=0; j<M_; j++)
                        if (domain_[idx][j])
                                label.push_back(j);
                return label;
        }
        else
                return label;
}

/* function:  sizeOfDomain
 * compute the size of the domain of the controller */
size_t sizeOfDomain(void) {
        size_t n=0;
        for(size_t i = 0; i < N_; i++) {
                if(domain_[i] && val_[i] < std::numeric_limits<double>::infinity())
                        n++;
        }
        return n;
}

/* function: getDomainIdx
 * copy the domain indices to an integer vector of length
 * sizeOfDomain() */
void getDomainIdx(size_t *idx) const {
        for(size_t j=0, i = 0; i < N_; i++) {
                if(domain_[i] && val_[i] < std::numeric_limits<double>::infinity()) {
                        idx[j]=i;
                        j++;
                }
        }
}

/* function: getDomainIdx
 * return the domain indices to an vector */
std::vector<size_t> getDomainIdx() const {
        std::vector<size_t> idx;
        for(size_t i = 0; i < N_; i++) {
                if(domain_[i] && val_[i] < std::numeric_limits<double>::infinity()) {
                        idx.push_back(i);
                }
        }
        return idx;
}

/* function: getValue
 * copy the values of the value function to a double vector of length
 * sizeOfDomain() */
void getValue(double *val) const {
        for(size_t j=0, i=0; i < N_; i++) {
                if(val_[i] < std::numeric_limits<double>::infinity()) {
                        val[j]=val_[i];
                        j++;
                }
        }
}

bool ifInDomain(size_t idx) {
        if(domain_[idx])
                return 1;
        else
                return 0;
}

void solve() {
}


}; /* close class def */
} /* close namespace */


#endif /* GAME_HH_ */
