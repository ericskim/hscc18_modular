#ifndef ABSTRACTIONGB_HH_
#define ABSTRACTIONGB_HH_

#include <iostream>
#include "UniformGrid.hh"
#include <limits>
#include <cstring>
#include "TransitionSystem.hh"


namespace scots {
/*
 * class: AbstractionGB
 *
 * Constructs an abstraction using a growth bound according to
 * the theory in http://arxiv.org/abs/1503.03715
 *
 */
template<class stateType, class inputType>
class AbstractionGB {
  private:
    /* var: stateSpace_
     */
    UniformGrid<stateType>* stateSpace_;
    /* var: inputSpace_
     */
    UniformGrid<inputType>* inputSpace_;
    /* var: transitionsSystem_
    */
    TransitionSystem* transitionSystem_;



  public:
    /* constructor:  AbstractionGB
     *
     * initialize with the state space <UniformGrid>  and input space  <UniformGrid>
     *
     *
     */
    AbstractionGB(UniformGrid<stateType> *stateSpace, 
                  UniformGrid<inputType> *inputSpace,
                  TransitionSystem *transitionSystem)
                : stateSpace_(stateSpace),
                  inputSpace_(inputSpace),
                  transitionSystem_(transitionSystem) {
        
    }
    ~AbstractionGB(){
    }


    /* function:  computeTransitionRelation
     *
     * provide:
     * system_post(x,u)   -   computes the center of the hyper rectangle that
     *                        containes the attainable set 
     * radius_post(x,u,z) -   computes the radius of the hyper rectangle contianing
     *                        the attainabel set 
     * overflow(x,r,z)    -   check if the hyper rectangel with center x, and radius r+z 
     *                        intersects with an overflow symbol
     *
     */

    template<class F1, class F2, class F3>
    void computeTransitionRelation(F1 &system_post, F2 &radius_post, F3 &overflow, char* cmd="default") {

      /* set number of states and labels for the transistion system */
      size_t N=stateSpace_->getN();
      size_t M=inputSpace_->getN();
      transitionSystem_->N_=N;
      transitionSystem_->M_=M;
      size_t dim=stateSpace_->getDimension();
      size_t T=0;
      size_t mode=0;
      if(strcmp(cmd,"preOnly")==0)
        mode=1;
      else if(strcmp(cmd,"postOnly")==0)
        mode=2;

      /* for out of bounds check */
      std::vector<size_t> npoints(dim);
      npoints = stateSpace_->getNofGridPoints();
      std::vector<size_t> NN(dim);
      NN=stateSpace_->getNN();
      /* variables for managing the post */
      int* lb = (int*) malloc(dim*sizeof(int));
      int* ub = (int*) malloc(dim*sizeof(int));
      int* nn = (int*) malloc(dim*sizeof(int));
      stateType first=stateSpace_->getFirstGridPoint();

      /* cell radius */
      stateType eta=stateSpace_->getEta();
      stateType z=stateSpace_->getZ();
      stateType r;

      /* state and input variables */
      stateType x;
      inputType u;

      state_size_t*** boundPointer=(state_size_t***)malloc(N*sizeof(state_size_t**));
      state_size_t** prePointer=(state_size_t**)malloc(N*sizeof(state_size_t*));
      state_size_t** postPointer=(state_size_t**)malloc(N*sizeof(state_size_t*));
      short** noPre=(short**)malloc(N*sizeof(short*));
      short** noPost=(short**)malloc(N*sizeof(short*));
      short** tempPointer=(short**)malloc(N*sizeof(short*));
      for(size_t i=0; i<N; i++) {
        prePointer[i]=NULL;
        tempPointer[i]=NULL;
        noPre[i]=NULL;
        noPost[i]=NULL;
        boundPointer[i]=NULL;
        postPointer[i]=NULL;
      }

      /* first loop */
      for(size_t i=0; i<N; i++) {
        /* is x an element of the overflow symbols ? */
        stateSpace_->itox(i,x);
        if(overflow(x))
          continue;
        /* loop over all inputs */
        for(size_t j=0; j<M; j++) {
          /* current state */
          stateSpace_->itox(i,x);
          /* current input */
          inputSpace_->itox(j,u);
          /* cell radius (including measurement errors) */
          for(size_t k=0; k<dim; k++)
            r[k]=eta[k]/2.0+z[k];
          /* integrate system and radius growth bound */
          /* the result is stored in x and r */
          system_post(x,u);
          radius_post(r,u);

          /* determine the cells which intersect with the attainable set*/
          int out_of_bounds=0;
          size_t npost=1;
          for(size_t k=0; k<dim; k++) {
            lb[k] = std::lround(((x[k]-r[k]-z[k]-first[k])/eta[k]));
            ub[k] = std::lround(((x[k]+r[k]+z[k]-first[k])/eta[k]));
            if(lb[k]<0 || ub[k]>=(int)npoints[k]){
              out_of_bounds=1;
              break;
            }
            /* how many post are there */
            nn[k]=npost;
            npost*=(ub[k]-lb[k]+1);
          }
          if(out_of_bounds) {
            if(!boundPointer[i])
              boundPointer[i]=(state_size_t**)malloc(M*sizeof(state_size_t*));
            boundPointer[i][j]=NULL;
            continue;
          }
          if(!boundPointer[i])
            boundPointer[i]=(state_size_t**)malloc(M*sizeof(state_size_t*));
          boundPointer[i][j]=(state_size_t*)malloc(2*sizeof(state_size_t));
          /* calculate k and add the lower-left and upper-right indices to the array */
          for(size_t k=0; k<npost; k++) {
            int t=(int) k;
            int s;
            size_t idx=0;
            for(int l=(dim-1); l>=0; l--) {
              s=t/nn[l];
              t=t%nn[l];
              idx+=(lb[l]+s)*NN[l];
            }
            if(!noPre[idx]) {
              noPre[idx]=(short*)calloc(M,sizeof(short));
              prePointer[idx]=(state_size_t*)calloc(M,sizeof(state_size_t));
              tempPointer[idx]=(short*)calloc(M,sizeof(short));
            }
            noPre[idx][j]+=1;
            if(k==0)
              boundPointer[i][j][0]=idx;
            if(k==npost-1)
              boundPointer[i][j][1]=idx;
          }
          T+=npost;
          if(!noPost[i]) {
            noPost[i]=(short*)calloc(M,sizeof(short));
            postPointer[i]=(state_size_t*)calloc(M,sizeof(state_size_t));
          }
          noPost[i][j]+=npost;
        }
      }

      /* check overflow */
      if(T>std::numeric_limits<state_size_t>::max()) {
        std::ostringstream os;
        os << "AbstractionGB.hh: Error: Too many transitions, the data type of list could not store so many transitions. Please change the definition of state_size_t.";
        throw std::runtime_error(os.str().c_str());
      }

      /* integrate noPre and noPost */
      size_t sumPre=0;
      size_t sumPost=0;
      for(size_t i=0; i<N; i++) {
        for(size_t j=0; j<M; j++) {
          if(noPre[i] && noPre[i][j]!=0) {
            prePointer[i][j]=sumPre;
            sumPre=sumPre+noPre[i][j];
          }
          if(noPost[i] && noPost[i][j]!=0) {
            postPointer[i][j]=sumPost;
            sumPost=sumPost+noPost[i][j];
          }
        }
      }

      /* define pre/post list */
      state_size_t *pre=NULL;
      state_size_t *post=NULL;
      switch(mode) {
      case 0: {
        pre=(state_size_t*)malloc(T*sizeof(state_size_t));
        post=(state_size_t*)malloc(T*sizeof(state_size_t));
        break;
      }
      case 1: {
        pre=(state_size_t*)malloc(T*sizeof(state_size_t));
        break;
      }
      case 2: {
        post=(state_size_t*)malloc(T*sizeof(state_size_t));
        break;
      }
      }

      /* second loop */
      for(size_t i=0; i<N; i++) {
        /* is x an element of the overflow symbols ? */
        if(boundPointer[i]==NULL)
          continue;
        /* loop over all inputs */
        for(size_t j=0; j<M; j++) {
          if(boundPointer[i][j]==NULL)
            continue;
          /* extract lower-left and upper-bound points*/
          size_t k_lb=boundPointer[i][j][0];
          size_t k_ub=boundPointer[i][j][1];
          size_t npost=1;
          /* determin nn */
          for(int v=(dim-1); v>=0; v--) {
            if(v) {
              lb[v]=k_lb/NN[v];
              ub[v]=k_ub/NN[v];
              k_lb=k_lb%NN[v];
              k_ub=k_ub%NN[v];
            } else {
              lb[v]=k_lb;
              ub[v]=k_ub;
            }
            nn[v]=ub[v]-lb[v]+1;
            npost*=ub[v]-lb[v]+1;
          }
          /* recalculate nn */
          size_t temp1,temp2;
          for(size_t v=0; v<dim; v++) {
            if(v==0) {
              temp2=nn[v];
              nn[v]=1;
            }
            else {
              temp1=nn[v];
              nn[v]=temp2;
              temp2=temp1*temp2;
            }
          }
          /* for each*/
          for(size_t k=0; k<npost; k++) {
            int t=(int) k;
            int s;
            size_t idx=0;
            for(int l=(dim-1); l>=0; l--) {
              s=t/nn[l];
              t=t%nn[l];
              idx+=(lb[l]+s)*NN[l];
            }
            switch(mode) {
            case 0: {
              state_size_t posPost=postPointer[i][j]+k;
              post[posPost]=idx;
              state_size_t posPre=prePointer[idx][j]+tempPointer[idx][j];
              pre[posPre]=i;
              break;
            }
            case 1: {
              state_size_t posPre=prePointer[idx][j]+tempPointer[idx][j];
              pre[posPre]=i;
              break;
            }
            case 2: {
              state_size_t posPost=postPointer[i][j]+k;
              post[posPost]=idx;
              break;
            }
            }
            tempPointer[idx][j]++;
          }
        }
      }

      /* free tempPointer and boundPointer */
      for(size_t i=0;i<N;i++) {
        if(tempPointer[i])
          free(tempPointer[i]);
        if(boundPointer[i]) {
          for(size_t j=0;j<M;j++) {
            if(boundPointer[i][j])
              free(boundPointer[i][j]);
          }
          free(boundPointer[i]);
        }
      }
      free(tempPointer);
      free(boundPointer);

      /* transitionsystem initialization*/
      transitionSystem_->noPost_=noPost;
      transitionSystem_->noPre_=noPre;
      transitionSystem_->T_=T;
      transitionSystem_->postPointer_=postPointer;
      transitionSystem_->prePointer_=prePointer;
      transitionSystem_->pre_=pre;
      transitionSystem_->post_=post;
      /* point relation of the transition system to the computed relation */
      free(lb);
      free(ub);
      free(nn);

    }

    /* function getPre */
    std::vector<stateType> getPre(stateType x, inputType u) {
      size_t k,j;
      std::vector<stateType> pre;
      std::vector<size_t> i;
      pre.clear();
      i.clear();
      stateSpace_->xtoi(k,x);
      inputSpace_->xtoi(j,u);
      i=transitionSystem_->getPre(k,j);
      for(int v=0;v<i.size();v++) {
        stateType s;
        stateSpace_->itox(i[v],s);
        pre.push_back(s);
      }
      return pre;
    }

    /* function getPost */
    std::vector<stateType> getPost(stateType x, inputType u) {
      size_t i,j;
      std::vector<stateType> post;
      std::vector<size_t> k;
      post.clear();
      k.clear();
      stateSpace_->xtoi(i,x);
      inputSpace_->xtoi(j,u);
      k=transitionSystem_->getPost(i,j);
      for(int v=0;v<k.size();v++) {
        stateType s;
        stateSpace_->itox(k[v],s);
        post.push_back(s);
      }
      return post;
    }
}; /* close class def */
} /* close namespace */

#endif /* ABSTRACTIONGB_HH_ */
