/*
 * Heap.hh
 *
 *  created on: 08.01.2016
 *      author: rungger, rigas
 */
#ifndef HEAP_HH_
#define HEAP_HH_

#include <iostream>
#include <fstream>
#include <limits>
#include <sstream>
#include <cassert>
#include <vector>
#include <functional>

/* class to manage the heap */
template<class T>
class Heap {
  private:
    /* the heap data */
    size_t* heap_;
    size_t size_;
    /* array containing the data which defines the order */
    T *val_;
  public:
    Heap(T* val, size_t n) : val_(val) {
      heap_=(size_t*)malloc(n*sizeof(size_t));
      size_=0;
    }
    ~Heap() {
      free(heap_);
    }
    size_t size() { 
      return size_;
    }
    bool empty(void) {
      if(size_)
        return false;
      else
        return true;
    }
    void push(size_t idx) {
      heap_[size_]=idx;
      size_++;
      heapifyup(size_ - 1);
    }
    size_t pop() {
      size_t min = heap_[0];
      heap_[0] = heap_[size_ - 1];
      size_--;
      heapifydown(0);
      return min;
    }
    void heapifyup(size_t index) {    
      while((index>0) && 
            (parent(index)<=size_) &&
            (val_[heap_[parent(index)]] > val_[heap_[index]])) {
        /* swap */
        size_t tmp = heap_[parent(index)];
        heap_[parent(index)] = heap_[index];
        heap_[index] = tmp;
        index = parent(index);
      }
    }
    void heapifydown(size_t index) {     
      size_t child = left(index);
      if((child>0) && 
         (right(index)>0) &&
         (val_[heap_[child]] > val_[heap_[right(index)]])) {
        child = right(index);
      }
      if(child>0) {
        /* swap */
        size_t tmp = heap_[index];
        heap_[index] = heap_[child];
        heap_[child] = tmp;
        heapifydown(child);
      }
    }

    size_t left(size_t parent) {
      size_t i=(parent<<1)+1; // 2 * parent + 1
      return (i<size_) ? i : 0;
    }
    size_t right(size_t parent) {
      size_t i=(parent<<1) + 2; // 2 * parent + 2
      return (i<size_) ? i : 0;
    }
    size_t parent(size_t child) {
      if (child != 0) {
        size_t i = (child - 1) >> 1;
        return i;
      }
      return size_+1;
    }
};

#endif
