#ifndef FILE_OPTMEM
#define FILE_OPTMEM

/**************************************************************************/
/* File:   optmem.hh                                                      */
/* Author: Joachim Schoeberl                                              */
/* Date:   04. Apr. 97                                                    */
/**************************************************************************/

#include <mydefs.hpp>

#include "ngarray.hpp"

namespace netgen
{

/** 
    Optimized Memory allocation classes
*/

class BlockAllocator
{
private:
  ///
  unsigned size, blocks;
  ///
  void * freelist;
  ///
  NgArray<char*> bablocks;
  // Intentional leak: heap-allocated so the mutex outlives static destruction.
  // Other static destructors (e.g. via operator delete on classes using this
  // allocator) may call Free() after this BlockAllocator has been destroyed;
  // with a member std::mutex, that would crash in MSVCP140!mtx_do_lock during
  // DLL unload (static destruction order fiasco).
  mutex* block_allocator_mutex;
public:
  ///
  DLL_HEADER BlockAllocator (unsigned asize, unsigned ablocks = 100);
  ///
  DLL_HEADER ~BlockAllocator ();
  ///
  DLL_HEADER void * Alloc ();
  ///
  DLL_HEADER void Free (void * p);
};

}

#endif
