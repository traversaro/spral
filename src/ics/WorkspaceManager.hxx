#pragma once

#include <stdexcept>

namespace spral {
namespace ics {

class WorkspaceManager {
public:
   WorkspaceManager(size_t maxsz)
   : workspace_(nullptr), head_(0), maxsz_(maxsz)
   {
      workspace_ = malloc(maxsz);
      if(!workspace_) throw std::bad_alloc();
   }
   ~WorkspaceManager() {
      if(workspace_) free(workspace_);
   }

   /** Returns an uninitialized workspace of size num*sizeof(T).
    * IMPORTANT: Workspace must be released in reverse order of aquisition!
    */
   template<typename T>
   T *get(size_t num) {
      T *ptr = (T *) ( ((char *) workspace_) + head_ );
      head_ += num*sizeof(T);
      return ptr;
   }

   template<typename T>
   void release(T const* ptr, size_t num) {
      head_ -= num*sizeof(T);
      if(ptr != (T *) ( ((char *) workspace_) + head_ ))
         throw std::runtime_error("Memory released in incorrect order!");
   }

private:
   void *workspace_;
   size_t head_; // current pointer to next allocation
   size_t maxsz_; // Maximum size (in bytes)
};

} /* namespace spral */
} /* namespace ics */
