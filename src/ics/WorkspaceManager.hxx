#pragma once

#include <except>

namespace spral {
namespace ics {

class WorkspaceManager {
public:
   WorkspaceManager(size_t maxsz)
   : workspace(nullptr), head(0), maxsz_(maxsz)
   {
      workspace = malloc(maxsz);
      if(!workspace) throw std::bad_alloc();
   }
   ~WorkspaceManager() {
      if(workspace) free(workspace);
   }

   /** Returns an uninitialized workspace of size num*sizeof(T).
    * IMPORTANT: Workspace must be released in reverse order of aquisition!
    */
   template<typename T>
   T *get(size_t num) {
      T *ptr = (T *) ( ((char *) workspace) + head );
      head += num*sizeof(T);
      return ptr;
   }

   template<typename T>
   void release(T const* ptr, size_t num) {
      head -= num*sizeof(num);
      if(ptr != (T *) ( ((char *) workspace) + head ))
         throw std::runtime_error("Memory released in incorrect order!");
   }

private:
   void *workspace;
   size_t head; // current pointer to next allocation
   size_t maxsz_; // Maximum size (in bytes)
};

} /* namespace spral */
} /* namespace ics */
