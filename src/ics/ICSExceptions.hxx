#pragma once

#include <exception>

namespace spral {
namespace ics {

class NotPositiveDefiniteException : public std::exception {
public:
   NotPositiveDefiniteException(int diagonal_index)
   : diagonal_index_(diagonal_index)
   {}

   virtual const char* what() const throw() {
      return "Matrix not positive definite!";
   }

private:
   int diagonal_index_;
}

} /* namespace ics */
} /* namespace spral */
