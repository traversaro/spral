#pragma once

#include "WorkspaceManager.hxx"

namespace spral {
namespace ics {

/** Base class for anything that offers a factor() call */
template <typename T>
class Factorizable {
public:

virtual
void factor(T const* aval, T* lval, WorkspaceManager &memhandler) const=0;

};

} /* namespace spral */
} /* namespace ics */
