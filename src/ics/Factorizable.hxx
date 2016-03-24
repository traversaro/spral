#pragma once

#include "WorkspaceManager.hxx"

namespace spral {
namespace ics {

/** Base class for anything that offers a factor() call */
template <typename T>
class Factorizable {
public:

/** Return number of child chunks */
virtual
int get_nchild() const=0;

/** Return true if this chunk has a parent */
virtual
bool has_parent() const=0;

/** Get chunk index of parent */
virtual
int get_parent_idx() const=0;

/** Get chunk index of this node */
virtual
int get_idx() const=0;

/** Factorize all nodes in this chunk */
virtual
void factor(T const* aval, T* lval, WorkspaceManager &memhandler) const=0;

};

} /* namespace spral */
} /* namespace ics */
