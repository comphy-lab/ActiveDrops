/**
# Chemically active drop between embedded planar walls

A two-dimensional unit circle starts at `(0,dropOffset)`, between the
no-slip walls `y=+/-wallHalfWidth`. The offset allows unequal interaction
with the two walls. This is a planar drop, not a three-dimensional sphere.
See [the shared driver](../src-local/dropMove-embed.h) for parameters,
wall chemistry, diagnostics and the resolved-gap restriction.

## Author
Vatsal Sanjay (vatsal.sanjay@comphy-lab.org)
*/

#define ACTIVE_DROP_PIPE 0
#include "dropMove-embed.h"
