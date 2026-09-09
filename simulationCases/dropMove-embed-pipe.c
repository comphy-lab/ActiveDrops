/**
# Chemically active drop on the axis of an embedded pipe

An axisymmetric unit sphere initially lies on `y=0`; `x` is axial and `y`
is radial. The embedded cylindrical wall is at `y=wallHalfWidth`.
Axisymmetry excludes transverse migration and non-axisymmetric modes.
See [the shared driver](../src-local/dropMove-embed.h) for parameters,
wall chemistry, diagnostics and the resolved-gap restriction.

## Author
Vatsal Sanjay (vatsal.sanjay@comphy-lab.org)
*/

#define ACTIVE_DROP_PIPE 1
#include "dropMove-embed.h"
