/**
# Straight embedded pipe and channel geometry

Reconstruct stationary wall fractions and metrics at every stored tree
level, including neighbour storage, after adaptation. The cylindrical
construction follows `comphy-lab/bretherton-drops-bubbles`,
`src-local/embed-vof-tube.h`, commit `5fba5fa5949cd2ef1c02d387370841533dda82ee`.
The planar construction integrates the interval between the two walls.
Include after `embed.h` and, for a pipe, `axi.h`.
*/

#ifndef EMBED_CHANNEL_GEOMETRY_H
#define EMBED_CHANNEL_GEOMETRY_H

static void confined_geometry (double halfwidth)
{
#if AXI
  solid (cs, fs, halfwidth - y);
  cm_update (cm, cs, fs);
  fm_update (fm, cs, fs);
#else
  solid (cs, fs, intersection(halfwidth - y, halfwidth + y));
#endif
#if TREE
  foreach_cell() {
    double lo = y - Delta/2., hi = min(y + Delta/2., halfwidth);
#if !AXI
    lo = max(lo, -halfwidth);
#endif
    double fraction = max(hi - lo, 0.)/Delta;
    cs[] = fraction;
    fs.x[] = fraction;
    if (allocated(1)) fs.x[1] = fraction;
    double lower = y - Delta/2., upper = y + Delta/2.;
#if AXI
    fs.y[] = lower < halfwidth ? 1. : 0.;
    if (allocated(0,1)) fs.y[0,1] = upper < halfwidth ? 1. : 0.;
    double metric = hi > lo ? (sq(hi) - sq(lo))/(2.*Delta) : 0.;
    cm[] = metric;
    fm.x[] = metric;
    if (allocated(1)) fm.x[1] = metric;
    fm.y[] = lower < halfwidth ? max(lower, 1e-20) : 0.;
    if (allocated(0,1))
      fm.y[0,1] = upper < halfwidth ? max(upper, 1e-20) : 0.;
#else
    fs.y[] = lower > -halfwidth && lower < halfwidth ? 1. : 0.;
    if (allocated(0,1))
      fs.y[0,1] = upper > -halfwidth && upper < halfwidth ? 1. : 0.;
#endif
  }
  for (scalar s in {cs, fs, cm, fm}) set_dirty_stencil(s);
  restriction ({cs, fs, cm, fm});
#endif
}

#endif
