# Spherical_Path_Comparison
Spherical Path Comparison (spComparison) Package is developed for quantitatively
measuring similarity of spherical paths, especially the paleomagnetic apparent
polar wander paths (APWPs) of tectonic plates. It is powered by GMT and PmagPy.

# Work in Progress
1. IMPORTANT: Testing and debugging per-segment tests separately on angle dif
   and length dif
2. IMPORTANT: Because of complexness of the related calculations about
   spherical surface geometry, current angle dif for each segment is actually
   the dif between it and its previous segment (but if we think about it, this
   solution should be able to give close results to dif between it and the 1st
   segment); Trying to figure out dif between it and the 1st segment
   * Complexness 1: Angle between two DIRECTIONAL geodesics (i.e. segments
     which are with DIRECTIONs in the order of poles' ages, but not neccessarily
     successive displacement ones); Angle between two geodesics (no constraints
     on their directions) could be 2 solutions at both intersections of the
     great circles that the two geodesics are on.
   * Complexness 2: For per-segment test on seg angle dif, it is hard to be
     implemented for the present ang dif solution (i.e. dif from the previous
     seg), so completing ang dif of every segment from the 1st seg is
     neccessary now.
