# Spherical_Path_Comparison
Spherical Path Comparison (spComparison) Package is developed for quantitatively
measuring similarity of spherical paths, especially the paleomagnetic apparent
polar wander paths (APWPs) of tectonic plates. It is powered by GMT and PmagPy.

# Work in Progress
1. IMPORTANT: Testing and debugging per-segment tests separately on angle
   difference (dif) and length dif
2. IMPORTANT: Because of complexity of the related algorithms about spherical
   surface geometry, current directional change for each segment (seg) is
   actually the dif between its previous seg and itself (if we think about it,
   this solution should be able to give close results to dif between the 1st seg
   and itself). This solution is also faster. I've already found a solution for
   calculating dif between always the 1st seg and each segment, but the
   execution time is much more than the previous solution. Maybe we need to
   balance complexity and time.
   * Complexity 1: Angle between two DIRECTIONAL geodesics (i.e. segments which
     are with DIRECTIONs in the order of poles' ages, but not neccessarily
     successive displacement ones); Angle between two geodesics (no constraints
     on their directions) could be 2 solutions at both intersections of the
     great circles that the two geodesics are on.
   * Complexity 2: For per-seg test on seg angle dif, it is hard to be
     implemented for the present ang dif solution (i.e. dif from the previous
     seg), so completing angle change of every seg compared to the 1st seg
     is neccessary now (achieved in v0.4.4).
