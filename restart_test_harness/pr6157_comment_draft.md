We ran into the same divergence in a regional run with the open boundaries on the south and north faces and the ranks split in x only, so the `Partition(4)` case is affected as well when the open boundaries are perpendicular to the partitioned direction. The current testset has its open boundaries in x, where `Partition(4)` never extends halos across them, so this case is not covered.

The testset with the boundaries moved to y (open `v`/`V` on south and north, `(Periodic, Bounded, Bounded)`, `η` varying in `y`) fails on `main` for `Partition(4)` and `Partition(2, 2)` with `extend_halos = true` and passes once the open-boundary fills in the substeps cover the barotropic kernel range. Adding the `y` orientation to the loop would cover it:

```julia
for arch in open_archs, extend_halos in (true, false), direction in (:x, :y)
```

with `build_open` choosing `(u, U)` on west/east or `(v, V)` on south/north from `direction`.
