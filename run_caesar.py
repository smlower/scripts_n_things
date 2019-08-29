import caesar

caesar.driver.drive(
    snapdirs=["/ufrc/narayanan/s.lower/pd_runs/simba_m25n512/caesarsnaps_newFOF"],
    snapname="snapshot_",
    snapnums=list(range(40, 306)),
    progen_rad=True,
    extension="hdf5",
    skipran=True,
)
