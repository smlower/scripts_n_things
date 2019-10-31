import caesar

caesar.driver.drive(
    snapdirs=["/ufrc/narayanan/s.lower/galaxies_project1/output/"],
    snapname="snapshot_",
    snapnums=list(range(20, 70)),
    progen_rad=False,
    extension="hdf5",
    skipran=False,
)
