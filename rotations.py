import spglib



code = 492



p212121_gr = spglib.get_symmetry_from_database(code)
p212121_gr2 = spglib.get_spacegroup_type(code)

p212121_pg_num = spglib.get_pointgroup(p212121_gr['rotations'])[1]





