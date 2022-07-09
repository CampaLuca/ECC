from factordb.factordb import FactorDB

value = 2465654761454364423890284910294453
fac = []

f = FactorDB(value)
f.connect()
fac += f.get_factor_list()