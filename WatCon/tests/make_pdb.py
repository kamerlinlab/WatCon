from test_WatCon import create_water_universe


u = create_water_universe(10)


ag = u.select_atoms('all')
ag.write('waters.pdb')
