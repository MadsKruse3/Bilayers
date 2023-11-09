import sys
from ase.db import connect

full_lowdim_db =  connect('full_lowdim.db')
with open('COD_ICSD_tocheck_source.txt') as f:
    lines = f.readlines()

lowdim_ids = []
#for line in lines:
    #mat = full_lowdim_db.get(f'dbid={line}')
    #print('lowdim database id:', mat.id)
    #lowdim_ids.append(mat.id)

#with full_lowdim_db:
#    for i in lowdim_ids:
#        full_lowdim_db.update(i, data={'experimental': True})

#with full_lowdim_db:
#    full_lowdim_db.update(int(sys.argv[1]), data={'experimental': True})

with full_lowdim_db:
    full_lowdim_db.update(int(sys.argv[1]), data={'experimental': False})




#monolayer_db = connect('monolayers_icsd_cod.db')
#lst = []
#for i in range(344):
#    lst.append(i+1)

#with monolayer_db:
#    for i in lst:
#        monolayer_db.update(i, data={'experimental': True})

#with monolayer_db:
#    monolayer_db.update(int(sys.argv[1]), data={'experimental': True})

#with monolayer_db:
#    monolayer_db.update(int(sys.argv[1]), data={'experimental': False})

#with monolayer_db:
#    monolayer_db.update(int(sys.argv[1]), data={'experimental': None})

#with monolayer_db:
#    monolayer_db.update(int(sys.argv[1]), data={'experimental': sys.argv[2]})
