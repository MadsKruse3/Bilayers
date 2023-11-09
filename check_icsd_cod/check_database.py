from ase.db import connect

monolayer_db = connect('full_lowdim.db')
rows = monolayer_db.select('')
for row in rows:
    number_of_materials.append(row)
    experimental = row.data.experimental


    
    


