import MDAnalysis as mda

u = mda.Universe('topol.tpr')
try:
    b = mda.Universe('topol.top')
    print(b.atoms)
except AttributeError:
    print('Not supported format.')

print(u.atoms)
ulabels = u.atoms.names
labels = [value.decode('utf-8') for value in ulabels]
print('Labels:')
print(labels)
utypes = u.atoms.types
types = [value.decode('utf-8') for value in utypes]
print('Types:')
print(types)
ures = u.residues
#ures = u.residues.resname
#res = [value.decode('utf-8') for value in ures]
print('Residue Names:')
print(ures)
#for ag in u.atoms:
#    print(ag.names)
