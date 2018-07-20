# Copyright 2016-2018 Berk Onat, Fawzi Mohamed
# 
#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
# 
#     http://www.apache.org/licenses/LICENSE-2.0
# 
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

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
