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


inputString = ( "\n" 
        + "   Energies (kJ/mol)\n"
        + "           Bond            U-B    Proper Dih.  Improper Dih.      CMAP Dih.\n"
        + "    8.49566e+02    2.38781e+03    1.86323e+03    1.57601e+02   -3.65167e+02\n"
        + "          LJ-14     Coulomb-14        LJ (SR)   Coulomb (SR)      Potential\n"
        + "    9.32330e+02    1.40818e+04   -1.22302e+03   -2.00907e+04   -1.40653e+03\n"
        + "    Kinetic En.   Total Energy  Conserved En.    Temperature Pressure (bar)\n"
        + "    3.19548e+03    1.78895e+03    1.78895e+03    3.03697e+02    0.00000e+00\n"
        + "   Constr. rmsd\n"
        + "    7.91659e-07\n")

textfilter = { 
        "Proper Dih." : "Proper-Dih.",
        "Improper Dih." : "Improper-Dih.",
        "CMAP Dih." : "CMAP-Dih.",
        "LJ (SR)" : "LJ-(SR)",
        "Coulomb (SR)" : "Coulomb-(SR)",
        "Kinetic En." : "Kinetic-En.",
        "Total Energy" : "Total-Energy",
        "Conserved En." : "Conserved-En.",
        "Pressure (bar)" : "Pressure-(bar)",
        "Constr. rmsd" : "Constr.-rmsd",
        }

count = 0
for line in inputString.splitlines():
    if "Energies (kJ/mol)" in line:
        count = 0
    else:
        count += 1
    for k, v in textfilter.items():
        line = line.replace(k, v)
    if count%2 == 1:
        storedline = line.split()
    else:
        newline = line.split()
        delimeter1 = ["=" for i in range(len(newline))]
        delimeter2 = ["," for i in range(len(newline))]
        mixedline = zip(storedline, delimeter1, newline, delimeter2)
        print("".join([j + " " for i in mixedline for j in i]))


