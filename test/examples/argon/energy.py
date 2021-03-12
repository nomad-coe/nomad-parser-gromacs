#!/usr/local/bin/python3
#
# Copyright The NOMAD Authors.
#
# This file is part of NOMAD.
# See https://nomad-lab.eu for further info.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#

import panedr

# Read the EDR file
path = 'ener.edr'
df = panedr.edr_to_df(path)

# The `verbose` optional parameter can be set to True to display the
# progress on stdout
df = panedr.edr_to_df(path, verbose=True)

# Get the average pressure after the first 10 ns
pressure_avg = df[u'Pressure'][df[u'Time'] > 1000].mean()
print(df[u'Pressure'][df[u'Time']])
print(df.keys())
