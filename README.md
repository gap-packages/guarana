[![CI](https://github.com/gap-packages/guarana/actions/workflows/CI.yml/badge.svg)](https://github.com/gap-packages/guarana/actions/workflows/CI.yml)
[![Code Coverage](https://codecov.io/github/gap-packages/guarana/coverage.svg?branch=master&token=)](https://codecov.io/gh/gap-packages/guarana)

# The GAP 4 package `Guarana`

In this package we demonstrate the algorithmic usefulness of the
so-called Mal'cev correspondence for computations with infinite
polycyclic groups; it is a correspondence that associates to every
$\mathbb{Q}$-powered nilpotent group $H$ a unique rational nilpotent Lie algebra
$L_H$ and vice-versa. The Mal'cev correspondence was discovered by
Anatoly Mal'cev in 1951.

## Installing this package

The Guarana package is part of the standard distribution of GAP and
so normally there should be no need to install it separately.

If by any chance it is not part of your GAP distribution, then
this package is available at

     https://gap-packages.github.io/guarana/

in form of a gzipped tar archive. The standard method is to unpack the
package into the `pkg` directory of your GAP distribution. This will
create a `guarana` subdirectory.
For other non-standard options please see the chapter "Installing a GAP
Package" in the GAP Reference Manual.

Note that the GAP packages nq, Polycyclic and Polenta are needed for this
package. The example groups `GUARANA.Tr_n_O1` and `GUARANA.Tr_n_O2` need
Alnuth, and `GUARANA.F_2c_Aut1` and `GUARANA.F_3c_Aut2` need FGA; the
test suite needs both. Normally all of them are contained in your
distribution. If not, they can be obtained at

    https://www.gap-system.org/packages/


## Loading the Guarana package

If the package is not already loaded
then you have to request it explicitly.
This can be done by `LoadPackage("guarana")`.
The `LoadPackage` command is described in the GAP Reference Manual.


## Running the test suite

Once the package is installed, it is possible to check the correct
installation by running the test suite of the package. GAP exits once
the tests are finished.

    gap> TestPackage("guarana");

For more details on test files see the section "Test Files" of the GAP
Reference Manual.

If the test suite runs into an error, then please submit
a report to our issue tracker at
  <https://github.com/gap-packages/guarana/issues>

## Copyright & License

Copyright (C) 2007  Björn Assmann <bjoern@mcs.st-and.ac.uk>

Distributed under the terms of the GNU General Public License (GPL) v2

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

                https://www.gnu.org/licenses/gpl-2.0.html

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301,
USA.
