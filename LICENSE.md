# License

Copyright (C) 2025-2026 Pascal Obstetar

`treesatnemeton` is free software: you can redistribute it and/or modify
it under the terms of the **GNU General Public License version 3** as
published by the Free Software Foundation.

This program is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
Public License for more details: <https://www.gnu.org/licenses/gpl-3.0.html>.

## Relicensing note (2026-07-01)

`treesatnemeton` was previously distributed under the MIT license. It has
been relicensed to **GPL-3** for consistency across the Néméton platform,
which moved to GPL-3 after the core package `nemeton` (v0.110.0) adopted
GPL-3 (transitive dependency on the GPL-3 `biodivMapR`, spec 028).

Note: this relicensing is a **deliberate platform-consistency choice**,
not a legal obligation — `treesatnemeton` does not import the `nemeton`
core package (it is a standalone Sentinel-1/2 species classifier whose
raster outputs `nemeton` consumes as data, not as linked code). See
ADR-015 in `nemetonplateform`.
