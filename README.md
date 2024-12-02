mdmagma
=======

Magma code for computing with modular curves by Maarten Derickx and Andrew V. Sutherland.

# Requirements

This requires [Andrew V. Sutherlands magma package](https://github.com/AndrewVSutherland/Magma). It is known to work with the version [commit 6a5a68a8c5ba...](https://github.com/AndrewVSutherland/Magma/tree/6a5a68a8c5ba526b707f746d559c8f23cbe641b1), and should also work with later versions.

# Research Papers

The code in this project has been used for the following research papers:

[Gonality of the modular curve X<sub>1</sub>(N)](https://arxiv.org/abs/1307.5719)

[Torsion subgroups of elliptic curves over quintic and sextic number fields](http://arxiv.org/abs/1608.07549)


# Version 2

A version 2 that is more consistent and organized is under development. 
V2 will use intrinsics, user defined types and have automated tests.
This is currently done in the v2 subfolder.

## Using V2

To use V2 run 

```shell
cd v2
magma
```

And then from withing magma run

```magma
AttachSpec("mdmagma.spec");
```

### Dependencies

mdmagma depends on https://github.com/AndrewVSutherland/Magma .
Note that `AttachSpec("mdmagma.spec");` also makes sure that all external dependencies
are attached. If you want to attach the files from mdmagma without also attaching the
dependencies then you need to attach `mdmagma_internal.spec` instead. This is usefull if
you already have the dependencies of mdmagma attached by default or if you are using 
[CHIMP](https://github.com/edgarcosta/CHIMP/). There is no guarantee that mdmagma will
work with dependencies other than those listed by running:

```bash
git submodule
```




## Testing

To run the tests for this repository execute. At the momemnt it is only possible to run the tests from the tests folder. Running them from any other folder will result in an error.

```shell
cd tests
magma v2/modular_curves/test_all.m
```

# Copyright

    copyright (C) 2025 Maarten Derickx, Andrew Sutherland
    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.
    
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
    
    You should have received a copy of the GNU General Public License
    along with this program.  If not, see [http://www.gnu.org/licenses/].
