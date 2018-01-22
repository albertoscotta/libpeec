Libpeec
=======

The library provides core functions for the development of a
triangular-mesh PEEC code.
It has been written in the context of my Master's degree thesis in electrical
engineering.

Project structure
-----------------
Two main folders

* 'src', where all the source code is stored.
* 'examples', where some examples can be found.

Building and installing
-----------------------

A makefile is provided in the the src directory. This makefile has the
ability to build and install the target libpeec.a. In the src directory, type

	make

This will build the target in the project root directory, which can be installed
via the rule install, type

	make install

The default install path for the static library is /usr/local/lib, while
for include files is /usr/local/include, but both can easily be changed by
editing the variables LIBDIR and INCLUDEDIR in Makefile.

The project documentation can be generated, provided that Doxygen is correctly
installed, with the command

	doxygen

Usage
-----

Usage information can be found in the Doxygen generated documentation.
