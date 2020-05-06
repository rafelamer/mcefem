# Implementation of the FEM in dimension 1 and 2  for educational purposes 
I'm teaching a course about mathematics and computer science at the Technical University of Catalonia BarcelonaTech and I started a project about the Finite Element Method in dimension 1.

One of the goals of the project was that my students could study the code and that this was as simple as possible, so I discarde more complex and powerfull programs and libraries. Then I wrote my own library for educational purposes.

As the project evolved, I wanted to add the two dimensiol method. This was done in a Final Degree Project in ESEIIAT of the student *Eduard Gómez Escandell*.

We used an adaptation Rouben Rostamian for his book *Programming Projects in C for Students of Engineering, Science and Mathematics* of triangle program by Jonathan Richard Shewchuk [https://www.cs.cmu.edu/~quake/triangle.html](https://www.cs.cmu.edu/~quake/triangle.html).

## Getting Started

When I started the project I was also interest in
1. How to read XML files in a C program.
2. How to use different linear solvers
3. How to embed Lua and Python interpreters
so I used I lot of libraries

### Prerequisites

To compile and install the library you need a Unix-like computer (I have tested it in Debian and Ubuntu) with a compiler (GCC) and the following packages or libraries
- lua5.2 and liblua5.2-dev
- xml2-dev and libxml2-dev
- python-dev
- libmatheval-dev
- libsuitesparse-dev
- libgsl-dev
- uthash-dev
- libpetsc3.10-dev

### Installing

To install the library, you have to run the following commands
```
~$ git clone https://github.com/rafelamer/mcefem.git
~$ cd mcefem
~$ make
~$ sudo make install
```
To install the programs, you have to do
```
~$ cd tfgfem
~$ make
~$ sudo make install
~$ cd ../tggmesh
~$ make
~$ sudo make install
```

## Running the demos

The folders *demos*, there are different folders with tests  and examples of use of the programs.

## Authors

* **Rafel Amer** and **Eduard Gómez Escandell**
ESEIAAT
Technical University of Catalonia BarcelonaTech
rafel.amer@upc.edu


## License

This project is licensed under the GNU Lesser General Public License.  See the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

I want to acknowledge the author of the triangle program  Jonathan Richard Shewchuk and the authors of all librarries that we have used in this project.

And also the author of the excellent book

- Programming Projects in C for Students of
Engineering, Science and Mathematics
Rouben Rostamian
SIAM
Computational Science & Engineering (2014)
ISBN 978-1-611973-49-5

