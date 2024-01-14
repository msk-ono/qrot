# Approximating Qubit Z-Rotation with H, S, T {#mainpage}

This C++ implementation is based on the algorithm presented in the paper [Optimal ancilla-free Clifford+T approximation of z-rotations](https://arxiv.org/abs/1403.2975).
This program is developed as part of my study.
Currenthly this program demonstrates notably poorer performance in comparison to the Haskell program [GridSynth](https://www.mathstat.dal.ca/~selinger/newsynth/) developed by the authors of the paper.

While GridSynth is distributed under the GPL license, my implementation has been developed independently without referencing the original code. Therefore, I am releasing this implementation under the MIT license.

## Building and Testing

To build this library, use Boost.Multiprecision, Boost.Program_options, and GoogleTest. Follow these commands:

```sh
$ cmake -S . -B build
$ cmake --build build
$ ctest --test-dir build
```

## Implemented Algorithms

* [Optimal ancilla-free Clifford+T approximation of z-rotations](https://arxiv.org/abs/1403.2975)
* [Fast and efficient exact synthesis of single qubit unitaries generated by Clifford and T gates](https://arxiv.org/abs/1206.5236v4)
* [Representation of Quantum Circuits with Clifford and π/8 Gates](https://arxiv.org/abs/0806.3834)

## TODO

* Enhance Speed
  * The prime factorization process currently consumes a significant amount of time
  * Implement the Quadratic Sieve algorithm
  * Use OpenMP
* Address numerical errors arising from floating-point calculations

## Run Results

The execution time of gridsynth_cpp is notably influenced by the success or failure of prime factorization.

* CPU: AMD Ryzen 7 3700X 8-Core Processor

```sh
$ time ./gridsynth_cpp pi/128 -d 10
TCount = 104
TSHTSHTHTHTHTSHTSHTSHTHTHTHTHTSHTHTSHTHTSHTSHTSHTHTHTSHTHTHTHTHTHTSHTSHTHTHTHTSHTSHTSHTHTHTSHTHTHTHTSHTSHTSHTHTHTHTSHTHTHTSHTSHTHTHTSHTSHTHTHTSHTSHTSHTSHTSHTSHTHTHTHTSHTSHTHTHTSHTSHTSHTSHTSHTSHTSHTHTSHTHTSHTHTHTSHTSHTSHTSHTSHTHTSHTSHTSHTHTSHTSHTSHTHTHTSHTSHTHTSHTSSSXS

real    0m0.288s
user    0m0.258s
sys     0m0.030s
```

```sh
$ time ./gridsynth_cpp pi/128 -d 20
TCount = 202
THTSHTHTSHTHTSHTSHTSHTHTHTSHTHTSHTHTHTSHTSHTSHTHTHTHTSHTSHTHTHTSHTSHTHTSHTHTHTSHTHTHTHTSHTSHTSHTSHTSHTSHTSHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTHTHTSHTSHTHTHTHTSHTSHTHTHTSHTHTHTHTHTHTSHTHTHTHTSHTHTSHTHTSHTSHTHTHTSHTSHTHTSHTHTHTSHTHTSHTHTHTSHTHTHTHTSHTHTHTSHTSHTHTSHTHTSHTHTSHTHTSHTHTSHTSHTHTHTSHTHTSHTSHTHTSHTHTHTHTSHTHTSHTSHTSHTSHTSHTHTSHTSHTSHTHTSHTHTSHTHTSHTSHTSHTSHTSHTSHTHTHTHTSHTHTHTHTHTHTHTSHTHTSHTSHTHTSHTHTHTHTSHTHTHTSHTHTSHTSHTHTSHTHTHTHTSHTHTSHTHTHTSHTHTSHTSHTHTSHTHTHTSHTHTHTHTSHTSHTSHTHTSHTHTHTSHWWWWW

real    0m9.688s
user    0m9.688s
sys     0m0.000s
```

```sh
$ time ./gridsynth_cpp pi/128 -d 30
TCount = 298
SHTSHTSHTSHTSHTHTHTHTHTHTSHTSHTSHTHTHTHTSHTSHTHTHTSHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTHTSHTSHTHTSHTSHTHTHTSHTSHTHTSHTHTSHTHTHTHTHTSHTHTHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTHTSHTHTSHTHTSHTHTHTHTHTSHTSHTSHTHTHTHTSHTSHTHTHTSHTSHTHTHTSHTSHTHTSHTHTSHTSHTSHTHTHTSHTSHTHTSHTHTHTHTHTHTSHTHTHTSHTSHTSHTSHTSHTSHTHTHTHTHTHTSHTSHTHTHTSHTHTHTSHTHTHTHTHTHTHTSHTSHTSHTHTSHTSHTSHTHTSHTSHTSHTHTSHTSHTSHTHTHTHTSHTSHTSHTHTSHTSHTSHTHTSHTSHTSHTSHTHTSHTSHTHTHTHTSHTHTSHTSHTSHTSHTSHTSHTSHTHTSHTHTSHTSHTHTHTSHTHTSHTHTSHTSHTSHTSHTHTSHTHTSHTHTHTHTHTSHTSHTSHTHTHTSHTSHTHTSHTSHTSHTSHTHTSHTHTHTSHTHTSHTHTSHTHTSHTHTHTHTHTSHTSHTHTSHTHTSHTHTSHTHTSHTSHTHTSHTSHTSHTHTHTHTHTSHTHTHTSHTHTSHTSHTHTSHTSHTHTSHTHTHTHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTSHTSHTSHTSHTHTSHTHTHTHTSHTHTHTHTHTSHTHTSHTHTHTHTSHTSHTHTSHTW

real    0m3.736s
user    0m3.726s
sys     0m0.010s
```

```sh
$ time ./gridsynth_cpp pi/128 -d 40
TCount = 400
TSHTHTSHTSHTHTHTSHTSHTSHTSHTSHTSHTHTHTSHTHTHTHTHTHTHTHTSHTHTSHTSHTSHTHTSHTHTHTSHTSHTSHTHTHTSHTSHTSHTSHTSHTHTSHTHTSHTSHTHTSHTHTSHTSHTHTHTSHTHTHTSHTHTHTHTHTSHTSHTHTSHTSHTSHTSHTSHTSHTSHTHTHTSHTHTHTHTSHTSHTHTSHTHTSHTSHTSHTHTSHTHTHTHTHTSHTSHTHTSHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTSHTSHTSHTHTHTSHTSHTHTHTSHTHTHTHTSHTHTHTHTHTSHTSHTSHTSHTSHTSHTHTSHTSHTHTHTSHTHTHTHTHTHTHTSHTHTSHTHTSHTSHTSHTSHTHTSHTSHTSHTSHTHTSHTHTHTHTHTHTHTHTSHTHTHTHTSHTHTHTSHTSHTHTSHTSHTHTHTHTSHTHTHTHTSHTHTHTSHTHTSHTHTSHTSHTHTSHTSHTSHTHTSHTSHTHTHTHTSHTSHTSHTHTHTHTSHTSHTHTHTHTHTHTSHTSHTHTHTSHTSHTSHTSHTHTHTSHTSHTSHTHTSHTHTHTHTHTSHTHTSHTSHTSHTHTHTHTSHTSHTHTSHTHTHTHTSHTSHTHTSHTHTSHTSHTHTSHTSHTSHTHTHTHTSHTHTSHTHTSHTSHTSHTSHTSHTHTSHTSHTHTHTHTSHTHTHTSHTHTHTSHTHTSHTSHTSHTHTHTSHTSHTSHTSHTSHTHTSHTSHTHTHTHTSHTSHTHTSHTHTHTSHTHTSHTHTSHTSHTSHTHTSHTSHTHTHTSHTHTHTHTSHTSHTHTHTHTHTSHTHTHTHTSHTSHTSHTHTHTHTSHTHTSHTHTHTSHTHTHTHTSHTSHTSHTSHTHTSHTHTSHTHTSHTHTHTSHTHTSHTHTHTHTHTSHTHTHTHTHTHTHTSHTHTSHTHTSHTSHTHTHTSHTHTSHTHTSHTHTSHTSHTSHTSHTHTSHTHTHTSHTHTSXSSSXW

real    0m4.377s
user    0m4.355s
sys     0m0.020s
```

```sh
$ time ./gridsynth_cpp pi/128 -d 50
TCount = 500
SHTHTHTHTSHTSHTSHTSHTSHTSHTSHTHTHTSHTSHTSHTHTHTHTHTSHTSHTHTSHTSHTSHTSHTSHTSHTSHTHTHTHTHTHTSHTSHTSHTSHTSHTHTSHTSHTSHTHTHTSHTSHTSHTHTHTSHTHTHTHTHTSHTSHTSHTHTSHTSHTSHTSHTSHTSHTHTHTHTSHTHTHTSHTSHTSHTSHTSHTHTSHTHTSHTHTSHTHTHTSHTSHTHTSHTSHTHTSHTHTSHTHTSHTSHTSHTSHTHTSHTSHTHTHTHTSHTHTSHTHTSHTSHTHTHTHTSHTHTSHTSHTHTHTSHTSHTHTSHTHTHTSHTSHTSHTHTSHTSHTSHTHTHTHTSHTSHTSHTSHTSHTSHTSHTHTSHTSHTSHTSHTSHTSHTHTHTHTHTHTHTSHTHTHTSHTHTHTSHTHTHTSHTHTHTSHTSHTSHTHTSHTHTSHTSHTSHTHTHTHTSHTHTSHTSHTHTHTHTHTSHTSHTSHTSHTSHTSHTHTSHTSHTHTHTHTSHTSHTSHTSHTHTHTHTSHTHTHTSHTHTSHTSHTSHTHTSHTHTHTHTHTHTSHTHTHTSHTSHTSHTHTSHTSHTHTHTHTHTHTHTHTSHTSHTHTHTSHTSHTSHTSHTSHTHTSHTHTHTHTHTSHTHTHTHTHTSHTHTSHTHTHTHTSHTHTHTHTHTSHTHTHTSHTSHTSHTSHTHTSHTHTHTSHTHTSHTSHTSHTHTHTSHTHTSHTHTSHTHTSHTSHTHTSHTHTSHTHTHTSHTSHTHTSHTHTHTSHTHTHTSHTSHTHTSHTHTSHTHTHTHTSHTSHTHTHTSHTSHTHTHTSHTSHTSHTSHTHTHTSHTSHTHTHTHTHTSHTHTHTSHTHTSHTSHTHTHTHTHTHTHTHTSHTHTHTHTSHTSHTHTHTHTHTSHTSHTHTSHTHTHTHTSHTSHTSHTHTHTHTSHTSHTHTHTSHTSHTSHTHTHTHTHTSHTHTSHTSHTSHTHTSHTSHTSHTHTHTSHTHTHTSHTHTSHTHTHTSHTSHTHTHTSHTHTSHTHTSHTHTSHTSHTSHTSHTSHTHTHTSHTHTSHTHTSHTHTSHTSHTHTHTHTSHTHTHTHTSHTHTSHTHTHTSHTSHTHTHTSHTSHTSHTSHTSHTSHTSHTHTSHTSHTHTHTHTHTSHTHTHTHTHTSHTSHTSHTHTHTSHTHTSHTHTHTSHTHTHTSHTHTHTSHTHTHTHTSHTHTSHTSHTHTHTSHTSHTHTHTSHTHTSHTSHTXSW

real    0m27.724s
user    0m27.724s
sys     0m0.000s
```

```sh
$ time ./gridsynth_cpp pi/128 -d 60
TCount = 598
HTSHTHTSHTHTHTHTSHTSHTHTSHTHTHTSHTHTSHTSHTHTSHTSHTHTHTSHTSHTHTSHTHTSHTSHTSHTSHTSHTHTHTHTHTSHTHTHTHTHTHTSHTHTSHTHTSHTSHTHTSHTSHTSHTSHTSHTHTSHTHTSHTSHTSHTHTHTHTHTHTSHTHTHTSHTHTHTSHTSHTSHTHTHTHTSHTHTSHTSHTHTHTHTSHTSHTSHTHTSHTHTSHTHTHTHTSHTHTSHTHTHTSHTHTHTHTHTSHTSHTSHTSHTSHTHTHTHTHTSHTHTHTHTHTSHTHTHTHTSHTSHTHTSHTSHTHTSHTSHTHTHTSHTSHTHTHTHTSHTSHTHTSHTHTHTSHTHTSHTHTHTHTHTHTHTHTHTSHTHTHTSHTSHTSHTSHTHTHTHTHTHTSHTSHTSHTSHTHTSHTSHTSHTSHTHTSHTHTSHTSHTHTSHTHTHTSHTHTSHTSHTHTHTSHTSHTSHTSHTHTSHTSHTSHTHTHTHTHTSHTHTSHTSHTSHTHTSHTSHTSHTSHTSHTHTSHTSHTSHTHTSHTHTSHTHTSHTHTSHTSHTSHTHTHTHTSHTSHTSHTSHTHTSHTSHTHTHTSHTSHTHTSHTSHTHTSHTSHTHTSHTSHTSHTHTSHTSHTHTSHTSHTHTHTHTSHTSHTHTHTHTSHTSHTSHTHTSHTSHTSHTHTHTHTSHTHTHTSHTHTSHTHTSHTHTSHTSHTSHTSHTSHTHTSHTHTSHTHTSHTHTSHTHTHTSHTHTSHTHTSHTHTHTHTSHTHTSHTSHTSHTHTSHTHTHTHTSHTSHTSHTHTSHTSHTHTHTSHTHTHTHTSHTHTSHTSHTHTHTHTSHTHTHTSHTSHTSHTSHTHTSHTHTSHTHTSHTHTHTHTSHTSHTSHTSHTSHTSHTSHTSHTHTHTHTSHTHTSHTSHTSHTHTSHTHTSHTSHTHTSHTSHTHTHTHTHTSHTHTHTSHTHTHTHTHTHTSHTSHTHTHTSHTHTHTSHTSHTSHTHTSHTHTSHTHTSHTHTSHTHTSHTSHTSHTHTSHTSHTHTHTHTSHTSHTHTHTSHTHTHTHTSHTHTSHTHTHTHTHTHTHTHTSHTSHTSHTHTHTSHTHTHTSHTHTSHTHTHTHTSHTSHTSHTSHTHTHTHTHTHTSHTHTSHTHTSHTHTSHTHTSHTSHTSHTHTSHTHTSHTHTHTHTHTSHTHTSHTSHTSHTHTHTHTSHTSHTHTHTSHTSHTHTHTHTSHTSHTSHTHTHTHTHTSHTHTHTSHTSHTHTSHTHTHTSHTHTSHTSHTHTSHTHTHTSHTSHTSHTHTHTHTSHTHTHTHTSHTSHTSHTHTSHTHTHTSHTSHTSHTSHTHTHTHTSHTHTSHTHTSHTHTSHTSHTSHTHTSHTSHTSHTHTHTHTSHTSHTHTSHTSHTSHTSHTHTHTHTHTHTHTSHTSHTHTSHTHTSHTSHTSHTHTSHTSHTHTSHTHTSHTSHTHTHTHTSHTSHTHTHTSHTHTHTHTSHXSXW

real    0m11.066s
user    0m11.046s
sys     0m0.020s
```
