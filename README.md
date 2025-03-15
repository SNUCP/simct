# Ciphertext Simulatable BFV

This repository provides an implementation of "Ciphertext-Simulatable HE from BFV with Randomized Evaluation"(https://eprint.iacr.org/2025/203), based on HIENAA library (https://github.com/snu-lukemin/HIENAA.jl).

Before you run the code, please make sure to install the HIENAA package.
To install, you can open the REPL and type the following command.

<pre>
<code>
]
add https://github.com/snu-lukemin/HIENAA.jl
</code>
</pre>

To run the test code, type the following command in the terminal,

<pre>
<code>
julia distance.jl
</code>
</pre>

or type the following command in the REPL, within the same path in the distance file.

<pre>
<code>
include("distance.jl")
</code>
</pre>