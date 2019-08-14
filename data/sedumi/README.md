# Problems expressed in the SEDUMI format, stored in `.MAT` files
Store in this folder any problems stored in `.MAT` files in the [SEDUMI format](https://github.com/sqlp/sedumi/blob/master/sedumi.m#L49-L92).

Some of the examples in this problem use some problems provided by Giovanni Fantuzzi (via an email on 4th of June 2019) described below:
- `ks_D2L30N240.mat`: This SDP solves in about 3 mins with `MOSEK` on my desktop, and the optimal solution is low rank (two PSD matrix variables of rank 1 or 2). I can generate larger problems with the same property if needed.
- `ks_D4L7N21.mat`: This SDP also solves in about 3 mins with `MOSEK` on my desktop, and has optimal solutions of relatively low rank. However, in this case I cannot say what the rank is. I am including this SDP because I would be interested in solving (much) larger instances of this kind.
- `ks_D4L7N28.mat`: A larger instance of the previous problem, which solves in approx. 15 mins with `MOSEK`.

All problems relate to bounding the average energy of solution of the Kuramoto-Sivashinsky equation, see this recent paper by David Goluskin and Giovanni Fantuzzi: https://doi.org/10.1088/1361-6544/ab018b.