# Class groups of imaginary quadratic points on $X_1(16)$

This folder contains the computations needed to reproduce the results of the paper: 
"Class groups of imaginary quadratic points on $X_1(16)$".


# Reproduction

Before trying to reproduce the results in this paper one should make sure that magma is installed.
The code in this directory has been verified to work with Magma V2.28-3. But later versions should
probably also work. Reproducing the results of this paper should not take much longer then 5 minutes 
and 600 MB of free RAM should be more then enough.

In order to verify the computations in the paper first download the code and make sure you are in 
the correct directory.

```bash
git clone --branch v0.2.1 https://github.com/koffie/mdmagma.git
cd mdmagma/computations/derickx-X1_16_classgroups/
```

## Reproducing single computations

In order to manually verify the computations in the paper run

```bash
magma -n order_2.m
magma -n order_5.m
```

This will verify the divisbility by 2 and 5 parts of the code.

## Reproducing everything and inpecting logs

This repository also contains a script that will reproduce all the results and make it
easy to compare the logs of this new run with my original computations. In order to do
this simply run:

```bash
./verify_all.sh
```

this requires GNU parallel to be installed and will verify all computations. The results
of all these computatios will be logged in the `logs` directory. The file 
`logs/verify_joblog.txt` contains an overview of the completed computations. The files of 
the form `logs/order_n.txt` contain the logs of the individual computations. After 
`verify_all.sh` is finished one can run
```bash
git diff
```
to see the differences between my original logs and the reproduction run.
If the reproduction was successful then the only differences reported by `git diff` 
should be differences in
1. the magma version and seed used
2. total time and total memory reported 
3. order of points
4. order of the elliptic curves used in elliptic curve chabauty.

And example of such acceptable differences can be seen at the bottom of this page

# Hardware and software specification

The computations have been run on a server at the University of Zagreb with an Intel Xeon
W-2133 CPU @ 3.60GHz with 12 cores and 64GB of RAM running Ubuntu 18.04.6 LTS and Magma V2.28-3.
The server used was overkill, and single core and 600 MB of RAM should be more than enough to 
reproduce the results.

# Acceptable differences in the logs


```diff
diff --git a/computations/derickx-X1_16_classgroups/logs/order_5.txt b/computations/derickx-X1_16_classgroups/logs/order_5.txt
index 59c30af..aa91377 100644
--- a/computations/derickx-X1_16_classgroups/logs/order_5.txt
+++ b/computations/derickx-X1_16_classgroups/logs/order_5.txt
@@ -1,4 +1,4 @@
-Magma V2.28-3     Mon Jul  7 2025 01:15:00 on euler    [Seed = 2677850875]
+Magma V2.28-3     Mon Jul  7 2025 02:07:23 on euler    [Seed = 2627718649]
 Type ? for help.  Type <Ctrl>-D to quit.
 
 Loading file "order_5.m"
@@ -37,57 +37,57 @@ Starting two cover descent + elliptic curve chabauty
 doing hk 1 0
 gamma 1
 finding points
-Time: 8.030
-Elliptic curve Elliptic Curve defined by y^2 + 6*x*y + (2*alpha^2 - 8*alpha + 
-    6)*y = x^3 + (alpha^2 - 3)*x^2 + (-3*alpha^2 + 22*alpha - 11)*x + 
-    (-13*alpha^2 + 12*alpha - 23) over K
-Time: 2.200
+Time: 8.050
+Elliptic curve Elliptic Curve defined by y^2 + 6*x*y + (4*alpha^2 + 12)*y = x^3 
++ (alpha^2 - 3)*x^2 + (-9*alpha^2 - 2*alpha - 29)*x + (-20*alpha^2 + 16*alpha - 
+    28) over K
+Time: 1.980
 starting chabauty
 13120 8
 chabauty result true <1, z^4 + (alpha^2 - alpha + 2)*z^3 + (alpha^2 - 3*alpha + 
     3)*z^2 + (alpha^2 - alpha + 2)*z + 1>
 {@ (1 : -1 : 0), (1 : 1 : 0), (0 : -1 : 1), (0 : 1 : 1) @}
 found all points
-doing hk 2 $.1 + $.3
+doing hk 2 $.1
+gamma 1
+finding points
+Time: 8.120
+Elliptic curve Elliptic Curve defined by y^2 + 6*x*y + (4*alpha^2 + 12)*y = x^3 
++ (alpha^2 - 3)*x^2 + (-9*alpha^2 - 2*alpha - 29)*x + (-20*alpha^2 + 16*alpha - 
+    28) over K
+Time: 1.970
+starting chabauty
+13120 8
+chabauty result true <1, z^4 + (alpha^2 - alpha + 2)*z^3 + (alpha^2 - 3*alpha + 
+    3)*z^2 + (alpha^2 - alpha + 2)*z + 1>
+{@ (1 : -1 : 0), (1 : 1 : 0), (0 : -1 : 1), (0 : 1 : 1) @}
+found all points
+doing hk 3 $.1 + $.2 + $.3
 gamma alpha^2 + 2
 finding points
 Time: 7.160
-Elliptic curve Elliptic Curve defined by y^2 + 6*x*y + (8*alpha^2 - 12*alpha - 
-    4)*y = x^3 + (alpha^2 - 5*alpha - 8)*x^2 + (-186*alpha^2 + 314*alpha - 
-    496)*x + (-686*alpha^2 + 1144*alpha - 2298) over K
-Time: 2.950
+Elliptic curve Elliptic Curve defined by y^2 + 6*x*y + (6*alpha^2 + 10)*y = x^3 
++ (alpha^2 + alpha - 2)*x^2 + (-180*alpha^2 + 278*alpha - 530)*x + (-859*alpha^2
+    + 1428*alpha - 2593) over K
+Time: 1.680
 starting chabauty
 4 8
 chabauty result true <alpha^2 + 2, z^4 + (alpha^2 - alpha + 2)*z^3 + (alpha^2 - 
     3*alpha + 3)*z^2 + (alpha^2 - alpha + 2)*z + 1>
-{@ (1/3 : -29/27 : 1), (1/3 : 29/27 : 1), (3 : -29 : 1), (3 : 29 : 1) @}
-found all points
-doing hk 3 $.2 + $.3
-gamma -alpha^2 - 4*alpha - 2
-finding points
-Time: 8.150
-Elliptic curve Elliptic Curve defined by y^2 + 6*x*y + (2*alpha^2 - 8*alpha + 
-    6)*y = x^3 + (alpha^2 - 3)*x^2 + (-3*alpha^2 + 22*alpha - 11)*x + 
-    (-13*alpha^2 + 12*alpha - 23) over K
-Time: 2.140
-starting chabauty
-13120 8
-chabauty result true <-alpha^2 - 4*alpha - 2, z^4 + (alpha^2 - alpha + 2)*z^3 + 
-    (alpha^2 - 3*alpha + 3)*z^2 + (alpha^2 - alpha + 2)*z + 1>
-{@ (1 : -1 : 0), (1 : 1 : 0), (0 : -1 : 1), (0 : 1 : 1) @}
+{@ (3 : -29 : 1), (3 : 29 : 1), (1/3 : -29/27 : 1), (1/3 : 29/27 : 1) @}
 found all points
 Claim 12: the rational points on Y_{16}
-[ (1 : -1 : 0), (1 : 1 : 0), (0 : -1 : 1), (0 : 1 : 1), (1 : -29 : 3), (1 : 29 :
-3), (3 : -29 : 1), (3 : 29 : 1) ]
+[ (1 : -1 : 0), (1 : 1 : 0), (0 : -1 : 1), (0 : 1 : 1), (3 : -29 : 1), (3 : 29 :
+1), (1 : -29 : 3), (1 : 29 : 3) ]
 Claim 13: table of points where divisibility by 5 could potentially fail
 (1 : -1 : 0) -3 -60 -15 2
 (1 : 1 : 0) 1 4 1 1
 (0 : -1 : 1) 1/3 -20/243 -15 2
 (0 : 1 : 1) -1 4 1 1
-(1 : -29 : 3) 29/242 -75261560815/829997587232 -2030 40
-(1 : 29 : 3) 0 0 1 1
 (3 : -29 : 1) -242/29 -628044748870/20511149 -2030 40
 (3 : 29 : 1) Infinity Infinity 1 1
+(1 : -29 : 3) 29/242 -75261560815/829997587232 -2030 40
+(1 : 29 : 3) 0 0 1 1
 Claim 14 successfully verified
 
-Total time: 261.379 [261.418] seconds, Total memory usage: 280.81MB
+Total time: 259.339 [259.357] seconds, Total memory usage: 280.81MB
```
