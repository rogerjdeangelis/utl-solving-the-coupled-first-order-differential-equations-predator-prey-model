%let pgm=utl-solving-the-coupled-first-order-differential-equations-predator-prey-model;

%stop_submission;

Solving the coupled first order differential equations predator prey model

Please be skeptical, I do have degrees in math and stat(but thats no guarnatee).
I may be reaching on some
of the statements below. I did some research on the net in public resources.
Unfortunately I did not capture all the links. Wiki was very useful

github
https://tinyurl.com/589cm4cb
https://github.com/rogerjdeangelis/utl-solving-the-coupled-first-order-differential-equations-predator-prey-model

https://en.wikipedia.org/wiki/Lotka%E2%80%93Volterra_equations#Applications_to_economics_and_marketing

The best way to analyze the predator prey model is numerically.
There is no closed form solutions to the predator-prey model, therefore it is difficult
to generalize and define properties.

The analysis below only deals with one specific stable model.
I do not examine the unstable cases where an eigenvalue of the Jacobi
matrix has a real part.The imaginary component of an eigenvalue is related to
oscillations. It is interesting that
      ix
     e    = cos(x) + i*sin(x)

 Coupled non-linear fist degree disfferential equations, easily
 stated with non-trivial solutions.

   dx/dt =   1.5x  - 0.1xy   Parameters   Initial Conditions
   dy/dt = 0.075xy - 1.0y     a = 1.5      x0   = 10
                              b = 0.1      y0   = 5
                              c = 0.075
   x is rats                  d = 1.0      x0   = 5
   y owls                                  y0   = 2.5

                                           x0   = 2.5
                                           y0   = 1.25
 Sections
        1 x(t) rat population over time
          y(t) owl population over time
        2 effect of initial conditions on time series
        3 Phase space diagram
        4 effect of initial conditions on phase spce diagram
        5 fit sines to the first complete ras(t) x(t)
        6 solving fo equilibrium using sympy
        7 removing time from coupled-first-order-differential-equations
          cx+by-alog(y)-dlog(x)=c

 Observations that make generalizations dificult

  1  The non-linear parametric differential equations are
     trancendental with both
     logarithmic and algebraic and has no
     closed form solutions,
  2  Solution above is oscillatory.
  3  Solutions are not perfectly periodic or sinusoidal.
  4  Amplitudes are not constant and vary with initial conditions.
  5  Period is not constat are vary with initial conditions.
  6  Without owls rats grow exponentially.
  7  Without rats owls die(decay) exponentially.


/*               _     _
 _ __  _ __ ___ | |__ | | ___ _ __ ___
| `_ \| `__/ _ \| `_ \| |/ _ \ `_ ` _ \
| |_) | | | (_) | |_) | |  __/ | | | | |
| .__/|_|  \___/|_.__/|_|\___|_| |_| |_|
|_|
*/

/**************************************************************************************************************************/
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/*            -------------------------------------------------------------------------------------------------           */
/*            |                                                                                               |           */
/*            | Numerical Analysis  (not sines and cosines)                                                   |           */
/*            |                                                                                               |           */
/*            | PROPERTIES                                 Non-Linear                           Parameters    |           */
/*            |                                                                                 a = 1.5       |           */
/*            |  1. Solutions not perfectly periodic.      dx/dt =   1.5x(t)  - 0.1x(t)y(t)     b = 0.1       |           */
/*            |  2. Amplitudes can vary over time.         dy/dt = 0.075x(t)y(t) - 1.0y(t)      c = 0.075     |           */
/*            |  3. The shape of the oscillations                                               d = 1.0       |           */
/*            |     is not a perfect sine wave                                                                |           */
/*            |                                                                       Initial Conditions      |           */
/*            |                                                                        x0   = 10              |           */
/*            |                                                                        y0   = 5               |           */
/*            |                                                                                               ||          */
/*            |                                        TIME                                                   |           */
/*            |            1.59     2.43        3.85   ----        5.92        7.22    8.06       9.48        |           */
/*            |                |       |            |                 |           |       |           |       |           */
/*            ---+--------+----^---+--------+--------+--------+--------+--------+--------+--------+--------+---           */
/* POPULATION |  0  |     1        2   |    3       |4        5       |6        7        8        9       10    POPULATION*/
/* ---------- |     |         TIME    TIME        TIME               TIME        TIME    TIME       TIME      | --------- */
/* OWLS & RATS|     |         PEAK    PEAK       TROUGH             TROUGH       PEAK    PEAK       TROUGH    |           */
/*            |     |         RATS    OWLS        RATS               OWLS        RATS    OWLS       RATS      |           */
/* OWLS       |     |        1.59     2.43        3.85               5.92        7.22    8.06       9.48      |           */
/* LAG     40 +     |          |       |            |                 |           |       |           |       +           */
/* BEHIND     |     |          |       |            |                 |           |       |           |       |           */
/* RATS       |     |         POPULATIONS           |                           POPULATIONS                     POPULATIOM*/
/* THIS       |     |        35.5      |            |    POPULATIONS  |         35.5      |           |       | PEAKS     */
/* DRIVES 35.5|-----|------   -r       +   -----  3.12--            4.81   ---   -r       +   --------+-------| 35.55 Rats*/
/* THE     35 +     |         rrr    34.23          |                 |         rrrr    34.23         |       +           */
/* RESULT 34.2|-----+--------rr+-r----ooo-----------+-----------------+---------r-+rr---oooo----------+-------| 34.23 Owls*/
/*            |     |        r | r   o |oo          |                 |        rr | r   o |oo         |       |           */
/*            |     | More  rr | rr oo | oo         |                 |        r  | r  o  | o         |       |           */
/*            |     | RATS  r  |  r o  |  o         |                 |       rr  |  r o  | oo        |       |           */
/*         30 +     | Then rr  |  r o  |  oo        |                 |       r   |  roo  |  oo       |       +           */
/*            |     | Owls r   |  rro  |   o        |                 |      rr   |  ro   |   o       |       |           */
/*            |     |     rr   |   r   |   oo       |                 |      r    |  rr   |   oo      |       |           */
/*            |     |     r    |   r   |    o OWLS  |                 |      r    |  or   |    o      |       |           */
/*            |     |     r    |  or   |     o      |                 |     rr    |  or   |    oo     |       |           */
/*         25 +     |    rr    |  o r  |     oo     |                 |     r     |  or   |     o     |       +           */
/*            |     |    r^    |  o r  |      o     |                 |    rr     |  orr  |     oo    |       |           */
/*            |     |   rr|    |  o r  |      oo    |                 |    r      | o  r  |      o    |       |           */
/*            |     |   r I    | o  r  |       o    |                 |    r      | o  r  |      oo   |       |           */
/*            |     |  rr |    | o   r |       oo   |                 |   rr      | o  r  |       o   |       |           */
/*         20 +     |  r  |    | o   r |        o   |                 |   r       |oo  rr |        o  |       +           */
/*            |     | rr  25   |oo   r |         o  |                 |  rr       |o    r |        oo |       |           */
/*            |     | r  RATS  |o    r |         oo |                 |  r        |o    r |         o |       |           */
/*            |     | r   7    |o     r|          o |                 | rr        oo    r |         oo|       |           */
/*            |     |r   OWLS  oo     r|          oo| NOT ENOUGH      | r         o     rr|          oo       |           */
/*         15 +     |r    |    o      r|           oo   RATS          |rr         o      r|           oo      +           */
/*            |     r     |    o      rr            oo  OWLS          rr         oo      r|           |o      |           */
/*            |    r|     |   oo       r            |oo  DIE          r          o|      rr           |oo     |           */
/*            |   rr|     |   o|       r            | oo             rr         oo|       r           | oo    |           */
/*            |  rr |     |  oo|       |r           |  oo           rr|         o |       rr          |  oo   |           */
/*         10 +  [10]     |  o |       |r           |   oo         rr |        oo |       |r          |   oo  +           */
/*            | INITIAL   V oo |       |rr          |    oo       rr  |        o  |       |rr         |       |           */
/*            | CONDITIONS oo  |       | rr         |     ooo    rr   |       o   |       | r         |       |           */
/*            | [5,10]    oo   |       |  r         |       ooo rr    |     ooo   |       | rr        |       | POPULATION*/
/*            |     |   ooo    |       |  rr        |         rrro    |    oo     |       |  rr       |       | TROUGHS   */
/*          5 +--[5]|oooo------+-------+---rrr------+-------rrr--oooooooooo-------+-------+---rrr-----+-------+ 4.81      */
/*            |    4.81        |       |     rrrr   |   rrrrr        4.81         |       |     rrr   |   rr  |           */
/*        3.12|---TROUGH-------+-------+--------rrrrrrrrr-----------TROUGH--------+-------+-------rrrrrrrrr---| 3.12      */
/*            |     |          |       |          3.12                |           |       |          3.12     |           */
/*            |     |          |       |          TROUGH              |           |       |         TROUGH    |           */
/*          0 +     |          |       |            |                 |           |       |           |       +           */
/*            |     |          |       |            |                 |           |       |           |       |           */
/*            ---+---------+---|---+--------+--------+--------+--------+--------+--------+--------+--------+---           */
/*               0  |      1   |   2   |    3       |4        5       |6        7 |      8|       9   |   10              */
/*                0.29        1.59   2.43         3.85               5.92       7.22     0.06        9.38                 */
/*                                                           TIME                                                         */
/*                                                           ----                                                         */
/*                                                                                                                        */
/*                                                                                                                        */
/*               --------------------------------------------------------------------------------------------             */
/*               |                     Dynamic Relation of Rat and Owl Populat|ions                         |             */
/*               |                     Preditor Prey Model (no closed form solution)                        |             */
/*               |                                                                                          |             */
/*               |        dx/dt = ax   - bxy                                     Parameters                 |             */
/*               |        dy/dt = cxy  - dy                                       a = 1.5                   |             */
/*               |                                                                b = 0.1                   |             */
/*               |        dx/dt=0 and dy/dt=0 when (x,y) = (d/c,a/b) nullclines   c = 0.075                 |             */
/*               |        when y=a/b => dx/dt=ax-ax=>y(t) =ax - ax = 0            d = 1.0                   |             */
/*               |                                                                                          |             */
/*               |        dx/dt =   1.5x  - 0.1xy                                                           |             */
/*               |        dy/dt = 0.075xy - 1.0y                                 Initial Conditions         |             */
/*               |                                                                x0   = 10                 |             */
/*               |        Definition of Parameters                                y0   = 5                  |             */
/*               |        ------------------------                                                          |             */
/*               |                                                               NullClines                 |             */
/*               |        dRat/dt = aRat   - bRatOwl                              y_ab = a/b                |             */
/*               |        dOwl/dt = cOwlRat- cOwl                                 x_dc = d/c                |             */
/*               |                                                                                          |             */
/*               |        a - rat birth rate                                PROPERIES                       |             */
/*               |                                                          ---------                       |             */
/*               |        b - rate at which owls                            Tangents at the intersection    |             */
/*               |            encounter and consume rats                    of nullclines are               |             */
/*               |            affects both populations                      perpendicular to nullclines     |             */
/*               |            (coupling term)                                                               |             */
/*               |                                                          Signs change across nullclines  |             */
/*               |        c - owls capacity to reproduce adjusted                                           |             */
/*               |            proportionally to the product of owl          rat,owl nullcliness             |             */
/*               |            and rat encounters                            steady state for owls & rats    |             */
/*               |                                                          at y=a/b dx/dt=0 a/b nullcline  |             */
/*               |        d - natural death rate of owls in                    x=d/c dy/dt=0 d/c nullcline  |             */
/*               |            the absense of prey                                                           |             */
/*               |            -c*owl (no rat food worst case?)              (0,0) is also a solution                      */
/*               |                                                          to dx/dt and dy/dt extinction                 */
/*               |                                                                                 -at                    */
/*               |                                                          No rats (x=0) owl(t)=Ae                       */
/*               |                                                          Owl extinction      at                        */
/*               |                                                          No owls rats(t) = Ae                          */
/*               |                                                          Expontial growth                              */
/*               |                                        PREY                                                            */
/*               |---------------------+------+------+------+------+------+------+------+------+------------|             */
/*               |                     0      5     10     15     20     25     30     35     40            |             */
/*               |                                                                                          |             */
/*               |        Calculation of direction is      |                                                |             */
/*               |        based on the vector direction    |                                                |             */
/*               |                                                                                          |             */
/*               |        degrees=180*atan(dydt/dxdt)/P   a/b=15  d(rats)/dt =0 (steady state               |             */
/*            40 +                                         |                     rates do not increase      + 40          */
/*               |                                                               or decrease over time)     |             */
/*               |                                  <--------------                                         |             */
/*               |                                                                                          |             */
/*            35 +                                    dx/dt  dy/dt   o (degrees)                            + 35          */
/*               |                                    -25.6, -0.0  ,0                                       |             */
/*               |                                   -----[+]-----o                                         |             */
/*               |                                 ###   X   Y    #####                                     |             */
/*            30 +                               ###   13.3,34.2     ####                                   + 30          */
/*               |            Positive slope -->##         |      o      ###                                |             */
/*               | Owls decrease Rats Increase ##      <-------- 0          ##                              |             */
/*               |                            ##           |                  ### <- Negative slpoe         |             */
/*            25 +                           *#            |                    ###    Rats decrease        + 25          */
/*     OWLS      |                          ##                                    ##   Owls increase (lag)  |             */
/*     PREDATOR  |                          ##            d/c                       ##                      | R           */
/*               | d/c=13.3 (steady state  ##                                        ##                     |             */
/*            20 + owls do not increase    ##           o  |                          ##                    + 20          */
/*               | or decrease)            #         -89   |          ^   o            ##      ^            |             */
/*               |         |               #           |   |          | 89              ##     |            |             */
/*               |         |   Initial     #           |   |          |                  #     |            |             */
/*      d     15 +         |  Population   # dxdt dyst |   0          |      dxdt dydt   o#    |            + 15          */
/*      - = 13.3 |-d/c=13--|------------   # -0.2,-11.1|-89+----------|-----  0.5,24.7,89 #----|--- d/c ----|             */
/*      c        |         |               ##          |   | 15,13.3  |        X    Y    ##    |            |             */
/*               |         V                #          V   | coexist  |      35.5,14.9 [+}                  |             */
/*            10 +          Negative Slope  ##             |      o                   ###                   + 10          */
/*               |Rats decrease Owls Increase**        ---------> 0                  ###                    |             */
/*               |                          >##            |                       ###<--Positve slope      |             */
/*               |           Initial           #    dxdt dydt  o            #####                           |             */
/*             5 +           Population [10,5] #[+]   13.5,-0.0,0      ########                             +  5          */
/*               |                                -------[+]-----------                                     |             */
/*               |                                       X  Y                                               |             */
/*               |                                     13.3,4.8                                             |             */
/*             0 +                                         |                                                +  0          */
/*               |                                   ------------->                                         |             */
/*               |                                       a/b=15                                            |              */
/*               |------------------------------------------------------------------------------------------|             */
/*               |                                                                                          |             */
/*               |                               INITIAL VALUES                                             |             */
/*               |                                                                                          |             */
/*               |    X       Y      DIRECTION     X0     Y0      COEXIST     DXDT     DYDT  ROWS           |             */
/*               |                                                                                          |             */
/*               |   13.3     4.8       -0.1      10.0    5.0    15,13.3      13.5     -0.0    29           |             */
/*               |   35.5    14.9       88.8      10.0    5.0    15,13.3       0.5     24.7   158           |             */
/*               |   13.3    34.2        0.0      10.0    5.0    15,13.3     -25.6     -0.0   243           |             */
/*               |    3.1    14.5      -89.2      10.0    5.0    15,13.3       0.2    -11.1   389           |             */
/*               |                                                                                          |             */
/*               |                                         |                                                |             */
/*               |                                         |                                                |             */
/*               ---------------------+------+------+------+------+------+------+------+------+--------------             */
/*                                    0      5     10     15     20     25     30     35     40                           */
/*                                                        a/b                                                             */
/*                                                       PREY                                                             */
/*                                                                                                                        */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*             _        ___ __     ___                   _      ___ __
/ |  _ __ __ _| |_ ___ / / |\ \   ( _ )     _____      _| |___ / / |\ \
| | | `__/ _` | __/ __| || __| |  / _ \/\  / _ \ \ /\ / / / __| || __| |
| | | | | (_| | |_\__ \ || |_| | | (_>  < | (_) \ V  V /| \__ \ || |_| |
|_| |_|  \__,_|\__|___/ | \__| |  \___/\/  \___/ \_/\_/ |_|___/ | \__| |
                       \_\  /_/                                \_\  /_/
*/

%utl_rbeginx;
parmcards4;

# Load required library
library(deSolve)
source("c:/oto/fn_tosas9x.R");

# Define the Lotka-Volterra model
lotka_volterra <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- a * X - b * X * Y
    dY <- c * X * Y - d * Y
    return(list(c(dX, dY)))
  })
}

# Set parameters
parameters <- c(a = 1.5, b = 0.1, c = 0.075, d = 1.0)
parameters
# Set initial state
state <- c(X = 10, Y = 5)

# Set time steps
times <- seq(0, 10, by = 0.01)

# Solve the differential equations
out <- ode(y = state, times = times, func = lotka_volterra, parms = parameters)

# Convert output to data frame
df <- as.data.frame(out)
pdf("d:/pdf/prey,pdf")
# Plot the results
plot(df$time, df$X, type = "l", col = "blue", xlab = "Time", ylab = "Population",
     main = "Lotka-Volterra Predator-Prey Model")
lines(df$time, df$Y, col = "red")
legend("topright", legend = c("Prey (X)", "Predator (Y)"), col = c("blue", "red"), lty = 1)

# Phase plot
plot(df$X, df$Y, type = "l", xlab = "Prey (X)", ylab = "Predator (Y)",
     main = "Phase Plot: Predator vs Prey")
fn_tosas9x(
      inp    = df
     ,outlib ="d:/sd1/"
     ,outdsn ="xtyt"
     );

;;;;
%utl_rendx;

%utl_rbeginx;
parmcards4;
library(pracma)
library(haven)
df<-read_sas("d:/sd1/want.sas7bdat")
# Find peaks for prey population
peaks <- findpeaks(df$Y, minpeakdistance = 10)
peaks[,2]<-peaks[,2] *.01
peaks
troughs <- -1*findpeaks(-df$Y)
-1*troughs
periods <- diff(peaks[,2]) *.01
periods
peaks <- findpeaks(df$X, minpeakdistance = 10)
peaks[,2]<-peaks[,2] *.01
peaks;
troughs <- -1*findpeaks(-df$X)
-1*troughs
periods <- diff(peaks[,2]) *.01
periods
;;;;
%utl_rendx(return=fromr);

options ls=100 ps=60;
proc plot data=sd1.have;
 plot x*time='r' y*time='o' /
   overlay box vref=34.2 35.5 3.12 4.8 href=2.43 8.06 1.59 7.22  5.92 3.85 9.48
   haxis=0 to 10 by 1 vaxis=0 to 40 by 5;
run;quit;


/**************************************************************************************************************************/
/*                                                                                                                        */
/* Y OWLS PEDATORS                                                                                                        */
/* ---------------                                                                                                        */
/*                                                                                                                        */
/* > peaks                                                                                                                */
/*          [,1] [,2]  [,3]  [,4]                                                                                         */
/* [1,] 34.23318 2.43   .29  5.92                                                                                         */
/* [2,] 34.23316 8.06  5.92 10.01                                                                                         */
/*                                                                                                                        */
/*                                                                                                                        */
/* > troughs <- -1*findpeaks(-df$Y)                                                                                       */
/* > -1*troughs                                                                                                           */
/*           [,1]  [,2] [,3]  [,4]                                                                                        */
/* [1,]  4.816798   .29    1  2.43                                                                                        */
/* [2,]  4.816787  5.92  243  8.06                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/* > periods y                                                                                                            */
/* [1] 5.63                                                                                                               */
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/* X RATS PREY                                                                                                            */
/* -----------                                                                                                            */
/*                                                                                                                        */
/* > peaks                                                                                                                */
/*          [,1] [,2]  [,3]  [,4]                                                                                         */
/* [1,] 35.54977 1.59     1  3.85                                                                                         */
/* [2,] 35.54956 7.22  3.85  9.48                                                                                         */
/*                                                                                                                        */
/*                                                                                                                        */
/* > troughs <- -1*findpeaks(-df$X)                                                                                       */
/* > -1*troughs                                                                                                           */
/*           [,1]. [,2]. [,3]. [,4]                                                                                       */
/* [1,]  3.123426  3.85  1.59  7.22                                                                                       */
/* [2,]  3.123434  9.48  7.22 10.01                                                                                       */
/*                                                                                                                        */
/* > periods <- diff(peaks[,2]) *.01                                                                                      */
/* > periods                                                                                                              */
/* [1] 5.63                                                                                                               */
/*                                                                                                                        */
/*                                                                                                                        */
/**************************************************************************************************************************/

%utl_rbeginx;
parmcards4;
library(pracma)
library(haven)
df<-read_sas("d:/sd1/xtyt.sas7bdat")
# Find peaks for prey population
peaks <- findpeaks(df$Y, minpeakdistance = 10)
peaks[,2]<-peaks[,2] *.01
peaks
troughs <- -1*findpeaks(-df$Y)
-1*troughs
periods <- diff(peaks[,2]) *.01
periods
peaks <- findpeaks(df$X, minpeakdistance = 10)
peaks[,2]<-peaks[,2] *.01
peaks;
troughs <- -1*findpeaks(-df$X)
-1*troughs
periods <- diff(peaks[,2]) *.01
periods
;;;;
%utl_rendx(return=fromr);

options ls=100 ps=60;
proc plot data=sd1.have;
 plot x*time='r' y*time='o' /
   overlay box vref=34.2 35.5 3.12 4.8 href=2.43 8.06 1.59 7.22  5.92 3.85 9.48
   haxis=0 to 10 by 1 vaxis=0 to 40 by 5;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/* EXTREMA                                                                                                                */
/*                                                                                                                        */
/* Y OWLS PEDATORS                                                                                                        */
/* ---------------                                                                                                        */
/*                                                                                                                        */
/* > peaks                                                                                                                */
/*          [,1] [,2]  [,3]  [,4]                                                                                         */
/* [1,] 34.23318 2.43   .29  5.92                                                                                         */
/* [2,] 34.23316 8.06  5.92 10.01                                                                                         */
/*                                                                                                                        */
/*                                                                                                                        */
/* > troughs <- -1*findpeaks(-df$Y)                                                                                       */
/* > -1*troughs                                                                                                           */
/*           [,1]  [,2] [,3]  [,4]                                                                                        */
/* [1,]  4.816798   .29    1  2.43                                                                                        */
/* [2,]  4.816787  5.92  243  8.06                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/* > periods y                                                                                                            */
/* [1] 5.63                                                                                                               */
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/* X RATS PREY                                                                                                            */
/* -----------                                                                                                            */
/*                                                                                                                        */
/* > peaks                                                                                                                */
/*          [,1] [,2]  [,3]  [,4]                                                                                         */
/* [1,] 35.54977 1.59     1  3.85                                                                                         */
/* [2,] 35.54956 7.22  3.85  9.48                                                                                         */
/*                                                                                                                        */
/*                                                                                                                        */
/* > troughs <- -1*findpeaks(-df$X)                                                                                       */
/* > -1*troughs                                                                                                           */
/*           [,1]. [,2]. [,3]. [,4]   (remove - sign i used -values to get trough)                                        */
/* [1,] - 3.123426  3.85  1.59  7.22                                                                                      */
/* [2,] - 3.123434  9.48  7.22 10.01                                                                                      */
/*                                                                                                                        */
/* > periods <- diff(peaks[,2]) *.01                                                                                      */
/* > periods                                                                                                              */
/* [1] 5.63                                                                                                               */
/*                                                                                                                        */
/*                                                                                                                        */
/*                ---+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+---       */                                                                                                                      */
/*                |                |       |            |                 |           |       |           |       |       */
/*                |                |       |            |                 |           |       |           |       |       */
/*             40 +                |       |            |                 |           |       |           |       +       */
/*                |                |       |            |                 |           |       |           |       |       */
/*                |                |       |            |                 |           |       |           |       |       */
/*                |                |       |            |                 |           |       |           |       |       */
/*                |----------------r-------+------------+-----------------+-----------r-------+-----------+-------|       */
/*             35 +                rr      |            |                 |          r|r      |           |       +       */
/*                |---------------r+------ooo-----------+-----------------+---------r-+-r----ooo----------+-------|       */
/*                |              r | r     |o           |                 |           |     o | o         |       |       */
/*                |                |     o | o          |                 |        r  |    o  |           |       |       */
/*                |             r  |  r    |  o         |                 |           |  r    |  o        |       |       */
/*             30 +                |    o  |            |                 |       r   |   o   |   o       |       +       */
/*                |            r   |       |   o        |                 |           |       |           |       |       */
/*                |                |   o   |            |                 |      r    |  r    |   o       |       |       */
/*      X & Y     |           r    |   r   |    o       |                 |           |  o    |           |       |       */
/*                |                |       |     o      |                 |     r     |       |    o      |       |       */
/*             25 +                |  o    |            |                 |           |   r   |     o     |       +       */
/*                |          r     |    r  |      o     |                 |     r     |  o    |           |       |       */
/*                |                |       |            |                 |           |       |      o    |       |       */
/*                |         r      | o     |       o    |                 |    r      |       |           |       |       */
/*                |                |       |        o   |                 |           | o  r  |       o   |       |       */
/*             20 +        r       |     r |            |                 |   r       |       |        o  |       +       */
/*                |                |       |         o  |                 |           |       |           |       |       */
/*                |       r        |o      |         o  |                 |  r        |o    r |         o |       |       */
/*                |       r        |      r|          o |                 |           |       |          o|       |       */
/*                |                o       |            |                 | r         |       |           o       |       */
/*             15 +      r         |       |           o|                 |r          o      r|           |o      +       */
/*                |     r          |       r            o                 |           |       |           |o      |       */
/*                |                o       |            |o                r          o|       |           | o     |       */
/*                |    r           |       |            | o              r|           |       r           |  o    |       */
/*                |   r           o|       |r           |  oo           r |         o |       |           |   o   |       */
/*             10 +  r           o |       |            |    o         r  |           |       |r          |    o  +       */
/*                |                |       |r           |     o        r  |        o  |       |           |       |       */
/*                |             o  |       | r          |      o      r   |       o   |       | r         |       |       */
/*                |           oo   |       |  r         |       ooo rr    |      o    |       |  r        |       |       */
/*                |         oo     |       |   r        |         rroo    |    oo     |       |   r       |       |       */
/*              5 +--ooooooo-------+-------+----rr------+-------rr----ooooooooo-------+-------+---rr------+-------+       */
/*                |                |       |      rrr   |    rrr          |           |       |     rrr   |    r  |       */
/*                |----------------+-------+---------rrrrrrrr-------------+-----------+-------+--------rrrrrrrr---|       */
/*                |                |       |            |                 |           |       |           |       |       */
/*                |                |       |            |                 |           |       |           |       |       */
/*              0 +                |       |            |                 |           |       |           |       +       */
/*                |                |       |            |                 |           |       |           |       |       */
/*                ---+--------+--------+--------+--------+--------+--------+--------+--------+--------+--------+---       */
/*                   0        1        2        3        4        5        6        7        8        9       10          */
/*                                                                                                                        */
/*                                                              TIME                                                      */
/*                                                                                                                        */
/*------------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                        */
/*            -------------------------------------------------------------------------------------------------           */
/*            |                                                                                              |            */
/*            | Numerical Analysis  (not sines and cosines)                                                  |            */
/*            |                                                                                              |            */
/*            | PROPERTIES                                 Non-Linear                           Parameters   |            */
/*            |                                                                                 a = 1.5      |            */
/*            |  1. Solutions not perfectly periodic.      dx/dt =   1.5x(t)  - 0.1x(t)y(t)     b = 0.1      |            */
/*            |  2. Amplitudes can vary over time.         dy/dt = 0.075x(t)y(t) - 1.0y(t)      c = 0.075    |            */
/*            |  3. The shape of the oscillations                                               d = 1.0      |            */
/*            |     is not a perfect sine wave                                                               |            */
/*            |                                                                       Initial Conditions     |            */
/*            |                                                                        x0   = 10             |            */
/*            |                                                                        y0   = 5              |            */
/*            |                                                                                              |            */
/*            |                                        TIME                                                  |            */
/*            |            1.59     2.43        3.85   ----        5.92        7.22    8.06       9.48       |            */
/*            |                |       |            |                 |           |       |           |      |            */
/*            ---+--------+----^---+--------+--------+--------+--------+--------+--------+--------+--------+--            */
/* POPULATION |  0  |     1        2   |    3       |4        5       |6        7        8        9       10   POPULATION */
/* ---------- |     |         TIME    TIME        TIME               TIME        TIME    TIME       TIME     | ---------  */
/* OWLS & RATS|     |         PEAK    PEAK       TROUGH             TROUGH       PEAK    PEAK       TROUGH   |            */
/*            |     |         RATS    OWLS        RATS               OWLS        RATS    OWLS       RATS     |            */
/* OWLS       |     |        1.59     2.43        3.85               5.92        7.22    8.06       9.48     |            */
/* LAG     40 +     |          |       |            |                 |           |       |           |      +            */
/* BEHIND     |     |          |       |            |                 |           |       |           |      |            */
/* RATS       |     |         POPULATIONS           |                           POPULATIONS                    POPULATIOM */
/* THIS       |     |        35.5      |                 POPULATIONS            35.5      |           |      | PEAKS      */
/* DRIVES 35.5|-----|------   -r       +   -----  3.12--            4.81   ---   -r       +   --------+------| 35.55 Rats */
/* THE     35 +     |         rrr    34.23          |                 |         rrrr    34.23         |      +            */
/* RESULT 34.2|-----+--------rr+-r----ooo-----------+-----------------+---------r-+rr---oooo----------+------| 34.23 Owls */
/*            |     |        r | r   o |oo          |                 |        rr | r   o |oo         |      |            */
/*            |     | More  rr | rr oo | oo         |                 |        r  | r  o  | o         |      |            */
/*            |     | RATS  r  |  r o  |  o         |                 |       rr  |  r o  | oo        |      |            */
/*         30 +     | Then rr  |  r o  |  oo Rats   |                 |       r   |  roo  |  oo       |      +            */
/*            |     | Owls r   |  rro  |   o        |                 |      rr   |  ro   |   o       |      |            */
/*            |     |     rr   |   r   |   oo       |                 |      r    |  rr   |   oo      |      |            */
/*            |     |     r    |   r   |    o       |                 |      r    |  or   |    o      |      |            */
/*            |     |     r    |  or   |     o      |                 |     rr    |  or   |    oo     |      |            */
/*         25 +     |    rr    |  o r  |     oo     |                 |     r     |  or   |     o     |      +            */
/*            |     |    r     |  o r  |      o     |                 |    rr     |  orr  |     oo    |      |            */
/*            |     |   rr     |  o r  |      oo    |                 |    r      | o  r  |      o    |      |            */
/*            |     |   r      | o  r  |       o    |                 |    r      | o  r  |      oo   |      |            */
/*            |     |  rr      | o   r |       oo   |                 |   rr      | o  r  |       o   |      |            */
/*         20 +     |  r       | o   r |        o   |                 |   r       |oo  rr |        o  |      +            */
/*            |     | rr       |oo   r |         o  |                 |  rr       |o    r |        oo |      |            */
/*            |     | r        |o    r |         oo |                 |  r        |o    r |         o |      |            */
/*            |     | r        |o     r|          o |                 | rr        oo    r |         oo|      |            */
/*            |     |r         oo     r|          oo| NOT ENOUGH      | r         o     rr|          oo      |            */
/*         15 +     |r         o      r|           oo   RATS          |rr         o      r|           oo     +            */
/*            |     r          o      rr            oo  OWLS          rr         oo      r|           |o     |            */
/*            |    r|         oo       r            |oo  DIE          r          o|      rr           |oo    |            */
/*            |   rr|         o|       r            | oo             rr         oo|       r           | oo   |            */
/*            |  rr |        oo|       |r           |  oo           rr|         o |       rr          |  oo  |            */
/*         10 +  [10]        o |       |r           |   oo         rr |        oo |       |r          |   oo +            */
/*            | INITIAL     oo |       |rr          |    oo       rr  |        o  |       |rr         |      |            */
/*            | CONDITIONS oo  |       | rr         |     ooo    rr   |       o   |       | r         |      |            */
/*            | [5,10]    oo   |       |  r         |       ooo rr    |     ooo   |       | rr        |      | POPULATION */
/*            |     |   ooo    |       |  rr        |         rrro    |    oo     |       |  rr       |      | TROUGHS    */
/*          5 +--[5]|oooo------+-------+---rrr------+-------rrr--oooooooooo-------+-------+---rrr-----+------+ 4.81       */
/*            |    4.81        |       |     rrrr   |   rrrrr        4.81         |       |     rrr   |   rr |            */
/*        3.12|---TROUGH-------+-------+--------rrrrrrrrr-----------TROUGH--------+-------+-------rrrrrrrrr--| 3.12       */
/*            |     |          |       |          3.12                |           |       |          3.12    |            */
/*            |     |          |       |          TROUGH              |           |       |         TROUGH   |            */
/*          0 +     |          |       |            |                 |           |       |           |      +            */
/*            |     |          |       |            |<-~RATS PERIOD  9.38-3.85 = 5.63 --------------->|      |            */
/*            |     |          |       |            |                 |           |       |           |      |            */
/*            |     |<-------~OWLS PERIOD 5.92-.029 = 5.89 --------- >|           |       |           |      |            */
/*            |     |          |       |            |                 |           |       |           |      |            */
/*            |     |          |       |            |                 |           |       |           |      |            */
/*            ---+---------+---|---+--------+--------+--------+--------+--------+--------+--------+--------+--            */
/*               0  |      1   |   2   |    3       |4        5       |6        7 |      8|       9   |   10              */
/*                0.29        1.59   2.43         3.85               5.92       7.22     0.06        9.48                 */
/*                                                           TIME                                                         */
/*                                                           ----                                                         */
/*                                                                                                                        */
/**************************************************************************************************************************/


/*___          __  __           _               _        _______     ___                   _      _______
|___ \    ___ / _|/ _| ___  ___| |_   _ __ __ _| |_ ___ / / _ \ \   ( _ )     _____      _| |___ / / _ \ \
  __) |  / _ \ |_| |_ / _ \/ __| __| | `__/ _` | __/ __| | | | | |  / _ \/\  / _ \ \ /\ / / / __| | | | | |
 / __/  |  __/  _|  _|  __/ (__| |_  | | | (_| | |_\__ \ | |_| | | | (_>  < | (_) \ V  V /| \__ \ | |_| | |
|_____|  \___|_| |_|  \___|\___|\__| |_|  \__,_|\__|___/ |\___/| |  \___/\/  \___/ \_/\_/ |_|___/ |\___/| |
                                                        \_\   /_/                                \_\   /_/
*/

%utl_rbeginx;
parmcards4;

# Load required library
library(deSolve)
source("c:/oto/fn_tosas9x.R");

# Define the Lotka-Volterra model
lotka_volterra <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- a * X - b * X * Y
    dY <- c * X * Y - d * Y
    return(list(c(dX, dY)))
  })
}

# Set parameters
parameters <- c(a = 1.5, b = 0.1, c = 0.075, d = 1.0)
parameters
# Set initial state
state <- c(X = 5, Y =2.5 )

# Set time steps
times <- seq(0, 10, by = 0.01)

# Solve the differential equations
out <- ode(y = state, times = times, func = lotka_volterra, parms = parameters)

# Convert output to data frame
df <- as.data.frame(out)
pdf("d:/pdf/prey,pdf1")
# Plot the results
plot(df$time, df$X, type = "l", col = "blue", xlab = "Time", ylab = "Population",
     main = "Lotka-Volterra Predator-Prey Model")
lines(df$time, df$Y, col = "red")
legend("topright", legend = c("Prey (X)", "Predator (Y)"), col = c("blue", "red"), lty = 1)

# Phase plot
plot(df$X, df$Y, type = "l", xlab = "Prey (X)", ylab = "Predator (Y)",
     main = "Phase Plot: Predator vs Prey")
fn_tosas9x(
      inp    = df
     ,outlib ="d:/sd1/"
     ,outdsn ="xtyt1"
     );
;;;;
%utl_rendx;

options validvarname=upcase;
proc sql;
 create
    table sd1.rats as
 select
    l.time
   ,l.x  as rats10
   ,r.x  as rats5
 from
   sd1.xtyt as l, sd1.xtyt1 as r
 where
   l.time = r.time
;quit;

%utl_rbeginx;
parmcards4;
library(pracma)
library(haven)
df<-read_sas("d:/sd1/rats.sas7bdat")
# Find peaks for prey population
peaks <- findpeaks(df$RATS10, minpeakdistance = 10)
peaks[,2]<-peaks[,2] *.01
peaks
troughs <- -1*findpeaks(-df$RATS10)
-1*troughs
periods <- diff(peaks[,2]) *.01
periods

peaks <- findpeaks(df$RATS5, minpeakdistance = 10)
peaks[,2]<-peaks[,2] *.01
peaks;
troughs <- -1*findpeaks(-df$RATS5)
-1*troughs
periods <- diff(peaks[,2]) *.01
periods
;;;;
%utl_rendx(return=fromr);


libname sd1 "d:/sd1";

options ls=110 ps=55;
proc plot data=sd1.rats;
  plot rats10*time='t' rats5*time='5' / overlay box
   vref=35.5 56.5 35.5 56.5 href=1.59 2.13 7.22 8.60 3.85 4.38 9.48 ;
run;quit;


/**************************************************************************************************************************/
/*                                                                                                                        */
/*  EXTREMA                                                                                                               */
/*                                                                                                                        */
/*  RATS10 PREY  RATS(Time=0)=10                                                                                          */
/*  ----------------------------                                                                                          */
/*                                                                                                                        */
/*  > peaks                                                                                                               */
/*           [,1] [,2]  [,3]  [,4]                                                                                        */
/*  [1,] 35.54977 1.59     1  3.85                                                                                        */
/*  [2,] 35.54956 7.22  3.85  9.48                                                                                        */
/*                                                                                                                        */
/*                                                                                                                        */
/*  > troughs <- -1*findpeaks(-df$X)                                                                                      */
/*  > -1*troughs                                                                                                          */
/*            [,1]. [,2]. [,3]. [,4]                                                                                      */
/*  [1,]  3.123426  3.85  1.59  7.22                                                                                      */
/*  [2,]  3.123434  9.48  7.22 10.01                                                                                      */
/*                                                                                                                        */
/*  > periods <- diff(peaks[,2]) *.01                                                                                     */
/*  > periods                                                                                                             */
/*  [1] 5.63                                                                                                              */
/*                                                                                                                        */
/*  RATS5  RATS(Time=5)                                                                                                   */
/*  --------------------                                                                                                  */
/*                                                                                                                        */
/*  > peaks <- findpeaks(df$RATS5, minpeakdistance = 10)                                                                  */
/*  > peaks[,2]<-peaks[,2] *.01                                                                                           */
/*  > peaks;                                                                                                              */
/*           [,1] [,2] [,3] [,4]                                                                                          */
/*  [1,] 56.49897 2.13    1  438                                                                                          */
/*  [2,] 56.49869 8.60  438 1001                                                                                          */
/*                                                                                                                        */
/*  > troughs <- -1*findpeaks(-df$RATS5)  (remove - sign i used -values to get trough)                                    */
/*  > -1*troughs                                                                                                          */
/*             [,1] [,2] [,3] [,4]                                                                                        */
/*  [1,] -0.8712222  438  213  860                                                                                        */
/*                                                                                                                        */
/*  > periods <- diff(peaks[,2]) *.01                                                                                     */
/*  > periods                                                                                                             */
/*  [1] 0.0647                                                                                                            */
/*                                                                                                                        */
/*------------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                        */
/*           -+----------------------------------------------------------------------------------------------------|      */
/*           |                                                                                                     |      */
/*           |  Time Series of Rat Population (x) for two dets of intial conditions                                |      */
/*           |  Owl population is hidden                                                                           |      */
/*           |                                                                                                     |      */
/*           |  Staring populations are determative                                                                |      */
/*           |                                                                                                     |      */
/*           |  Alplitudes and Periods                                                                             |      */
/*           |                                                                                                     |      */
/*           |  dx/dt =   1.5x(t)  - 0.1x(t)y(t)                                                                   |      */
/*           |  dy/dt = 0.075x(t)y(t) - 1.0y(t)                                                                    |      */
/*           |                                                                                                     |      */
/*           |  Initial Conditions                                                                                 |      */
/*           |   x0   = 10                                                                                         |      */
/*           |   y0   = 5                                                                                          |      */
/*           |                                                                                                     |      */
/*           |   x0   = 5                                                                                          |      */
/*           |   y0   = 2.5                                                                                        |      */
/*           |                                                                                                     |      */
/*           -+---------+---------+---------+---------+---------+---------+---------+---------+---------+--------+-|      */
/*           |0         1         2         3         4         5         6         7         8         9        10|      */
/*           |                                                                                                     |      */
/*           |                                                                                                     |      */
/*        60 +                 RATS(0)=5                                                        RATS(0)=5          |      */
/*           |----------------+----5 <- 56.6 --------    56.5  -----------------------+-------------5-- 56.6 +-----|      */
/*           |                    555                                                              555             |      */
/*           |                    5|55        Plot symbol for RATS                                55|5             |      */
/*           |                   55| 5          t = (x0,y0) = (10,5)                              5 |55            |      */
/*           |                   5 | 5          5 = (x0,y0) = (5,,2.5)                            5 | 5            |      */
/*        50 +                   5 | 5                                                           55 | 5            +      */
/*           |                  5  | 55       Owl population is hidden                           5  | 5            |      */
/*           |                  5  |  5                                                          5  | 5            |      */
/*           |                  5  |  5                                                         55  |  5           |      */
/*           |                 5   |  5                                                         5   |  5           |      */
/*           |                 5   |  5                                                         5   |  5           |      */
/*           |                 5   |  5                                                        55   |  5           |      */
/*        40 +                55   |  55                                                       5    |  5           +      */
/*           |                5    |   5                                                       5    |  5           |      */
/*           |                5    |   5                                                      55    |   5          |      */
/*      RATS |--RATS(0)=10-> ttt---+---5----------35.5----------------  RATS(0)=10 ->ttt------5-----+---5----+-----|      */
/*           |              tt|tt  |   5                                            tt|tt     5     |   5    |     |      */
/*           |             t 5| t  |   5                                           tt | tt   55     |   5    |     |      */
/*           |            tt5 |  t |   5                                          tt  |  t   5      |   5    |     |      */
/*        30 +           tt 5 |  t |    5                                         t   |  tt  5      |   5    |     +      */
/*           |           t 55 |   t|    5                                        tt   |   t 5       |   55   |     |      */
/*           |          tt 5  |   t|    5                                       tt    |   tt5       |    5   |     |      */
/*           |         tt  5  |   tt    5                                       t     |    t5       |    5   |     |      */
/*           |         t  5   |    t    5                                      tt     |    t        |    5   |     |      */
/*           |        tt  5   |    t    5                                     tt      |   5tt       |    5   |     |      */
/*           |       tt  5    |    |t    5                                    t       |   5 t       |    5   |     |      */
/*        20 +       t   5    |    |t    5                                   tt       |  55 tt      |     5  |     +      */
/*           |      tt  5     |    |tt   5                                  tt        |  5   t      |     5  |     |      */
/*           |     tt  55     |    | t   5                                 tt         | 55   t      |     5  |     |      */
/*           |    tt   5      |    | tt  55                               tt          |55    tt     |     5  |     |      */
/*           |   tt   55      |    |  t   5                              tt           |5      t     |     55 |     |      */
/*           |  tt   55       |    |  tt  5                             tt            5       tt    |      5 |     |      */
/*  rats(0)  | tt   55        |    |   t  55                           tt            55        tt   |      5 |     |      */
/*      = 10 +tt   55         |    |   tt  5                          tt            55|         t   |      55|     +      */
/*           |    55          |    |    tt 5                        ttt           555 |          t  |       5|     |      */
/*  rats(0)  |  555           |    |     tt55                     ttt            55   |          ttt|       55     |      */
/*       = 5 +555             |    |      ttt5                 tttt            555    |            tt        5     |      */
/*           |                |    |        ttttt    |    tttttt            5555      |             |tttt    |55  t|      */
/*           |                |    |          55tt- 3.1 -tt             55555         |             |   ttttttttttt|      */
/*           |                |    |            55555  --- 0.67 -- 55555              |             |        |   55|      */
/*         0 +                |    |                                                  |             |        |            */
/*           |                |    |                 |    |                           |             |        |     |      */
/*           |                |    |                 |    |                           |             |        |     |      */
/*           |                |    |                 |    |                           |             |        |     |      */
/*           |               1.59 2.13             3.85  4.38                        7.22           |       9.48   |      */
/*           -+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+-      */
/*            0         1         2         3         4         5         6         7         8         9        10       */
/*                                                                                                                        */
/*                                                                    TIME                                                */
/*                                                                                                                        */
/**************************************************************************************************************************/


/*____         _                                                      _       _
|___ /   _ __ | |__   __ _ ___  ___   ___ _ __   __ _  ___ ___  _ __ | | ___ | |_
  |_ \  | `_ \| `_ \ / _` / __|/ _ \ / __| `_ \ / _` |/ __/ _ \| `_ \| |/ _ \| __|
 ___) | | |_) | | | | (_| \__ \  __/ \__ \ |_) | (_| | (_|  __/| |_) | | (_) | |_
|____/  | .__/|_| |_|\__,_|___/\___| |___/ .__/ \__,_|\___\___|| .__/|_|\___/ \__|
        |_|                              |_|                   |_|

*/

%utl_rbeginx;
parmcards4;

# Load required library
library(deSolve)
source("c:/oto/fn_tosas9x.R");

# Define the Lotka-Volterra model
lotka_volterra <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dX <- a * X - b * X * Y
    dY <- c * X * Y - d * Y
    return(list(c(dX, dY)))
  })
}

# Set parameters
parameters <- c(a = 1.5, b = 0.1, c = 0.075, d = 1.0)
parameters
# Set initial state
state <- c(X = 10, Y = 5)

# Set time steps
times <- seq(0, 10, by = 0.01)

# Solve the differential equations
out <- ode(y = state, times = times, func = lotka_volterra, parms = parameters)

# Convert output to data frame
df <- as.data.frame(out)
df;
pdf("d:/pdf/prey,pdf")
# Plot the results
plot(df$time, df$X, type = "l", col = "blue", xlab = "Time", ylab = "Population",
     main = "Lotka-Volterra Predator-Prey Model")
lines(df$time, df$Y, col = "red")
legend("topright", legend = c("Prey (X)", "Predator (Y)"), col = c("blue", "red"), lty = 1)

# Phase plot
plot(df$X, df$Y, type = "l", xlab = "Prey (X)", ylab = "Predator (Y)",
     main = "Phase Plot: Predator vs Prey")
fn_tosas9x(
      inp    = df
     ,outlib ="d:/sd1/"
     ,outdsn ="xtyt"
     );

;;;;
%utl_rendx;
libname sd1 "d:/sd1";

proc print data=sd1.xtyt;
run;quit;

/*----                                                                   ----*/
/*----  You may need to set the increment to .01 for x in r script above ----*/
/*----                                                                   ----*/

options ls=64 ps=44;
proc plot data=xtyt /*(where=(not missing(ltr)))*/;
 plot y*x='#' $ ltr/box  href=13.3 vref=15  haxis=0 to 40 by 5 vaxis=0 to 40 by 5;
run;quit;

options validvarname=upcase;
data xtyt /*(keep=x y rownames )*/;
  length ltr  $32;
  a = 1.5;
  b = 0.1;
  c = 0.075;
  d = 1.0;
  y_ab = a/b;
  x_dc = d/c;
  x0   = 10;
  y0   = 5;
  coexist=cats(put(x_dc,5.1),',',put(y_ab,5.1));
  set sd1.xtyt end=dne;;
  if rownames in  (29,158,243,389) then do;
   dxdt=a * X - b * X * Y;
   dydt=c * X * Y - d * Y;
   direction = atan((dydt)/(dxdt)) * 180/constant('pi');
   ltr1 = cats(put(dxdt,5.1),',',put(dydt,5.1),',',put(round(direction),4.)     );
   ltr = cats(put(x,5.1),',',put(y,5.1));
   ltr = cats('slope=',put(direction,4.));
   output;
 end;
 else output;
 drop ltr1;
run;quit;

proc print data=xtyt(where=(not missing(ltr))) ;
format _numeric_ 5.2;
run;quit;

options ls=64 ps=44;
proc plot data=xtyt /*(where=(not missing(ltr)))*/;
 plot y*x='*' $ ltr/box  href=13.3 vref=15  haxis=0 to 40 by 5 vaxis=0 to 40 by 5;
run;quit;

/**************************************************************************************************************************/
/*                                                                                                                        */
/*                                                                                                                        */
/*                                                                                                                                                                                                              */
/*        TIME     X     Y  ROW     LTR    A   B    C   D A/B D/CX X0 Y0  DXDT  DYDT DIRECTION                                                                                                                  */
/*                                                                                                                                                                                                              */
/*        0.28 13.28  4.82   29 slope=-0  1.5 0.1 0.075 1 15 13.33 10  5  3.52 -0.02   -0.09                                                                                                                    */
/*        1.57 35.55 14.85  158 slope=89  1.5 0.1 0.075 1 15 13.33 10  5  0.53 24.75   88.78                                                                                                                    */
/*        2.42 13.33 34.23  243 slope=0   1.5 0.1 0.075 1 15 13.33 10  5  25.6 -0.02    0.04                                                                                                                    */
/*        3.88  3.13 14.50  389 slope=-89 1.5 0.1 0.075 1 15 13.33 10  5  0.16 -11.1   -89.2                                                                                                                    */
/*                                                                                                                        */
/*               Plot of Y*X$LTR.  Symbol used is '#'.                                                                    */
/*                                                                                                                        */
/*     --+------+------+------+------+------+------+------+------+--                                                      */
/*  40 +                    |                                      +                                                      */
/*     |                    |                                      |                                                      */
/*     |                    |                                      |                                                      */
/*     |                    |                                      |                                                      */
/*  35 +                 slope= 0                                  +                                                      */
/*     |              **************                               |                                                      */
/*     |            ***     |      *****                           |                                                      */
/*     |          ***       |          ****                        |                                                      */
/*  30 +         **         |             ***                      +                                                      */
/*     |        **          |                **                    |                                                      */
/*   Y |        *           |                  ***                 |                                                      */
/*     |       **           |                    ***               |                                                      */
/*  25 +       *            |                      **              +                                                      */
/*     |      **            |                        **            |                                                      */
/*     |      *             |                         **           |                                                      */
/*     |      *             |                          **          |                                                      */
/*  20 +      *             |                           **         +                                                      */
/*     |      *             |                            **        |                                                      */
/*     |     **             |                             *        |                                                      */
/*     |     *              |                              *       |                                                      */
/*  15 +-----*-slope=--89---+--------------------slope=-89-*-------+                                                      */
/*     |     *              |                             **       |                                                      */
/*     |     **             |                             *        |                                                      */
/*     |      *             |                           ***        |                                                      */
/*  10 +      *             |                         ***          +                                                      */
/*     |      **            |                       ***            |                                                      */
/*     |       **           |                  *****               |                                                      */
/*     |        ****    slope= -0       ********                   |                                                      */
/*   5 +           **********************                          +                                                      */
/*     |                    |                                      |                                                      */
/*     |                    |                                      |                                                      */
/*     |                    |                                      |                                                      */
/*   0 +                    |                                      +                                                      */
/*     --+------+------+------+------+------+------+------+------+--                                                      */
/*       0      5     10     15     20     25     30     35     40                                                        */
/*                                                                                                                        */
/*------------------------------------------------------------------------------------------------------------------------*/
/*                                                                                                                        */
/*                                                                                                                        */
/*               --------------------------------------------------------------------------------------------             */
/*               |                     Dynamic Relation of Rat and Owl Populat|ions                         |             */
/*               |                     Preditor Prey Model (no closed form solution)                        |             */
/*               |                                                                                          |             */
/*               |        dx/dt = ax   - bxy                                     Parameters                 |             */
/*               |        dy/dt = cxy  - dy                                       a = 1.5                   |             */
/*               |                                                                b = 0.1                   |             */
/*               |        dx/dt=0 and dy/dt=0 when (x,y) = (d/c,a/b) nullclines   c = 0.075                 |             */
/*               |        when y=a/b => dx/dt=ax-ax=>y(t) =ax - ax = 0            d = 1.0                   |             */
/*               |                                                                                          |             */
/*               |        dx/dt =   1.5x  - 0.1xy                                                           |             */
/*               |        dy/dt = 0.075xy - 1.0y                                 Initial Conditions         |             */
/*               |                                                                x0   = 10                 |             */
/*               |        Definition of Parameters                                y0   = 5                  |             */
/*               |        ------------------------                                                          |             */
/*               |                                                               NullClines                 |             */
/*               |        dRat/dt = aRat   - bRatOwl                              y_ab = a/b                |             */
/*               |        dOwl/dt = cOwlRat- cOwl                                 x_dc = d/c                |             */
/*               |                                                                                          |             */
/*               |        a - rat birth rate                                PROPERIES                       |             */
/*               |                                                          ---------                       |             */
/*               |        b - rate at which owls                            Tangents at the intersection    |             */
/*               |            encounter and consume rats                    of nullclines are               |             */
/*               |            affects both populations                      perpendicular to nullclines     |             */
/*               |            (coupling term)                                                               |             */
/*               |                                                          Signs change across nullclines  |             */
/*               |        c - owls capacity to reproduce adjusted                                           |             */
/*               |            proportionally to the product of owl          rat,owl nullcliness             |             */
/*               |            and rat encounters                            steady state for owls & rats    |             */
/*               |                                                          at y=a/b dx/dt=0 a/b nullcline  |             */
/*               |        d - natural death rate of owls in                    x=d/c dy/dt=0 d/c nullcline  |             */
/*               |            the absense of prey                                                           |             */
/*               |            -c*owl (no rat food worst case?)              (0,0) is also a solution                      */
/*               |                                                          to dx/dt and dy/dt extinction                 */
/*               |                                                                                 -at                    */
/*               |                                                          No rats (x=0) owl(t)=Ae                       */
/*               |                                                          Owl extinction      at                        */
/*               |                                                          No owls rats(t) = Ae                          */
/*               |                                                          Expontial growth                              */
/*               |                                        PREY                                                            */
/*               |---------------------+------+------+------+------+------+------+------+------+------------|             */
/*               |                     0      5     10     15     20     25     30     35     40            |             */
/*               |                                                                                          |             */
/*               |        Calculation of direction is      |                                                |             */
/*               |        based on the vector direction    |                                                |             */
/*               |                                                                                          |             */
/*               |        degrees=180*atan(dydt/dxdt)/P   a/b=15  d(rats)/dt =0 (steady state               |             */
/*            40 +                                         |                     rates do not increase      + 40          */
/*               |                                                               or decrease over time)     |             */
/*               |                                  <--------------                                         |             */
/*               |                                                                                          |             */
/*            35 +                                    dx/dt  dy/dt   o (degrees)                            + 35          */
/*               |                                    -25.6, -0.0  ,0                                       |             */
/*               |                                   -----[+]-----o                                         |             */
/*               |                                 ###   X   Y    #####                                     |             */
/*            30 +                               ###   13.3,34.2     ####                                   + 30          */
/*               |            Positive slope -->##         |      o      ###                                |             */
/*               | Owls decrease Rats Increase ##      <-------- 0          ##                              |             */
/*               |                            ##           |                  ### <- Negative slpoe         |             */
/*            25 +                           *#            |                    ###    Rats decrease        + 25          */
/*     OWLS      |                          ##                                    ##   Owls increase (lag)  |             */
/*     PREDATOR  |                          ##            d/c                       ##                      | R           */
/*               | d/c=13.3 (steady state  ##                                        ##                     |             */
/*            20 + owls do not increase    ##           o  |                          ##                    + 20          */
/*               | or decrease)            #         -89   |          ^   o            ##      ^            |             */
/*               |         |               #           |   |          | 89              ##     |            |             */
/*               |         |   Initial     #           |   |          |                  #     |            |             */
/*      d     15 +         |  Population   # dxdt dyst |   0          |      dxdt dydt   o#    |            + 15          */
/*      - = 13.3 |-d/c=13--|------------   # -0.2,-11.1|-89+----------|-----  0.5,24.7,89 #----|--- d/c ----|             */
/*      c        |         |               ##          |   | 15,13.3  |        X    Y    ##    |            |             */
/*               |         V                #          V   | coexist  |      35.5,14.9 [+}                  |             */
/*            10 +          Negative Slope  ##             |      o                   ###                   + 10          */
/*               |Rats decrease Owls Increase**        ---------> 0                  ###                    |             */
/*               |                          >##            |                       ###<--Positve slope      |             */
/*               |           Initial           #    dxdt dydt  o            #####                           |             */
/*             5 +           Population [10,5] #[+]   13.5,-0.0,0      ########                             +  5          */
/*               |                                -------[+]-----------                                     |             */
/*               |                                       X  Y                                               |             */
/*               |                                     13.3,4.8                                             |             */
/*             0 +                                         |                                                +  0          */
/*               |                                   ------------->                                         |             */
/*               |                                       a/b=15                                            |              */
/*               |------------------------------------------------------------------------------------------|             */
/*               |                                                                                          |             */
/*               |                               INITIAL VALUES                                             |             */
/*               |                                                                                          |             */
/*               |    X       Y      DIRECTION     X0     Y0      COEXIST     DXDT     DYDT  ROWS           |             */
/*               |                                                                                          |             */
/*               |   13.3     4.8       -0.1      10.0    5.0    15,13.3      13.5     -0.0    29           |             */
/*               |   35.5    14.9       88.8      10.0    5.0    15,13.3       0.5     24.7   158           |             */
/*               |   13.3    34.2        0.0      10.0    5.0    15,13.3     -25.6     -0.0   243           |             */
/*               |    3.1    14.5      -89.2      10.0    5.0    15,13.3       0.2    -11.1   389           |             */
/*               |                                                                                          |             */
/*               |                                         |                                                |             */
/*               |                                         |                                                |             */
/*               ---------------------+------+------+------+------+------+------+------+------+--------------             */
/*                                    0      5     10     15     20     25     30     35     40                           */
/*                                                        a/b                                                             */
/*                                                       PREY                                                             */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*  _          __  __           _     _       _ _   _       _                       _         _
| || |    ___ / _|/ _| ___  ___| |_  (_)_ __ (_) |_(_) __ _| |   ___ ___  _ __   __| |  _ __ | |__   __ _ ___  ___
| || |_  / _ \ |_| |_ / _ \/ __| __| | | `_ \| | __| |/ _` | |  / __/ _ \| `_ \ / _` | | `_ \| `_ \ / _` / __|/ _ \
|__   _||  __/  _|  _|  __/ (__| |_  | | | | | | |_| | (_| | | | (_| (_) | | | | (_| | | |_) | | | | (_| \__ \  __/
   |_|   \___|_| |_|  \___|\___|\__| |_|_| |_|_|\__|_|\__,_|_|  \___\___/|_| |_|\__,_| | .__/|_| |_|\__,_|___/\___|
                                                                                       |_|
*/      */

%utl_rbeginx;
parmcards4;
library(deSolve)
library(haven)
source("c:/oto/fn_tosas9x.R");
# Load required library

# Define the Lotka-Volterra model
lotka_volterra <- function(t, state, parameters) {
  with(as.list(c(state, parameters)), {
    dP <- a * P - b * P * V
    dV <- -c * V + d * P * V
    list(c(dP, dV))
  })
}

# Set parameters
parameters <- c(a = 1.5, b = 0.1, c = 1, d = 0.075)

# Set initial state

# Set time steps
times <- seq(0, 40, by = 0.01)

pdf("d:/pdf/prey.pdf")
state <- c(P = 10, V = 5)
# Solve the differential equations
out <- ode(y = state, times = times, func = lotka_volterra, parms = parameters)
state <- c(P = 5, V =2.5)
# Solve the differential equations
out1 <- ode(y = state, times = times, func = lotka_volterra, parms = parameters)
state <- c(P =2.5, V =1.25 )
# Solve the differential equations
out2 <- ode(y = state, times = times, func = lotka_volterra, parms = parameters)
# Create the phase space plot

plot(out[, "P"], out[, "V"], type = "l",
     xlab = "Prey Population", ylab = "Predator Population",
     main = "Phase Space Plot: Lotka-Volterra Model")
points(state["P"], state["V"], col = "red", pch = 16)
fn_tosas9x(
      inp    = out
     ,outlib ="d:/sd1/"
     ,outdsn ="out"
     );
fn_tosas9x(
      inp    = out1
     ,outlib ="d:/sd1/"
     ,outdsn ="out1"
     );
fn_tosas9x(
      inp    = out2
     ,outlib ="d:/sd1/"
     ,outdsn ="out2"
     );
pdf()
;;;;
%utl_rendx;

proc sql;
  create
    table outs as
  select
    l.p as owl1
   ,l.v as rat1
   ,c.p as owl2
   ,c.v as rat2
   ,r.p as owl3
   ,r.v as rat3
 from
   sd1.out            as l
   left join sd1.out1 as c on l.time = c.time
   left join sd1.out2 as r on l.time = r.time
;quit;


options ls=100 ps=64;
proc plot data=outs;
  plot owl1 *rat1='1'  owl2 *rat2='2' owl3*rat3='3'/overlay box;
run;quit;


/**************************************************************************************************************************/
/*                                                                                                                        */
/*                                                  RATS                                                                  */
/*        0           10           20           30           40           50           60           70                    */
/*        -+------------+------------+------------+------------+------------+------------+------------+--                 */
/*        |                                                                                             |                 */
/*        |                                                                                             |                 */
/*     80 +               33333333333                             This shows birth and death cycles     + 80              */
/*        |            333           33333                        we saw inte the parametric time       |                 */
/*        |          333                 3333                     series plots                          |                 */
/*        |         33                       333                                                        |                 */
/*        |        33                           333               The larger the loops the larger the   |                 */
/*        |       33                              3333            extrema.                              |                 */
/*     70 +      33                                  333                                                + 70              */
/*        |      3      Initial Conditions             333                                              |                 */
/*        |     33       x0   = 2,5                       333                                           |                 */
/*        |     3        y0   = 1.25                        333                                         |                 */
/*        |    33                                             333                                       |                 */
/*        |    3                                                33                                      |                 */
/*     60 +    3                                                  33                                    + 60              */
/*        |   33                                                    33                                  |                 */
/*        |   3           2222222222                                  33                                |                 */
/*        |   3        2222        222222                               33                              |                 */
/*        |   3       22                2222                              33                            |                 */
/*        |   3     222                    2222                            33                           |                 */
/*     50 +  33    22                         2222                           33                         + 50              */
/*        |  3     2                             222                           33                       |                 */
/*        |  3    22   Initial Conditions           222                          33                     |                 */
/*        |  3   22     x0   = 5                     222                          33                    |                 */
/*        |  3   2      y0   = 2.5                    222                          33                   |                 */
/*        |  3  22                                       222                         33                 |                 */
/* OWL 40 +  3  2                                          222                         33               + 40  OWLS        */
/*        |  3  2                                            222                        33              |                 */
/*        | 33 22                                              22                         33            |                 */
/*        | 3  2         1111111111111                          222                        33           |                 */
/*        | 3  2       111           1111                         22                        33          |                 */
/*        | 3  2      11                1111                       222                       33         |                 */
/*     30 + 3 22     11                    111                       22                        33       + 30              */
/*        | 3 2     11                       111                      22                        33      |                 */
/*        | 3 2    11                          111                      22                       33     |                 */
/*        | 3 2    1                             11                      22                       33    |                 */
/*        | 3 2   11                              111                     22                       33   |                 */
/*        | 3 2   1                                 11                     22                       3   |                 */
/*     20 + 3 2   1    Initial Conditions            11                     2                        3  + 20              */
/*        | 3 2   1     x0   = 10                     11                    22                       33 |                 */
/*        | 322  11     y0   = 5                       1                     22                       3 |                 */
/*        | 32   1                                     1                      2                       3 |                 */
/*        | 32   1                                     11                     2                       3 |                 */
/*        | 32   1                                     1                     22                       3 |                 */
/*     10 + 3 2  11                                   11                     2                       33 + 10              */
/*        | 3 2   1                                 111                    222                      33  |                 */
/*        | 3 2   11                             1111                    222                      333   |                 */
/*        | 3 22   1111                    1111111                   22222                     3333     |                 */
/*        | 33 22     1111111111111111111111                 222222222                  3333333         |                 */
/*        |  33 2222222222222222222222222222222222222222222222        3333333333333333333               |                 */
/*      0 +   333333333333333333333333333333333333333333333333333333333                                 +  0              */
/*        |                                                                                             |                 */
/*        -+------------+------------+------------+------------+------------+------------+------------+--                 */
/*         0           10           20           30           40           50           60           70                   */
/*                                                  RATS                                                                  */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*___     __ _ _         _                   _                    _
| ___|   / _(_) |_   ___(_)_ __   ___  ___  | |_ ___    _ __ __ _| |_ ___
|___ \  | |_| | __| / __| | `_ \ / _ \/ __| | __/ _ \  | `__/ _` | __/ __|
 ___) | |  _| | |_  \__ \ | | | |  __/\__ \ | || (_) | | | | (_| | |_\__ \
|____/  |_| |_|\__| |___/_|_| |_|\___||___/  \__\___/  |_|  \__,_|\__|___/

*/

data addTerm;

   pi=constant('pi');

   set sd1.xtyt(where= (3.84 < time <9.49)); /* first complete cycle of rats */

   time=(time-3.84)*3.14/5.65;

   sin1 = sin(1*time);
   sin2 = sin(2*time);
   sin3 = sin(3*time);
   sin4 = sin(4*time);
   sin6 = sin(6*time);
   sin8 = sin(8*time);

   keep x time sin: ;

run;quit;

proc reg data=addTerm;
  model x = sin1 sin2 sin3 sin4  sin6 ;
  output out=res(keep=time x p resid) pred=p r=resid
run;quit;

data back;
  set res;
  time = 5.65/3.14*time + 3.84;
run;quit;

options ls=110 ps=64;
proc plot data=back;
 plot resid*time='*' /box;
 plot  p*time='.' x*time='o' / overlay box;
 plot  time*p='.' time*x='o' / overlay box;
run;quit;


/**************************************************************************************************************************/
/*                                                                                                                        */
/*                                                                                                                        */
/*       -+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+-          */
/*        |                                                                                                     |         */
/*        |                                                                                                     |         */
/*    2.0 +    RESIDUALS First complete cycle of RATS(t) X(t)                                                   +         */
/*        |                                                                                                     |         */
/*        |    RESIDUALS vs Time                                                                                |         */
/*        |                                                                                                     |         */
/*        |                                                                                                     |         */
/*        |                                                              ***                                    |         */
/*        |                                                              * *                                    |         */
/*    1.5 +                                                             ** *                                    +         */
/*        |                                                             *   *                                   |         */
/*        |                                                             *   *                                   |         */
/*        |                                                             *   *                                   |         */
/*        |                                                            *    *                                   |         */
/*        |                                                            *     *                                  |         */
/*        |                                      ***                   *     *                                  |         */
/*    1.0 +                                     ** **                  *     *                                  +         */
/*        |                                     *   *                  *     *                                  |         */
/*        |                                    **    *                *      *                                  |         */
/*        |            ***                     *     *                *       *                                 |         */
/*        |           ** **                   **     **               *       *                                 |         */
/*        |          **   *                   *       *               *       *               *****             |         */
/* R      |          *     *                  *       *               *       *               *   **            |         */
/* e  0.5 +         **     *                 **       **              *       *              **    *            +         */
/* s      |         *       *                *         *             *        *              *     **           |         */
/* i      |         *       *                *         *             *         *            **      *           |         */
/* d      |        *        **              **         *             *         *            *        *          |         */
/* u      |        *         *              *          **            *                      *        *          |         */
/* a      |        *         *              *           *            *         *           *          *         |         */
/* l      |       **         **            **           *            *         *           *          *         |         */
/*    0.0 +       *           *            *            *           *          *           *           *        +         */
/*        |       *           *            *            **          *          *          *            *        |         */
/*        |       *            *           *             *          *           *         *            **       |         */
/*        |      *             *          *              *          *           *         *             *       |         */
/*        |      *             *          *              *          *           *         *             **      |         */
/*        |      *              *        **              **         *           *        *               *      |         */
/*        |      *              *        *                *        *            *        *               **     |         */
/*   -0.5 +     *               **       *                *        *            *        *                *     +         */
/*        |     *                *      **                *        *             *      **                 *    |         */
/*        |     *                **     *                 **       *             *      *                  **   |         */
/*        |     *                 *    **                  *       *             *      *                   *   |         */
/*        |    *                  **   *                   *      *              *      *                   **  |         */
/*        |                        ****                    *      *              **    *                     *  |         */
/*        |                          *                     **     *               *    *                        |         */
/*   -1.0 +                                                 *     *               *    *                        +         */
/*        |                                                 *    *                *   *                         |         */
/*        |                                                 **   *                 *  *                         |         */
/*        |                                                  *  **                 * **                         |         */
/*        |                                                  ** *                   **                          |         */
/*        |                                                   ***                                               |         */
/*        |                                                                                                     |         */
/*   -1.5 +                                                                                                     +         */
/*        |                                                                                                     |         */
/*        -+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+-         */
/*        3.6       4.2       4.8       5.4       6.0       6.6       7.2       7.8       8.4       9.0       9.6         */
/*                                                                                                                        */
/*                                                         TIME                                                           */
/*                                                                                                                        */
/*      --+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+--         */
/*      |                                                                                                       |         */
/*      |                                                                                                       |         */
/*   40 +     Predicted vs Observed  first rats(t) period(cycle)                                                +         */
/*      |                                                                                                       |         */
/*      |     '.' ~Predicted & observed on top od each other                                                    |         */
/*      |     'o' ~Observed diff from pred ooo                                                                  |         */
/*      |                                                                                                       |         */
/*      |                                                            ooo                                        |         */
/*   35 +                                                           .. oo                                       +         */
/*      |                                                         ... .. o                                      |         */
/*      |                                                        ..     ..o                                     |         */
/*      |                                                       ..o      .oo                                    |         */
/*      |                                                      ..o        .o                                    |         */
/*      |                                                      .o         ..o                                   |         */
/*   30 +                                                     .o           .o                                   +         */
/*      |                                                    ..o            .                                   |         */
/* P    |                                                    .o             .o                                  |         */
/* r    |                                                   ..o              .                                  |         */
/* e    |                                                   .o               ..                                 |         */
/* d    |                                                  .o                 .                                 |         */
/* i 25 +                                                 ..o                 .                                 +         */
/* c    |                                                 .o                   .                                |         */
/* t    |                                                ..                    .                                |         */
/* e    |                                                .o                    o.                               |         */
/* d    |                                               ..                      .                               |         */
/*      |                                               .                       o.                              |         */
/* V 20 +                                              ..                        .                              +         */
/* a    |                                             o.                         ..                             |         */
/* l    |                                             .                          o.                             |         */
/* u    |                                            ..                           ..                            |         */
/* e    |                                           o.    R-Square     0.9949     o.                            |         */
/*      |                                          o..                             ..                           |         */
/* o 15 +                                         o..                              o.                           +         */
/* f    |                                        o..                                ..                          |         */
/*      |                                       o..                                 o.                          |         */
/* X    |                                      o..                                   ..                         |         */
/*      |                                     o..                                    o.                         |         */
/*      |                                    o..                                      ..                        |         */
/*   10 +                                  o...                                       o..                       +         */
/*      |                                ....                                          o.                       |         */
/*      |                             ....o                                             ..                      |         */
/*      |                           ...oo                                                ..                     |         */
/*      |                        ...ooo                                                   ..                    |         */
/*      |                      ...ooo                                                      ..o                  |         */
/*    5 +                    ...oo                                                          ..oo                +         */
/*      |     ..        ooo...                                                                ...ooo      ...   |         */
/*      |     o.....oooo....                                                                    ...........oo   |         */
/*      |          ......                                                                                       |         */
/*      |                                                                                                       |         */
/*      |                                                                                                       |         */
/*    0 +                                                                                                       +         */
/*      |                                                                                                       |         */
/*      --+---------+---------+---------+---------+---------+---------+---------+---------+---------+---------+--         */
/*       3.6       4.2       4.8       5.4       6.0       6.6       7.2       7.8       8.4       9.0       9.6          */
/*                                                                                                                        */
/*                                                        TIME                                                            */
/*                                                                                                                        */
/*       ---+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+---          */
/*  TIME |                                                                                                     |          */
/*   9.6 +                   Observed vs Predicted  first rats(t) period(cycle)                                +          */
/*       |         oo.                                                                                         |          */
/*       |          o.       '.' ~Predicted & observed on top of each other                                    |          */
/*       |          ..       'o' ~Observed diff from pred ooo                                                  |          */
/*       |          .                                                                                          |          */
/*   9.0 +          .o        Better visualization of residuals                                                +          */
/*       |          .oo                                                                                        |          */
/*       |           ..o                                                                                       |          */
/*       |            ...o                                                                                     |          */
/*       |              ....                                                                                   |          */
/*   8.4 +                 .....                                                                               +          */
/*       |                     .......                                                                         |          */
/*       |                          o........                                                                  |          */
/*       |                                oo.........                                                          |          */
/*       |                                       ooo.........                                                  |          */
/*   7.8 +                                                 oo.........                                         +          */
/*       |                                                           .........o                                |          */
/*       |                                                                    .......ooo                       |          */
/*       |                                                                           ......oooo                |          */
/*       |                                                                                .....ooo             |          */
/*   7.2 +                                                                                    .. o             +          */
/*       |                                                                                    ..oo             |          */
/*       |                                                                               ooo...                |          */
/*       |                                                                          ooo.....                   |          */
/*       |                                                                   oooo......                        |          */
/*   6.6 +                                                             ooo.......                              +          */
/*       |                                                       o........                                     |          */
/*       |                                                 ........                                            |          */
/*       |                                          .......o                                                   |          */
/*       |                                     ......o                                                         |          */
/*   6.0 +                                .....ooo                                                             +          */
/*       |                             ....ooo                                                                 |          */
/*       |                          ....oo                                                                     |          */
/*       |                        ...oo                                                                        |          */
/*       |                      o..o                                                                           |          */
/*   5.4 +                    oo..                                                                             +          */
/*       |                  oo..                                                                               |          */
/*       |                 oo..                                                                                |          */
/*       |               oo...                                                                                 |          */
/*       |              o...                                                                                   |          */
/*   4.8 +             ...                                                                                     +          */
/*       |           ...                                                                                       |          */
/*       |          ..o                                                                                        |          */
/*       |         ..o                                                                                         |          */
/*       |         .oo                                                                                         |          */
/*   4.2 +         .o                                                                                          +          */
/*       |         ..                                                                                          |          */
/*       |          ..                                                                                         |          */
/*       |         oo.                                                                                         |          */
/*       |                                                                                                     |          */
/*   3.6 +                                                                                                     +          */
/*       ---+-----------+-----------+-----------+-----------+-----------+-----------+-----------+-----------+---          */
/*          0           5          10          15          20          25          30          35          40             */
/*                                                                                                                        */
/*                                                Predicted Value of X                                                    */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*__                                                       _ _ _ _          _
 / /_    ___ _   _ _ __ ___  _ __  _   _   ___  __ _ _   _(_) (_) |__  _ __(_)_   _ _ __ ___
| `_ \  / __| | | | `_ ` _ \| `_ \| | | | / _ \/ _` | | | | | | | `_ \| `__| | | | | `_ ` _ \
| (_) | \__ \ |_| | | | | | | |_) | |_| ||  __/ (_| | |_| | | | | |_) | |  | | |_| | | | | | |
 \___/  |___/\__, |_| |_| |_| .__/ \__, | \___|\__, |\__,_|_|_|_|_.__/|_|  |_|\__,_|_| |_| |_|
             |___/          |_|    |___/          |_|
*/

Coexistance equlilibrium

   0  = a*rat  - b*rat*owl
   0  = -c*owl + d*rat*owl



%utl_pybeginx;
parmcards4;
import pyperclip;
from sympy import solve, symbols;
rat,owl,a,b,c,d = symbols("rat owl a b c d");
res=solve([a*rat  - b*rat*owl
         ,-c*owl + d*rat*owl] , rat, owl, dict=True);
print(res);
pyperclip.copy(str(res));
;;;;
%utl_pyendx(return=pts);

/**************************************************************************************************************************/
/*                                                                                                                        */
/* %put &=pts;                                                                                                            */
/*                                                                                                                        */
/* PTS=[{owl: 0, rat: 0}, {owl: a/b, rat: c/d}]                                                                           */
/*                                                                                                                        */
/* Check by substitution                                                                                                  */
/*                                                                                                                        */
/* Does                                                                                                                   */
/*   a*rat  - b*rat*owl      = 0                                                                                          */
/* Yes                                                                                                                    */
/*   a*(c/d) -b*(c/d)*(a/b)  = 0                                                                                          */
/* Likewise for                                                                                                           */
/* Does                                                                                                                   */
/*   -c*owl + d*rat*owl      = 0                                                                                          */
/* Yes                                                                                                                    */
/*   -c*(a/b) + d*(c/d)*(a/b = 0                                                                                          */
/*                                                                                                                        */
/**************************************************************************************************************************/

/*____                                 _               _   _
|___  |  _ __ ___ _ __ ___   _____   _(_)_ __   __ _  | |_(_)_ __ ___   ___
   / /  | `__/ _ \ `_ ` _ \ / _ \ \ / / | `_ \ / _` | | __| | `_ ` _ \ / _ \
  / /   | | |  __/ | | | | | (_) \ V /| | | | | (_| | | |_| | | | | | |  __/
 /_/    |_|  \___|_| |_| |_|\___/ \_/ |_|_| |_|\__, |  \__|_|_| |_| |_|\___|
                                               |___/
*/


 dx/dt = ax   - bxy
 dy/dt = cxy  - dy

 Lets divide

 dy/dt =   cxy  - dy      dy   y ( cx - d )
 ---------------------- = -- = -------------
 dx/dt =   ax   - bxy     dx   x ( a - by )


 x ( a - by ) dy  =  y ( cx - d ) dx

 x ( a - by ) dy  -  y ( cx - d ) dx



 Divide by xy

 x ( a - by ) dy  =  y ( cx - d ) dx
 ---------------     ----------------
      xy                  xy

 ( a - by )            ( cx - d )
 ---------- dy    =   ----------- dx
      y                    x

 ( a/y - b ) dy  =  ( c - d/x) dx


 a ( 1/y - b/a) dy  = d ( c/d - 1/x) dx

 Note the integral of 1/x = log(x)


  /
  | a ( 1/y - b/a) dy  = -by + alog(y) + c1
 /

  /
  | d ( c/d - 1/y) dx  = cx - d*log(x)  + c2
 /


  c*x + b*y - a*log(y) - d*log(x)



/*              _
  ___ _ __   __| |
 / _ \ `_ \ / _` |
|  __/ | | | (_| |
 \___|_| |_|\__,_|

*/
















































































































































































