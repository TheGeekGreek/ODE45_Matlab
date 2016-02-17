# Project: Differential Equations

Run `sample_direction_field_script.m` for a sample implementation. See `f.m` for some sample functions.

## Sample Output

For `f(t,x) = -0.5y`

![alt text](https://github.com/TheGeekGreek/ODE45_Matlab/blob/master/sample_direction_field_script_plot.png "vector field")


## Doc

To get detailed usage information: 

`help direction_field`;

```  
direction_field:  
        Draws a normed direction field using an ODE 
 
  parameters: ( function, timespan, y0, x)
    function -> f(t,y)
    timespan -> [start_time, end_time]
    y0       -> starting value
    x        -> lower_bound:step_size:upper_bound
 
  The vector field will be on a square grid on the interval of x, y0
  defines the starting value of the differential, timespan is the interval
  on which we will solve, and f is a function f(t,y)
```

`help explicit_euler`;

```
  explicit_euler:
        Solves a differential using the Explicit Euler
 
  parameters: ( function, timespan, y0, bound)
    function -> f(t,y)
    timespan -> [start_time, end_time]
    y0       -> starting value
    steps    -> the number of steps
```

`help ODE45`;

```
  ODE_45:
        Solves a differential using the Embbedded Runge-Kuta Method
 
  parameters: ( function, timespan, y0)
    function -> f(t,y)
    timespan -> [start_time, end_time]
    y0       -> starting value
```
