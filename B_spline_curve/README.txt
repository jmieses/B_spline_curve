This program creates a B-spline curve using the deBoor-Cox recursion formula. It implements a set of function to create the curve. 
In addition, the Bernstein polynomial that create a Bezier curve works in this program too. 

The actual B-spline curve presented in this program is not reaching the last control point as a clamped curve. The implementation 
of the curve algorithm included a break constraint for index vector out of range. Thus, one could speculate that the control
point is not being reached due to this code fixed. 

In the other hand, the curve works for an arbitrary set of control points, and for an arbitrary degree curve. 