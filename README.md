# reentry
A simple simulation of Apollo capsule reentry. 

reentry-ld-controlled.py: L/D (lift to drag ratio) is PID controlled to achieve the nominal range to splashdown. 

reentry-ballistic.py: no bank angle controll, i.e. fix L/D

Fun fact: If L/D is too high, the capsule skips back to space for a more brutal second g-peak.

rk45.py: Runge-Kutta ODE solvers of orders 4 and 5

usatmos.py: U.S. Standard Atmosphere, 1976

