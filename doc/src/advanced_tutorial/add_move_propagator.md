# Adding the move to the propagator

This section will talk about:

- How to add a move to the propagator
- As an example we write a NVT simulation of an atomic liquid using only translation.
(Comment: not sure if this causes confusion since it adds a lot of functions we did not talk about)
- We will look at the performance using unix `time` functionality and the output of how much moves will be accepted. 
(Comment: time to motivate using cache, acceptance to motivate updating the displacement)