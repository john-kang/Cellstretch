# Cellstretch

Code for publications in:

1. Response of an actin filament network model under cyclic stretching through a coarse grained Monte Carlo approach.
Kang J, Steward RL, Kim Y, Schwartz RS, LeDuc PR, Puskar KM.
J Theor Biol. 2011 Apr 7;274(1):109-19. doi: 10.1016/j.jtbi.2011.01.011. Epub 2011 Jan 15.
https://www.ncbi.nlm.nih.gov/pubmed/21241710
This paper describes Monte Carlo simulations of actin filaments under cyclic stretch conditions)

2. Structurally Governed Cell Mechanotransduction through Multiscale Modeling
John Kang, Kathleen M. Puskar, Allen J. Ehrlicher, Philip R. LeDuc & Russell S. Schwartz
Scientific Reports volume 5, Article number: 8622 (2015)
https://www.nature.com/articles/srep08622
This paper describes Monte Carlo simulations of actin filaments under both cyclic stretch and shear conditions to simulate mechanotransduction signaling

It's been a while since I wrote this code and the annotation/comments may be confusing, so please feel free to contact me regarding questions

The main function is "cellstretch_john_main.m" for input parameters
This main file calls "cellstretch_john_execute.m" which calls further helper functions
