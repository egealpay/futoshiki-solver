# Futoshiki Solver

Futoshiki is a board-based puzzle game, also known under the name Unequal. It is playable on a square board having a given fixed size (4x4 for example). Techniques for solving could be found at futoshiki.org. In this repository, you can find 3 different implementation for solving a Futoshiki puzzle. 

## Sequential Implementation
This implementation solves a Futoshiki puzzle in sequential order with a single thread.

You can find the detials from this [report](https://github.com/egealpay/futoshiki-solver/blob/master/sequential/sequential-report.pdf) 


## Parallel Implementation
This implementation solves a Futoshiki puzzle in parallel order with multiple threads. OpenMP is used to manage threads.
Backtracking algorithm was used in 3 different designs. 

You can find the detials from this [report](https://github.com/egealpay/futoshiki-solver/blob/master/parallel/parallel-report.pdf) 


## GPU Implementation
This implementation solves a Futoshiki puzzle on a GPU with CUDA. 

You can find the detials from this [report](https://github.com/egealpay/futoshiki-solver/blob/master/gpu/gpu-report.pdf) 
