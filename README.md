# WIP: TernaryLBM (workflow discussion)

## WP1: Modularisation, classes

What modules will be part of the software? Software requirements?

Note: Introduce goals and DLs
Note: Milestones (MS): MS may belong or may not belong to WP
1) Module Evaporation
  - Issue ## will be connected with a ML; when everything has been implemented and issues have been resolved --> ML is achieved
3) Module Condensation

Breaking into steps.
1) (One LP) Calculating on single lattice point. 2 weeks.
2) (Many LPs) Implementing different datatypes, the way you access data. Load balancing, how data accessed
2a) Sequential (one process without communication)
2b) Parallel
Note for (2) - inbetween: some test cases, benchmarking, in the next 2-3 weeks --> then we dıscuss and do how we benchmark ıt
3) Performance analysis and optimisation of SW as a whole.
4) While adding new modules, documentation, correctness, results.

## WP2: Performance analysis and optimisation

How do we measure (approaches, tools)? Where are the bottlenecks? How do we deal will them

## WP3: Results, correctness, demonstration

Which results are we interested? How do we check correctness? What will be demonstrated?
