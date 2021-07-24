This project contains implemetations of algorithms used to solve huge linear systems called Krylov methods.
In particular one can find here implementations of algorithms both with and without preconditioning:

- 'GMRES/' implements a GMRES(m) method with restarts 
- 'GCR/' implements a GCR(m) method with restarts
- 'ORTHODIR/' implements a ORTHODIR(m) method with restarts
- 'ORTHOMIN/' implements a ORTHOMIN(m,k) method with restarts

For improvement of algorithms, the conditioners are used:

- 'Conditioners/' implements some popular conditioners for improvement of methods

We are testing the methods and making plots which are later saved in the 'PLOTS/' directory in the TEST folders:

- 'TESTS_no_preconditioned/' tests methods without preconditioning
- 'TESTS_preconditioned/' tests methods with preconditioning
- 'TESTS_both/' tests both methods simultaniousely


