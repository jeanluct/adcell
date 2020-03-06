* Include `lecture28.tex` writeup as a doc?

* Make `private` folder?  But I think right now all the functions are
  used.  Function to return k vectors?  (Add to jlt.)

* `examples` folder?  Rename `tarek.m` to something more
  illustrative.  Maybe `onecell.m`?

* Check sign of u.  The operator should be -u.grad?

* Laplacian function to return diffusive part only.  This is useful so
  we can make many Ak's for different diffusivities with recomputing
  u.grad.

* Help messages!

* Write `adcell.effdiff` function.  Returns Deff and optional chi.
  Compare to large-D limiting form.

* To get ready for eventual distribution:

  * Add LICENSE and COPYING.  GPL v3 ok?  MIT to agree with *jlt lib*?

  * Add boilerplate to m-files.

  * Don't forget to credit `fourdif` authors.  No license?
