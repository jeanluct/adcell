* Include `lecture28.tex` writeup as a doc?

* Increase max before complaining about lack of MEX file?  It's
  hardwired in `fft2dotgrad.m`, so I'm loath to change it.

* Add `isempty` option for default arguments.

* Make `private` folder?  But I think right now all the functions are
  used.  Function to return k vectors?  (Add to jlt.)

* `examples` folder?  Rename `tarek.m` to something more
  illustrative.  Maybe `onecell.m`?

* Help messages!

* Write `adcell.effdiff` function.  Returns Deff and optional chi.
  Compare to large-D limiting form.

* To get ready for eventual distribution:

  * Add LICENSE and COPYING.  GPL v3 ok?  MIT to agree with *jlt lib*?

  * Add boilerplate to m-files.

  * Don't forget to credit `fourdif` authors.  No license?
