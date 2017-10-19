#####
Input
#####

All input files must contain an ``[input]`` section looking like this:

.. code::

    [input]
    version = 1

Introducing a ``version`` key helps us to make changes to the input file
format while keeping compatibility with previous formats. Please note
that Lumol is not in version 1.0 yet and we currently cannot guarantee
compatibility for input files.

The input files can also contain a ``[log]`` section to control where
should the code output be printed. Please see the corresponding
`documentation <input/log.html>`__ for more information.

.. toctree::
   :hidden:

   log
