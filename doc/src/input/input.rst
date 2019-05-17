``[input]`` section
*******************

All input files must contain an ``[input]`` section looking like this:

.. code::

    [input]
    version = 1

Introducing a ``version`` key helps us to make changes to the input file format
while keeping compatibility with previous formats. Please note that Lumol is not
in version 1.0 yet and we currently cannot guarantee compatibility for input
files.
