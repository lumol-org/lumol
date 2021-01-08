``[log]`` section
*****************

Lumol sends various logging messages while running a simulation. Some of them
are informational messages (``charge was set to 1.2 for 132 particles``), others
are warnings or error message (``infinite energy!``) and some are for debugging
purposes.

By default, Lumol prints all the informational, warning and error messages to
the standard terminal output. This allows to run the code and redirect the
output to a specific file in the usual UNIX way: ``lumol input.toml >
simulation.log``.

Lumol also offers more detailed configuration for logging output, for example if
you only want to print errors and warnings, and redirect everything else to a
file. This configuration happens in the ``[log]`` section of the input file.
This section can contain either a single output ``target``, or multiple
``targets``.

**Example**

.. code::

    # Single log target
    [log]
    target = "lumol.log"

    # Multiple log targets
    [log]
    targets = [
        {target = "lumol.log"},
        {target = "<stdout>", level = "warning"}
    ]

    # Multiple log targets, alternative syntax
    [[log.targets]]
    target = "lumol.log"

    [[log.targets]]
    target = "<stdout>"
    level = "warning"

For multiple targets, the ``targets`` key must be an array of tables (either
indicated by two brackets, i.e. ``[[log.targets]]`` or by defining the table
``targets = [{target = ...}, {target = ...}, ...]``), containing individual
target configuration with the same syntax as a single target.

The only required key is the ``target`` key, identifying where to write the
messages. It is interpreted as  path to a file, except for the two special cases
of ``<stdout>`` and ``<stderr>`` that are used to write messages to the standard
terminal output stream or error stream respectively.

.. code::

    # Write all messages to the standard error stream
    [log]
    target = "<stderr>"

Optional keys are ``level`` and ``append``. The ``level`` key controls which
messages are sent to the specified output and defaults to ``info``.  Available
levels are the following:

-  ``error``: only error messages;
-  ``warning``: error and warning messages;
-  ``info``: error, warnings and informational messages;

The ``debug`` (debug messages) and ``trace`` (very verbose debug) are also
available, but mainly intended for developers.

The ``append`` key is a boolean value (``true`` or ``false``) only used for
files, and controls whether to overwrite the file or append new messages at the
end of an existing file. It defaults to ``false``, meaning that the file is
overwritten by every simulation run.

.. code::

    [log]
    # Use multiple targets
    targets = [
        # Print warnings to the standard output stream
        {target = "<stdout>", level = "warning"},
        # Save all messages to the 'lumol.log' file
        {target = "lumol.log"},
        # Save debug messages to 'debug.log', keeping the file across simulation
        # runs.
        {target = "debug.log", level = "debug", append = true},
    ]
