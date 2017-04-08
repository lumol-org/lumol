# Logging configuration

Lumol sends various logging messages while running a simulation. Some of them
are informational messages (`charge was set to 1.2 for 132 particles`), others
are warnings or error message (`infinite energy!`) and some are for debug
purposes.

By default, lumol prints all the informational, warning and error messages to
the standard terminal output. This allows to run the code and redirect the
output to a specific file with the usual UNIX way: `lumol input.toml >
simulation.log`.

Lumol also offers finer-grain configuration for logging output, for example only
printing errors and warnings, and redirecting everything else to a file. This
configuration happens in the `[log]` section of the input file. This section can
contain either a single output `target`, or multiple `targets`.

```toml
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
```

In the multiple targets case, the `targets` key must be an array of tables,
containing individual target configuration with the same syntax as a single
target.

The only required key is the `target` key, identifying where to write the
messages. It is interpreted as the path to a file, expect for the two special
values of `<stdout>` and `<stderr>` that are used to write messages to the
standard terminal output stream or error stream respectively.

```toml
# Write all messages to the standard error stream
[log]
target = "<stderr>"
```

Optional keys are `level` and `append`. The `level` key controls which messages
are sent to this output, and default to `info`. Available levels are the
following:

- `error`: only error messages;
- `warning`: error and warning messages;
- `info`: error, warnings and informational messages;

The `debug` (debug messages) and `trace` (very verbose debug) are also
available, but mainly intended for developers.

The `append` key is a Boolean value (`true` or `false`) only used for files, and
controlling whether to overwrite the file or append new messages at the end of
the file. It defaults to `false`, meaning that the file is overwritten by every
simulation run.

```toml
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

```
