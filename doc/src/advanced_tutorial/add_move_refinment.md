# Refining the implementation

This section will talk about:

- Restrict `dr` to sane values (for example a box length) in the `setup` function.
- Don't store the system but only the molecule that is subject to change.
- Use the cache to compute only the incremental change to energy.
- Implement an update of the displacement and maybe talk about Markov Chains (microscopic reversibility) here?